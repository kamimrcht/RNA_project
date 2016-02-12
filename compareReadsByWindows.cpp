#include <fstream>
#include <cstring>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include "compareReadsByWindows.h"

using namespace std;

//~ struct window{
	//~ uint index;
	//~ uint read;

    //~ bool operator==(const window& a) const
	//~ {
		//~ return (index == a.index && read == a.read);
	//~ }
//~ };


//~ uint64_t transformWindowToHash(window win){
	//~ hash<int> winHash;
	//~ return winHash(win.read + win.index);
//~ }


//~ namespace std { template <> struct hash<window> {
	//~ typedef window argument_type;
	//~ typedef uint64_t result_type; uint64_t operator()(window key) const { return transformWindowToHash(key); } }; }


//~ struct compareWindow{
    //~ bool operator()(const window& win1, const window& win2){
        //~ return win1.index <win2.index;
    //~ }
//~ };



uint64_t transformWindowToHash(window win){
	hash<int> winHash;
	return winHash(win.read + win.index);
}


vector<window> removeDuplicates(const vector<window>& vect){
	vector<window> vectResult;
	for (uint i(0); i< vect.size(); ++i){
		if (i == vect.size()-1 or not (vect[i].index == vect[i+1].index and vect[i].read == vect[i+1].read)){
			vectResult.push_back(vect[i]);
		}
	}
	return vectResult;
}


string revComp(const string& seq){
	string revCompSeq = "";
	int pos = seq.size()-1;
	char nt;
	do{
		nt = seq[pos];
		switch (nt) {
			case 'A':
				revCompSeq += 'T';
				break;
			case 'T':
				revCompSeq += 'A';
				break;
			case 'C':
				revCompSeq += 'G';
				break;
			case 'G':
				revCompSeq += 'C';
				break;
		}
		--pos;
	} while (pos>=0);
	return revCompSeq;
}


string getCanonical(const string& seq){
	string revCompSeq = revComp(seq);
	return min(seq, revCompSeq); 
}


string getKmer(const string& sequence, int posi, int k){
	return getCanonical(sequence.substr(posi, k));
}


void setKmersToWindows(uint indexRead, uint indexWindow, string kmer, const unordered_map<string, uint>& solidKmers, unordered_map <string, vector <window>>& kmerToWindows){
	if (solidKmers.count(kmer)){
		window win({indexWindow, indexRead});
		if (kmerToWindows.count(kmer)){
			kmerToWindows[kmer].push_back(win);
			sort(kmerToWindows[kmer].begin(), kmerToWindows[kmer].end(), compareWindow());
			kmerToWindows[kmer] = removeDuplicates(kmerToWindows[kmer]);
		} else {
			vector <window> vectWin({win});
			kmerToWindows[kmer] = vectWin;
		}
	}
}



void getKmersinFromReadsInMap(uint k, const vector <string>& readSet, unordered_map <string, uint>& kmersFromFile){
	for (uint readIndex(0); readIndex < readSet.size(); ++ readIndex){
		string readSequence(readSet[readIndex]);
		if (not readSequence.empty()){
			uint position(0);
			string kmer;
			
			while (position + k - 1 < readSequence.size()){
				kmer = getKmer(readSequence, position, k);
				if (kmersFromFile.count(kmer)){
					++kmersFromFile[kmer];
				} else {
					kmersFromFile[kmer] = 1;
				}
				++position;
			}
		}
	}
}


void getSolidKmers(const unordered_map<string, uint>& kmersFromFile, unordered_map<string, uint>& solidKmers){
	 for (auto iter = kmersFromFile.begin(); iter != kmersFromFile.end(); ++iter ){
		if (iter->second > 1){
			solidKmers[iter->first] = iter->second;
		} 
	 }
}


void getKmersinWindowsFromReads(uint k, uint w, const vector <string>& readSet, const unordered_map<string, uint>& solidKmers, unordered_map <string, vector<window>>& kmerToWindows){
	for (uint indexRead(0); indexRead < readSet.size(); ++ indexRead){
		string readSequence(readSet[indexRead]);
		if (not readSequence.empty()) {
			int position(0);
			uint posiForKmer;
			string kmer;
			uint indexWindow(0);
			while (position != -1){
				if (position + w + k - 1 < readSequence.size()){
					posiForKmer = position;
					position += w;
				} else {
					posiForKmer = readSequence.size() - w - k + 1;
					position = -1;
				}
				for (uint iter(posiForKmer); iter < posiForKmer + w; ++iter){
					kmer = getKmer(readSequence, iter, k);
					setKmersToWindows(indexRead, indexWindow, kmer, solidKmers, kmerToWindows);
				}
				++indexWindow;
			}
			//~ for (auto iter = kmerToWindows.begin(); iter != kmerToWindows.end(); ++iter ){
				//~ cout << iter->first << endl;
				//~ for (uint i(0); i <iter->second.size();++i){
					//~ cout << "index:"<<iter->second[i].index <<  " read:" <<iter->second[i].read << endl;
				//~ }
			//~ }
		}
	}
}



void compareReadWindows(uint k, uint w, const vector<string>& readSet,const unordered_map<string, uint> solidKmers, const unordered_map<string, vector<window>> kmerToWindows){
	for (uint indexRead(0); indexRead < readSet.size(); ++ indexRead){
		string readSequence(readSet[indexRead]);
		if (not readSequence.empty()){
			int position(0);
			uint posiForKmer(0);
			uint indexWindow(0);
			while (position != -1){
				uint nbKmers(0);
				unordered_map<window, double> similarity;
				if (position + w + k - 1 < readSequence.size()){
					posiForKmer = position;
					position += w;
				} else {
					posiForKmer = readSequence.size() - w - k + 1;
					position = -1;
				}
				for (uint iter(posiForKmer); iter < posiForKmer + w; ++iter){
					string kmer = getKmer(readSequence, iter, k);
					if (solidKmers.count(kmer)){
						++ nbKmers;
						for (auto iter = kmerToWindows.begin(); iter != kmerToWindows.end(); ++iter ){
							for (uint i(0); i<iter->second.size(); ++i){
								if (iter->second[i].read != indexRead){
									if (kmer == iter->first){
										if (similarity.count(iter->second[i])){
											++similarity[iter->second[i]];
										} else {
											similarity[iter->second[i]] = 1;
										}
									}
								}
							}
						}
					}
				}
				for (uint rIndex(0); rIndex<readSet.size();++rIndex){  // to order the output in cout
					//~ cout << rIndex << endl;
					for (auto iter = similarity.begin(); iter != similarity.end(); ++iter ){
						if (iter->first.read == rIndex){  // to order the output in cout
							iter->second /= nbKmers;
							if (iter->second >= 0.7){
								cout << "reads:" << iter->first.read << " " << indexRead << " windows:" << iter->first.index << " " << indexWindow <<  " score:" << iter->second <<endl;
							}
						}
					}
				}
				++indexWindow;
			}
		}
		break; // REMOVE IF all v all
	}
}
