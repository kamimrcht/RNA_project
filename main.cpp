#include <fstream>
#include <cstring>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>

using namespace std;

struct window{
	int index;
	int read;

    bool operator==(const window& a) const
	{
		return (index == a.index && read == a.read);
	}
};


uint64_t transformWindowToHash(window win){
	hash<int> winHash;
	return winHash(win.read + win.index);
}


namespace std { template <> struct hash<window> {
	typedef window argument_type;
	typedef uint64_t result_type; uint64_t operator()(window key) const { return transformWindowToHash(key); } }; }


struct compareWindow{
    bool operator()(const window& win1, const window& win2){
        return win1.index <win2.index;
    }
};



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


void setKmersToWindows(int indexRead, int indexWindow, string kmer, const unordered_map<string, int>& solidKmers, unordered_map <string, vector <window>>& kmerToWindows){
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


void getKmersinWindowsFromReads(int k, int w, const vector <string>& readSet, const unordered_map<string, int>& solidKmers, unordered_map <string, vector<window>>& kmerToWindows){
	for (uint indexRead(0); indexRead < readSet.size(); ++ indexRead){
		string readSequence(readSet[indexRead]);
		if (not readSequence.empty()) {
			int position(0);
			int posiForKmer;
			string kmer;
			int indexWindow(0);
			while (position != -1){
				if (position + w + k - 1 < (int)readSequence.size()){
					posiForKmer = position;
					position += w;
				} else {
					posiForKmer = readSequence.size() - w - k + 1;
					position = -1;
				}
				for (int iter(posiForKmer); iter < posiForKmer + w; ++iter){
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


void getKmersinFromReadsInMap(int k, const vector <string>& readSet, unordered_map <string, int>& kmersFromFile){
	for (uint readIndex(0); readIndex < readSet.size(); ++ readIndex){
		string readSequence(readSet[readIndex]);
		if (not readSequence.empty()){
			int position(0);
			string kmer;
			
			while (position + k - 1 < (int)readSequence.size()){
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


void getSolidKmers(const unordered_map<string, int>& kmersFromFile, unordered_map<string, int>& solidKmers){
	 for ( unordered_map<string, int>::const_iterator iter = kmersFromFile.begin(); iter != kmersFromFile.end(); ++iter ){
		if (iter->second > 1){
			solidKmers[iter->first] = iter->second;
		} 
	 }
}


void compareReadWindows(int k, int w, const vector<string>& readSet,const unordered_map<string, int> solidKmers, const unordered_map<string, vector<window>> kmerToWindows){
	for (uint indexRead(0); indexRead < readSet.size(); ++ indexRead){
		string readSequence(readSet[indexRead]);
		if (not readSequence.empty()){
			int position(0), posiForKmer(0);
			uint indexWindow(0);
			while (position != -1){
				uint nbKmers(0);
				unordered_map<window, double> similarity;
				if (position + w + k - 1 <(int) readSequence.size()){
					posiForKmer = position;
					position += w;
				} else {
					posiForKmer = readSequence.size() - w - k + 1;
					position = -1;
				}
				for (int iter(posiForKmer); iter < posiForKmer + w; ++iter){
					string kmer = getKmer(readSequence, iter, k);
					if (solidKmers.count(kmer)){
						++ nbKmers;
						for (auto iter = kmerToWindows.begin(); iter != kmerToWindows.end(); ++iter ){
							for (uint i(0); i<iter->second.size(); ++i){
								if (iter->second[i].read != (int)indexRead){
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


int main(int argc, char ** argv){
	if (argc < 4){
		cout << "command line: ./rnaLR reads.fasta k w" << endl;
	} else {
		string fileName = argv[1];
		uint k = stoi(argv[2]);
		uint w = stoi(argv[3]);
		ifstream readFile(fileName);
		vector <string> readSet;
		string sequence;
		while (not readFile.eof()){
            getline(readFile, sequence);
			getline(readFile, sequence);
			readSet.push_back(sequence);
		}
		unordered_map <string, int> kmersFromFile; //  TODO: destroy
		unordered_map <string, int> solidKmers;
		unordered_map <string, vector<window>> kmerToWindows;
		getKmersinFromReadsInMap(k, readSet, kmersFromFile);
		getSolidKmers(kmersFromFile, solidKmers);
		getKmersinWindowsFromReads(k, w, readSet, solidKmers, kmerToWindows);
		compareReadWindows(k, w, readSet, solidKmers, kmerToWindows);
	}
	return 0;
}

