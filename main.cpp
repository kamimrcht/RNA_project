#include <fstream>
#include <cstring>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>

using namespace std;

//~ string revComp(const string& seq){
	//~ string revCompSeq = "";
	//~ int pos = seq.size()-1;
	//~ char nt;
	//~ do{
		//~ nt = seq[pos];
		//~ switch (nt) {
			//~ case 'A':
				//~ revCompSeq += 'T';
				//~ break;
			//~ case 'T':
				//~ revCompSeq += 'A';
				//~ break;
			//~ case 'C':
				//~ revCompSeq += 'G';
				//~ break;
			//~ case 'G':
				//~ revCompSeq += 'C';
				//~ break;
		//~ }
		//~ --pos;
	//~ } while (pos>=0);
	//~ return revCompSeq;
//~ }


//~ string getCanonical(const string& seq){
	//~ string revCompSeq = revComp(seq);
	//~ return min(seq, revCompSeq); 
//~ }


//~ string getKmer(const string& sequence, int posi, int k){
	//~ return getCanonical(sequence.substr(posi, k));
//~ }


//~ void getKmersFromRead(const string& readSeq, int k, unordered_map<string, int>& kmersFromFile){
	//~ for (uint posi(0); posi < readSeq.size()-k+1; ++posi){
		//~ string kmer = getKmer(readSeq, posi, k);
		//~ if (not kmersFromFile.count(kmer)){
			//~ kmersFromFile[kmer] = 1;
		//~ } else {
			//~ kmersFromFile[kmer] += 1;
		//~ }
	//~ }
//~ }


//~ void getSolidKmers(const unordered_map<string, int>& kmersFromFile, unordered_map<string, int>& solidKmers){
	 //~ for ( unordered_map<string, int>::const_iterator iter = kmersFromFile.begin(); iter != kmersFromFile.end(); ++iter ){
		//~ if (iter->second >1){
			//~ solidKmers[iter->first] = iter->second;
		//~ } 
	 //~ }
//~ }


//~ vector<int> removeDuplicates(const vector<int>& vect){
	//~ vector<int> vectResult;
	//~ for (uint i(0); i< vect.size(); ++i){
		//~ if (i == vect.size()-1 or vect[i]!=vect[i+1]){
			//~ vectResult.push_back(vect[i]);
		//~ }
	//~ }
	//~ return vectResult;
//~ }


//~ void putKmersToWindowsTarget(const unordered_map <string, int>& solidKmers, unordered_map <string, vector<int>>& kmerToWindows, const string& target, int indexWindow, uint posiOnTarget, int w, int k){
	//~ for (uint posInW(0); posInW < posiOnTarget+w; ++posInW){
		//~ string kmer = getKmer(target, posInW, k);
		//~ if (solidKmers.count(kmer)){
			//~ if (not kmerToWindows.count(kmer)){
				//~ vector<int> vecIndexes;
				//~ vecIndexes.push_back(indexWindow);
				//~ kmerToWindows[kmer] = vecIndexes;
			//~ } else {
				//~ kmerToWindows[kmer].push_back(indexWindow);
				//~ sort(kmerToWindows[kmer].begin(), kmerToWindows[kmer].end());
				//~ vector <int> vect = removeDuplicates(kmerToWindows[kmer]);
				//~ kmerToWindows[kmer] = vect;
			//~ }
		//~ } else {
			//~ cout << kmer << endl;
		//~ }
	//~ }
//~ }

//~ void getSimilarityWindowQuery(int posiOnQuery, string readSeq, int k, int w, unordered_map <string, int> solidKmers, unordered_map <string, vector<int>> kmerToWindowsTarget, unordered_map<int,double>& similarity){
	//~ int nbKmersinWindowQuery(0);
	//~ for (int posInW(0); posInW < posiOnQuery+w; ++posInW){
		//~ string kmer = getKmer(readSeq, posInW, k);
		//~ if (solidKmers.count(kmer)){
			//~ ++ nbKmersinWindowQuery;
			//~ if (kmerToWindowsTarget.count(kmer)){ //  kmer is in target
				//~ for (uint i(0); i<kmerToWindowsTarget[kmer].size();++i){
					//~ int winT(kmerToWindowsTarget[kmer][i]);
					//~ if (similarity.count(winT)==1){
						//~ ++ similarity[winT];
					//~ } else {
						//~ similarity[winT] = 1;
					//~ }
				//~ }
			//~ }
		//~ }
	//~ }
	//~ for (auto iter(similarity.begin()); iter != similarity.end(); ++iter ){
		//~ if (iter->second != 0){
			//~ iter->second /= nbKmersinWindowQuery;
		//~ }
	//~ }
//~ }


//~ int main(int argc, char ** argv){
	//~ if (argc < 4){
		//~ cout << "command line: ./rnaLR reads.fasta k w" << endl;
	//~ } else {
		//~ int k = stoi(argv[2]);
		//~ int w = stoi(argv[3]);
		//~ string target = "AATCGATTCTT"; 
		//~ vector<string> query = {"AATCGATTCTTGTGGGCCCTGAGATCGATTCTT", revComp(target)};
		//~ unordered_map <string, int> kmersFromFile;
		//~ unordered_map <string, int> solidKmers;
		//~ unordered_map <string, vector<int>> kmerToWindows;
		//~ getKmersFromRead(target, k, kmersFromFile);
		//~ for (uint s(0); s < query.size(); ++s){
			//~ getKmersFromRead(query[s], k, kmersFromFile);
		//~ }
		//~ getSolidKmers(kmersFromFile, solidKmers);
		//~ uint posi(0);
		//~ int indexWindow(0);
		//~ do{
			//~ putKmersToWindowsTarget(solidKmers, kmerToWindows, target, indexWindow, posi, w, k);
			//~ ++ indexWindow;
			//~ posi += w;
		//~ } while (posi < target.size()-w-k+2);
		
		//~ for (uint reads(0); reads < query.size(); ++reads){
			//~ cout << "read n" <<  reads << endl;
			//~ uint posiOnQuery(0);
			//~ int indexWindowQuery(0);
			//~ do{
			//~ unordered_map<int,double> similarity;  //key: window on target/value: score
			//~ getSimilarityWindowQuery(posiOnQuery, query[reads], k, w, solidKmers, kmerToWindows, similarity);
			//~ for (auto iter = similarity.begin(); iter != similarity.end(); ++iter ){
				//~ if (iter->second >= 0.5){
					//~ cout << iter->first << " " << indexWindowQuery << " " << iter->second << endl;
				//~ }
			//~ }
			//~ ++ indexWindowQuery;
			//~ posiOnQuery += w;
			//~ } while (posiOnQuery < query[reads].size()-w-k+2);
		//~ }
		//~ return 0;
	//~ }
//~ }


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
				for (auto iter = similarity.begin(); iter != similarity.end(); ++iter ){
					iter->second /= nbKmers;
					if (iter->second >= 0.7){
						cout << "reads:" << iter->first.read << " " << indexRead << " windows:" << iter->first.index << " " << indexWindow <<  " score:" << iter->second <<endl;
					}
				}
				++indexWindow;
			}
		}
		break;
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
        //~ ofstream out("out.fa");
		//~ vector <string> readSet ({"AATCGATTCTT", "AATCGATTCTT", "AATCGATTCTT"});
		//~ vector <string> readSet ({"AATCGATTCTT", revComp("AATCGATTCTT")});
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

