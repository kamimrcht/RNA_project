#include <fstream>
#include <cstring>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include "compareReadsByWindows.h"
#include "utils.h"
#include "detectJumps.h"

using namespace std;

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
		//~ string readSequence(getCanonical(readSet[readIndex]));
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


void getStartKmerEnum(int& position, uint& startForKmer, uint sizeSeq, uint k, uint w){
	if (position + w + k - 1 < sizeSeq){
		startForKmer = position;
		position += w;
	} else {
		startForKmer = sizeSeq - w - k + 1;
		position = -1;
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
				getStartKmerEnum(position, posiForKmer, readSequence.size(), k, w);
				for (uint iter(posiForKmer); iter < posiForKmer + w; ++iter){
					kmer = getKmer(readSequence, iter, k);
					setKmersToWindows(indexRead, indexWindow, kmer, solidKmers, kmerToWindows);
				}
				++indexWindow;
			}
		}
	}
}

//   computes similirity for each window against all other windows from other reads, and return the number of solid kmers in the window to compute the final similarity
uint getSimilarity(uint posiForKmer, uint& nbKmers, uint indexReadTarget, const string& readSequence, const unordered_map<string, uint> solidKmers, const unordered_map<string, vector<window>> kmerToWindows, unordered_map<window, double>& similarity, uint k, uint w){
	for (uint iter(posiForKmer); iter < posiForKmer + w; ++iter){ // for each kmer in target window
		string kmer = getKmer(readSequence, iter, k);
		if (solidKmers.count(kmer)){
			++ nbKmers;
			for (auto iter = kmerToWindows.begin(); iter != kmerToWindows.end(); ++iter ){  //looking for a corresponding kmer in indexed windows
				for (uint i(0); i<iter->second.size(); ++i){
					if (iter->second[i].read != indexReadTarget){ // if we are not in the same read
						if (similarity.count(iter->second[i])){
							if (kmer == iter->first){
								++similarity[iter->second[i]];
							}
						} else {
							if (kmer == iter->first){
								similarity[iter->second[i]] = 1;
							} else {
								similarity[iter->second[i]] = 0;
							}
						}
					}
				}
			}
		}
	}
	return nbKmers;
}

void getSimilarWindowsPairs(uint readSetSize, uint indexWindowTarget, unordered_map<window, double>& similarity, uint nbKmers, unordered_map <uint, vector <pairOfIndexWindow>>& pairsOfaRead, uint w){
	for (uint rIndex(0); rIndex<readSetSize; ++rIndex){  // to order the output in cout
		for (auto iter = similarity.begin(); iter != similarity.end(); ++iter ){
			if (iter->first.read == rIndex){  // to order the output in cout
				iter->second /= nbKmers;
				if (iter->second > 0.7 and nbKmers > w/2){
					pairOfIndexWindow pair({indexWindowTarget, iter->first.index});
					if (pairsOfaRead.count(rIndex)){
						pairsOfaRead[rIndex].push_back(pair);
					} else {
						vector <pairOfIndexWindow> vec;
						vec.push_back(pair);
						pairsOfaRead[rIndex] = vec;
					}
				}
				cout << indexWindowTarget <<  " " << iter->first.index <<  " " << iter->second << endl; // for heatmaps
			}
		}
	}
}

void compareAndCorrectReadWindows(uint k, uint w, vector<string>& readSet,const unordered_map<string, uint> solidKmers, const unordered_map<string, vector<window>> kmerToWindows){
	for (uint indexReadTarget(0); indexReadTarget < readSet.size(); ++ indexReadTarget){
		unordered_map <regionInRead, vector<regionInRead>> correspondance;
		unordered_map <uint, vector <pairOfIndexWindow>> pairsOfaRead;
		string readSequence(readSet[indexReadTarget]);
		if (not readSequence.empty()){
			int position(0);
			uint posiForKmer(0);
			uint indexWindowTarget(0);
			while (position != -1){ // for each window in target read
				uint nbKmers(0);
				unordered_map<window, double> similarity;
				getStartKmerEnum(position, posiForKmer, readSequence.size(), k, w);
				nbKmers = getSimilarity(posiForKmer, nbKmers, indexReadTarget, readSequence, solidKmers, kmerToWindows, similarity, k, w);
				getSimilarWindowsPairs(readSet.size(), indexWindowTarget, similarity, nbKmers, pairsOfaRead, w);
				++indexWindowTarget;
			}
			for (auto iter = pairsOfaRead.begin(); iter != pairsOfaRead.end(); ++iter){
				if (not iter->second.empty()){
					sort(iter->second.begin(), iter->second.end(), comparePairOfIndexWindow());
					detectJumps(iter->second, indexReadTarget, iter->first, correspondance);
				}
			}
			//~ cout << "*** correction phase" << endl;
			//~ consensusBetweenRegions(correspondance, readSet, w, k);
			//~ for (auto iter=correspondance.begin(); iter!=correspondance.end(); ++iter){
				//~ for (uint i(0); i < iter->second.size(); ++i){
					//~ cout << iter->first.firstWindow  <<  " " << iter->first.lastWindow << " : " << iter->second[i].firstWindow << " " << iter->second[i].lastWindow << endl;
				//~ }
			//~ }
		}
		break; // REMOVE IF all v all
	}
}
