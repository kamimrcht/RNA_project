#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <unordered_set>
#include <stdio.h>
#include <string>

using namespace std;

uint THRESHOLD1(1);
uint THRESHOLD2(3);
uint THRESHOLD3(1);
uint D(2);
uint SIZE_INSERTION(2);

void fillSignatures(vector<string>& vectSequences, uint k, unordered_map<string, vector<int>>& kmerToSignature){
	string kmer;
	vector<string> toDelete;
	unordered_map<string, vector<uint>> kmerToIndex;
	/* keep only kmers that appear once in a read. Link kmer -> read index in a map */
	for (uint r(0); r < vectSequences.size(); ++r){
		for (uint p(0); p < vectSequences[r].size() - k + 1; ++p){
			kmer =  vectSequences[r].substr(p, k);
			if (kmerToIndex.count(kmer)){
				if (kmerToIndex[kmer].back() == r){
					toDelete.push_back(kmer);
				} else {
					kmerToIndex[kmer].push_back(r);
				}
			} else {
				kmerToIndex[kmer].push_back(r);
			}
		}
	}
	for (uint k(0); k < toDelete.size(); ++k){
		kmerToIndex.erase(toDelete[k]);             // erasing by key
	}
	/* fill signature for each kmer of the prevous map */
	for (auto kmer(kmerToIndex.begin()); kmer != kmerToIndex.end(); ++kmer){
		vector<int> signature(vectSequences.size(), 0);
		for (uint r(0);  r < vectSequences.size(); ++r){
			for (uint rr(0); rr < kmer->second.size(); ++rr){
				if (r == kmer->second[rr]){
					signature[r] = 1;
					break;
				}
			}
		}
		kmerToSignature.insert({kmer->first, signature});
	}
}


void dispatchErrors(vector<string>& vectSequences, uint k, unordered_map<string, vector<int>>& kmerToSignature, unordered_map<string, vector<pair<uint, uint>>>& kmerToReadPosi, unordered_map<uint, vector<uint>>& readToErrorPosition){
	string kmer;
	uint precP(0);
	for (uint r(0); r < vectSequences.size(); ++r){
		for (uint p(0); p < vectSequences[r].size() - k + 1; ++p){
			kmer =  vectSequences[r].substr(p, k);
			auto got = kmerToSignature.find(kmer);
			if (got != kmerToSignature.end()){
				uint sumSignature(0);
				for (uint n : got->second){
					sumSignature += n;
				}
				if (sumSignature > THRESHOLD1) { // we consider the kmer is errorless
					pair<uint, uint> rp({r, p});
					auto got2 = kmerToReadPosi.find(kmer);
					if (got2 != kmerToReadPosi.end()){
						got2->second.push_back(rp);
					} else {
						kmerToReadPosi.insert({kmer, {rp}});
					}
				} else { // kmer with error
					readToErrorPosition[r].push_back(p);
				}
			}
		}
	}
}

vector<string> generateAltKmer(string& kmer ){
	vector<string> ntInt({"A", "C", "G", "T"});
	vector<string> result;
	char ntc(kmer.back());
	string nt; nt.push_back(ntc);
	for (uint t(0); t < kmer.size(); ++t){
		for (string i : ntInt){
			if (i != nt){
				result.push_back(kmer.substr(0, t) + i + kmer.substr(t + 1));
				//~ result.push_back(kmer.substr(0, kmer.size()-1) + i);
			}
		}
	}
	return result;
}


string makeConsensus(vector<string>& seqVec){
	string nt;
	vector<string> ntInt({"A", "C", "G", "T"});
	string nuc;
	string result;
	for (uint i(0); i < seqVec[0].size(); ++i){
		vector<uint> score(4, 0);
		for (string s : seqVec){
			nt = s.substr(i, 1);
			if (nt == "A"){
				++score[0];
			} else if (nt == "C"){
				++score[1];
			} else if (nt == "G"){
				++score[2];
			} else {
				++score[3];
			}
		}
		uint max(0);
		for (uint n(0); n < score.size(); ++n){
			if (score[n] > max){
				max = score[n];
				nuc = ntInt[n];
			}
		}
		result += nuc;
	}
	return result;
}


string returnSeqBetween2Kmers(uint k, uint read, uint start, string& nextKmer, vector<string>& vectSequences){
	string result(vectSequences[read].substr(start, k));
	uint incr(k + 1);
	uint size(vectSequences[read].size() - start);
	while (result.substr(result.size()-k) != nextKmer and result.size() < size){
		result = vectSequences[read].substr(start, incr);
		++incr;
	}
	if (result.substr(result.size()-k) == nextKmer){
		return result.substr(0, result.size()-1);
	} else {
		return "";
	}
}



void followingKmers(vector<string>& vectSequences, uint k, unordered_map<string, vector<pair<uint, uint>>>& kmerToReadPosi, unordered_set<string>& branchingKmersLeft, unordered_set<string>& branchingKmersRight){
	string kmer, prev, next;
	unordered_map<string, vector<string>> kmerLeft, kmerRight;
	for (uint r(0); r < vectSequences.size(); ++r){
		for (uint p(0); p < vectSequences[r].size() - k + 1; ++p){
			kmer = vectSequences[r].substr(p, k);
			if (kmerToReadPosi.count(kmer)){ // only on errorless kmers
				if (p == 0){
					prev = "";
					next = vectSequences[r].substr(p + 1, k);
				} else if (p == vectSequences[r].size() - k){
					prev = vectSequences[r].substr(p - 1, k);
					next = "";	
				} else {
					prev = vectSequences[r].substr(p - 1, k);
					next = vectSequences[r].substr(p + 1, k);
				}
				auto gotn(kmerToReadPosi.find(next));
				auto gotp(kmerToReadPosi.find(prev));
				auto got(kmerLeft.find(kmer));
				auto got2(kmerRight.find(kmer));
				if (gotn != kmerToReadPosi.end()){
					if (got2 != kmerRight.end()){
						got2->second.push_back(next);
					} else {
						vector<string> v = {next};
						kmerRight[kmer] = v;
					}
				}
				if (gotp != kmerToReadPosi.end()){
					if (got != kmerLeft.end()){
						got->second.push_back(prev);
					} else {
						vector<string> v = {prev};
						kmerLeft[kmer] = v;
					}
				}
			}
		}
	}
	uint s;
	for (auto p(kmerLeft.begin()); p != kmerLeft.end(); ++p){
		s = 0;
		vector<string> newVec;
		sort(p->second.begin(), p->second.end());
		bool same(true);
		while (s < p->second.size() - 1){
			if (p->second[s] != p->second[s + 1]){
				newVec.push_back(p->second[s]);
				same = false;
				
			} else {
				if (same == false){
					newVec.push_back(p->second[s]);
				}
				same = true;
			}
			++s;
		}
		if (not same and s == p->second.size() - 1){ // last element
			newVec.push_back(p->second[s]);
		}
		if (newVec.size() > 1){
			branchingKmersLeft.insert(p->first);
		}
	}
	
	for (auto n(kmerRight.begin()); n != kmerRight.end(); ++n){
		s = 0;
		vector<string> newVec;
		sort(n->second.begin(), n->second.end());
		bool same(true);
		while (s < n->second.size() - 1){
			if (n->second[s] != n->second[s + 1]){
				newVec.push_back(n->second[s]);
				same = false;
			} else {
				if (same == false){
					newVec.push_back(n->second[s]);
				}
				same = true;
			}
			++s;
		}
		if(not same and s == n->second.size() - 1){ // last element
			newVec.push_back(n->second[s]);
		}
		if (newVec.size() > 1){
			branchingKmersRight.insert(n->first);
		}
	}
}

// gets the most frequent sequence length
uint sizeToKeep(vector<uint>& vecSizes){ 
	uint maxVal, maxFreq(0), freq(0);
	sort(vecSizes.begin(), vecSizes.end());
	for (uint s(1); s < vecSizes.size(); ++s){
		if (vecSizes[s] == vecSizes[s-1]){
			++freq;
			if (s == vecSizes.size() - 1 and maxFreq == 0){
				maxFreq = vecSizes.size();
				maxVal = vecSizes[s];
			}
		} else {
			if (freq > maxFreq){
				maxFreq = freq;
				maxVal = vecSizes[s-1];
			}
		}
	}
	if (maxFreq > 1){
		return maxVal;
	} else {
		return 0;
	}
}


vector<int> sumOfProfiles(vector<int>& p1, vector<int>& p2){
	vector<int> result;
	for (uint i(0); i < p1.size();  ++i){
		result.push_back(p1[i] + p2[i]);
	}
	return result;
}


unordered_set<int> getZeros(vector<int>& sum){
	unordered_set<int> result;
	for (uint i(0); i < sum.size(); ++i){
		if (sum[i] == 0){
			result.insert(i);
			//~ cout << "zero at position " << i << endl;
		}
	}
	return result;
}


unordered_set<int> zerosOfInterval(vector<string>& vectSequences, uint read, vector<uint>& positionsVec, uint k, unordered_map<string, vector<pair<uint, uint>>>& kmerToReadPosi,  unordered_map<string, vector<int>>& kmerToSignature){
	//~ uint posi(position1 + 1), nbKmers(0);
	uint posi(1), nbKmers(0);
	unordered_set<int> result;
	string kmer;
	vector<int> profile, sum;
	//~ for (uint i(0); i < positionsVec.size(); ++i){
		//~ cout << "gg "<<i << endl;
	//~ }
	//~ cout << positionsVec.size() << endl;
	//~ cout << "summed " << endl;
	kmer = vectSequences[read].substr(positionsVec[0], k);
	sum =  kmerToSignature[kmer];
	while (posi < positionsVec.size()){ //TODO : maybe less than k -1 ?
		kmer = vectSequences[read].substr(positionsVec[posi], k);
		if (kmerToReadPosi.count(kmer)){
			profile = kmerToSignature[kmer];
			//~ cout << "SUM" << endl; 
			//~ for (uint i : sum){
				//~ cout << i;
			//~ }
			
			//~ cout << endl;
			sum = sumOfProfiles(profile, sum);
			//~ cout << "SUM" << endl; 
			//~ for (uint i : sum){
				//~ cout << i;
			//~ }
			//~ cout << endl;
			++nbKmers;
		}
		++posi;
	}
	if (sum.size() > 2){ // TODO: bad threshold
		result = getZeros(sum);
	}
	
	return result;
}

vector<uint> positionsRight(uint start, uint k){
	vector<uint> result;
	for (uint i(start); i < start + k - 1; ++i){
		result.push_back(i);
	}
	return result;
}

vector<uint> positionsLeft(uint start, uint k){
	vector<uint> result;
	for (uint i(start); i > start - k + 1; --i){
		result.push_back(i);
	}
	return result;
}



void propagateSignature(vector<uint>& signatures, uint position1, uint position2){
	uint lastSign;
	if (position1 == 0){
		lastSign = 1;
	} else {
		lastSign = signatures[position1 - 1] + 1;
	}
	for (uint i(position1); i < position2; ++i){
		signatures[i] = lastSign;
	}
}


void computeExons(vector<string>& vectSequences, uint k,  unordered_map<string, vector<pair<uint, uint>>>& kmerToReadPosi,  unordered_map<string, vector<int>>& kmerToSignature, unordered_set<string>& kmerLeft, unordered_set<string>& kmerRight){
	uint posi, posiPrec, posiNext, posiPrev, lastJunction, lastSignature;
	string kmer1, kmer2, currentKmer, lastTrueKmer("");
	bool found, foundk2, confirmed, confirmed2;
	//~ vector<int> profile1, profile2;
	unordered_set<int> zerosInPrec, zerosInNext;
	for (uint r(0); r < vectSequences.size(); ++r){
		cout << r << endl;
		posi = k - 1; lastJunction = 0;  lastSignature = 0; confirmed = false;
		vector<uint> signatures(vectSequences[r].size(), 0);
		while (posi < vectSequences[r].size() - k - 1){
			vector<int> sum, sum2;
			found = false;
			//// *** CHUNK 1: try and find a first branching k-mer k1 ////
			kmer1 = vectSequences[r].substr(posi, k);
			auto gotNE(kmerToReadPosi.find(kmer1));
			if (gotNE != kmerToReadPosi.end()){ // errorless kmer1
				lastTrueKmer = kmer1;
			} else { //kmer1 with error : we try to mutate it to see if we have registered  an branching kmer that is nearly similar
				vector<string> altk1V = generateAltKmer(kmer1);
				for (string altkmer1 : altk1V){
					if (kmerRight.count(altkmer1)){ // we found a nearly similar branching kmers that opens a bubble
						kmer1 = altkmer1;
						break;
					} else if (kmerLeft.count(altkmer1)) { // we found a nearly similar branching kmers that closes a bubble
						kmer1 = altkmer1;
						break;
					} else {
						kmer1 = ""; // we found nothing we will just continue
					}
				}
			}
			//// *** CHUNK 2: if a branching k-mer k1 is found, try to find a know neighbour k2 ////
			if (not kmer1.empty()){ // else we were unable to find a kmer1 without error, we continue
				auto got1(kmerRight.find(kmer1));
				//// ** Beginning of alternative exons (opening a bubble) ////
				if (got1 != kmerRight.end()){
					// try and get an errorless kmer2
					posiNext = posi + 1;
					kmer2 = vectSequences[r].substr(posiNext, k);
					auto gotNE2(kmerToReadPosi.find(kmer2));
					bool foundk2 = false;
					if (gotNE2 != kmerToReadPosi.end()){
						foundk2	= true;
						cout << "k2 opening ok" << endl;
					}
					if (not foundk2){
						while (posiNext < posi + 1 +  k + SIZE_INSERTION){ // TODO : is it a good interval for the research
							++posiNext;
							kmer2 = vectSequences[r].substr(posiNext, k);
							auto gotNE2 = kmerToReadPosi.find(kmer2);
							if ( gotNE2 != kmerToReadPosi.end() ){
								foundk2 = true;
								cout << "new k2 opening ok" << endl;
								break;
							}
						}
					}
					if (foundk2){ // a valid next kmer has been found /!\ sometimes the sequencing error re-creates a kmer present elsewhere and we are misled
						lastTrueKmer = kmer2;
						vector<uint> vpositionsR;
						if (posiNext < posi + k){
							vpositionsR = positionsRight(posi + 1 , k); // k2 not too far from k1 (still a junction kmer) // filled with k-1 positions
						} else {
							//~ cout << "R error " << posiNext << " " << vectSequences[r].substr(posiNext, k) << " " << vectSequences[r].substr(posiNext + k -1 - 1, k) << endl;
							vpositionsR = positionsRight(posiNext, k); // in case there are errors and k2 away from k1, not a junction kmer /!\ WARNING: not sure to detect the junction in this case, because we compare the signatures of 2 exons that might be used in the same quantity in the data set 
						}
						vector<uint> vpositionsL = positionsLeft(posi, k); // filled with k-1 positions
						unordered_set<int> zerosOfIntervalR,  zerosOfIntervalL;
						if (not vpositionsR.empty() and not vpositionsL.empty()){
							//~ cout << "oui" << endl;
							zerosOfIntervalR = zerosOfInterval(vectSequences, r, vpositionsR, k, kmerToReadPosi, kmerToSignature);
							zerosOfIntervalL =  zerosOfInterval(vectSequences, r, vpositionsL, k, kmerToReadPosi, kmerToSignature);
							if (zerosOfIntervalL != zerosOfIntervalR){
								cout << "junction !! at position "  << posi + k - 1 << " (after kmer "<< kmer1  <<")" <<endl;
								propagateSignature(signatures, lastSignature, posi + k);
								lastSignature = posi + k;
								confirmed = true;
							}
						}
						cout << kmer1 << " " << kmer2 << endl;
					} // else we can't do nothing, there are too many errors after kmer1
					else { cout << "no k2 opening" << endl;}
				} else {
					auto got1 = kmerLeft.find(kmer1);
					if (got1 != kmerLeft.end()){
						//// ** End of alternative exons (closing a bubble) ////
						// try and get an errorless kmer2 which is a previous kmer
						posiPrev = posi - 1;
						kmer2 = vectSequences[r].substr(posiPrev, k);
						auto gotNE2(kmerToReadPosi.find(kmer2));
						bool foundk2 = false;
						if (gotNE2 != kmerToReadPosi.end()){
							foundk2	= true;
							cout << "k2 closing ok" << endl;
						}
						if (not foundk2){
							while (posiPrev > posi - 1 - k - SIZE_INSERTION){
								--posiPrev;
								kmer2 = vectSequences[r].substr(posiPrev, k);
								auto gotNE2 = kmerToReadPosi.find(kmer2);
								if ( gotNE2 != kmerToReadPosi.end() ){
									foundk2 = true;
									cout << "new k2 closing ok" << endl;
									break;
								}
							}
						}
						if (foundk2){ // a valid previous kmer has been found /!\ sometimes the sequencing error re-creates a kmer present elsewhere and we are misled
							vector<uint> vpositionsL;
							if (posiNext < posi + k){
								vpositionsL = positionsLeft(posi - 1 , k); // k2 not too far from k1 (still a junction kmer) // filled with k-1 positions
							} else {
								vpositionsL = positionsLeft(posiNext, k); // in case there are errors and k2 away from k1, not a junction kmer /!\ WARNING: not sure to detect the junction in this case, because we compare the signatures of 2 exons that might be used in the same quantity in the data set 
							}
							vector<uint> vpositionsR = positionsRight(posi, k); // filled with k-1 positions
							unordered_set<int> zerosOfIntervalR,  zerosOfIntervalL;
							if (not vpositionsR.empty() and not vpositionsL.empty()){
								zerosOfIntervalR = zerosOfInterval(vectSequences, r, vpositionsR, k, kmerToReadPosi, kmerToSignature);
								zerosOfIntervalL =  zerosOfInterval(vectSequences, r, vpositionsL, k, kmerToReadPosi, kmerToSignature);
								if (zerosOfIntervalL != zerosOfIntervalR){
									cout << "junction !! at position "  << posi << " (before kmer "<< kmer1  <<")" <<endl;
									propagateSignature(signatures, lastSignature, posi);
									lastSignature = posi;
									confirmed = true;
								}
							}
							cout << kmer2 << " " << kmer1 << endl;
						} // else we can't do nothing, there are too many errors after kmer1
						else { cout << "no k2 closing" << endl;}
					}
				}
			}
			if (confirmed){
				posi += k;
				confirmed = false;
			} else {
				++posi;
			}
		}
		cout << "read " << r << endl;
		for (uint s : signatures){
			cout << s;
		}
		cout << endl;
	}
}


int main(int argc, char **argv){
	vector<string> ref({"TATTCTTTTCTAGATATATTTTTTTGTACACTGGTTTTTTTTCAACATTTGCGGACTCTGTTTTTATTCGATACCTTAGCTCTACGACACTCAATATCCAGCGGTAGTTTACCATTTGCCAGCAGGTTTTGGAATGGCTTGACAGGACCCAAAGTGTGGGATCTGTCTTGCCCAGTATCTGGATATTGTGGGTCATGTATCGAAATCATTTGCTGAGCAGAAGATGGATTTCTATGTTCCTCGCTGATCAAAACCTGTGGAAATGCCAGACTTGTCGCAAGAAATACAGTAGCACGAGATGTTTGCCAGAACTCCTGTGAGAGGACAGAATCTGGCGGAATAAATAAAACTGGAAACCCC",
	"TATTCTTTTCTAGATATATTTTTTTGTACACTGGTTTTTTTTCAACATTTGCGGACTCTGTTTTTATTCGATACAAGAAATACAGTAGCACGAGATGTTTGCCAGAACTCCTGTGAGAGGACAGAATCTGGCGGAATAAATAAAACTGGAAACCCC"});
	
	vector<string> vectSequences({			"TATTCTTTTCTAGATATATTTTTTTGTACACTGGTTTTTTTTCAACATTTGCGGACTCTGTTTTTATTCGATACTTTAGCTCTACGACACTCAATATCCAGCGGTAGTTTACCATTTGCCAGCAGGTTTTGGAATGGCTTGACAGGACCCAAAGTGTGGGATCTGTCTTGCCCAGTATCTGGATATTGTGGGTCATGTATCGAAATCATTTGCTGAGCAGAAGATGGATTTCTATGTTCCTCGCTGATCAAAACCTGTGGAAATGCCAGACTTGTCGCAAGAAATACAGTAGCACGAGATGTTTGCCAGAACTCCTGTGAGAGGACAGAATCTGGCGGAATAAATAAAACTGGAAACCCC","TATTCTTTTCTAGATATATTTTTTTGTACACTGGTTTTTTTTCAACATTTGCGGACTCTGTTTTTATTCGATACCTTAGCTCTACGACACTCAATATCCAGCGGTAGTTTACCATTTGCCAGCAGGTTTTGGAATGGCTTGACAGGACCCAAAGTGTGGGATCTGTCTTGCCCAGTATCTGGATATTGTGGGTCATGTATCGAAATCATTTGCTGAGCAGAAGATGGATTTCTATGTTCCTCGCTGATCAAAACCTGTGGAAATGCCAGACTTGTCGCAAGAAATACAGTAGCACGAGATGTTTGCCAGAACTCCTGTGAGAGGACAGAATCTGGCGGAATAAATAAAACTGGAAACCCC","TATTCTTTTCTAGATATATTTTTTTGTACACTGGTTTTTTTTCAACATTTGCGGACTCTGTTTTTATTCGATACCTTAGCTCTACGACACTCAATATCCAGCGGTAGTTTACCATTTGCCAGCAGGTTTTGGAATGGCTTGACAGGACCCAAAGTGTGGGATCTGTCTTGCCCAGTATCTGGATATTGTGGGTCATGTATCGAAATCATTTGCTGAGCAGAAGATGGATTTCTATGTTCCTCGCTGATCAAAACCTGTGGAAATGCCAGACTTGTCGCAAGAAATACAGTAGCACGAGATGTTTGCCAGAACTCCTGTGAGAGGACAGAATCTGGCGGAATAAATAAAACTGGAAACCCC","TATTCTTTTCTAGATATATTTTTTTGTACACTGGTTTTTTTTCAACATTTGCGGACTCTGTTTTTATTCGATACCTTAGCTCTACGACACTCAATATCCAGCGGTAGTTTACCATTTGCCAGCAGGTTTTGGAATGGCTTGACAGGACCCAAAGTGTGGGATCTGTCTTGCCCAGTATCTGGATATTGTGGGTCATGTATCGAAATCATTTGCTGAGCAGAAGATGGATTTCTATGTTCCTCGCTGATCAAAACCTGTGGAAATGCCAGACTTGTCGCAAGAAATACAGTAGCACGAGATGTTTGCCAGAACTCCTGTGAGAGGACAGAATCTGGCGGAATAAATAAAACTGGAAACCCC","TATTCTTTTCTAGATATATTTTTTTGTACACTGGTTTTTTTTCAACATTTGCGGACTCTGTTTTTATTCGATCCCTTAGCTCTACGACACTCAATATCCAGCGGTAGTTTACCATTTGCCAGCAGGTTTTGGAATGGCTTGACAGGACCCAAAGTGTGGGATCTGTCTTGCCCAGTATCTGGATATTGTGGGTCATGTATCGAAATCATTTGCTGAGCAGAAGATGGATTTCTATGTTCCTCGCTGATCAAAACCTGTGGAAATGCCAGACTTGTCGCAAGAAATACAGTAGCACGAGATGTTTGCCAGAACTCCTGTGAGAGGACAGAATCTGGCGGAATAAATAAAACTGGAAACCCC",
	

	"TATTCTTTTCTAGATATATTTTTTTGTACACTGGTTTTTTTTCAACATTTGCGGACTCTGTTTTTATTCGATACAAGAAATACAGTAGCACGAGATGTTTGCCAGAACTCCTGTGAGAGGACAGAATCTGGCGGAATAAATAAAACTGGAAACCCC",
	"TATTCTTTTCTAGATATATTTTTTTGTACACTGGTTTTTTTTCAACATTTGCGGACTCTGTTTTTATTCGATACAAGAAATACAGTAGCACGAGATGTTTGCCAGAACTCCTGTGAGAGGACAGAATCTGGCGGAATAAATAAAACTGGAAACCCC",
	"TATTCTTTTCTAGATATATTTTTTTGTACACTGGTTTTTTTTCAACATTTGCGGACTCTGTTTTTATTCGATACAAGAAATACAGTAGCACGAGATGTTTGCCAGAACTCCTGTGAGAGGACAGAATCTGGCGGAATAAATAAAACTGGAAACCCC",
	"TATTCTTTTCTAGATATATTTTTTTGTACACTGGTTTTTTTTCAACATTTGCGGACTCTGTTTTTATTCGATACTAGAAATACAGTAGCACGAGATGTTTGCCAGAACTCCTGTGAGAGGACAGAATCTGGCGGAATAAATAAAACTGGAAACCCC",
	

	});
	unordered_map<string, vector<int>> kmerToSignature;
	unordered_map<string, vector<pair<uint, uint>>> kmerToReadPosi;
	unordered_map<uint, vector<uint>> readToErrorPosition;
	uint k(9);
	for (uint i(2); i < 3; ++i){
		cout << i << endl;
		fillSignatures(vectSequences, k, kmerToSignature);
		dispatchErrors(vectSequences, k, kmerToSignature, kmerToReadPosi, readToErrorPosition);
	}

	
	unordered_set<string> branchingKmersLeft, branchingKmersRight;
	followingKmers(vectSequences, k, kmerToReadPosi, branchingKmersLeft, branchingKmersRight);
	//~ cout << "iiiiiiiiiiiiiii" << endl;
	//~ for (string i : branchingKmersLeft){
		//~ cout << i << endl;
	//~ }
	//~ cout << "iiiiiiiiiiiiiii" << endl;
	//~ for (string i : branchingKmersRight){
		//~ cout << i << endl;
	//~ }
	//~ cout << "iiiiiiiiiiiiiii" << endl;
	computeExons(vectSequences, k, kmerToReadPosi,  kmerToSignature,  branchingKmersLeft, branchingKmersRight);
	return 0;
}





