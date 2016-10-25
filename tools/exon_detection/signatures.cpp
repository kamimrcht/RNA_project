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

uint THRESHOLD1(2);
uint THRESHOLD2(3);
uint THRESHOLD3(1);
uint D(2);

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

vector<string> generateAltKmer(string& kmer){
	vector<string> ntInt({"A", "C", "G", "T"});
	vector<string> result;
	char ntc(kmer.back());
	string nt; nt.push_back(ntc);
	for (string i : ntInt){
		if (i != nt){
			result.push_back(kmer.substr(0, kmer.size()-1) + i);
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






void followingKmers(vector<string>& vectSequences, uint k, unordered_map<string, vector<pair<uint, uint>>>& kmerToReadPosi, unordered_map<string, vector<string>>& kmerToPrevious, unordered_map<string, vector<string>>& kmerToNext){
	string kmer, prev, next;
	for (uint r(0); r < vectSequences.size(); ++r){
		for (uint p(0); p < vectSequences[r].size() - k + 1; ++p){
			kmer = vectSequences[r].substr(p, k);
			if (kmerToReadPosi.count(kmer)){
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
				auto got(kmerToPrevious.find(kmer));
				auto got2(kmerToNext.find(kmer));
				if (gotn != kmerToReadPosi.end()){
					if (got2 != kmerToNext.end()){
						got2->second.push_back(next);
					} else {
						vector<string> v={next};
						kmerToNext[kmer] = v;
					}
				}
				if (gotp != kmerToReadPosi.end()){
					if (got != kmerToPrevious.end()){
						//~ cout << "here" << endl;
						got->second.push_back(prev);
					} else {
						vector<string> v={prev};
						kmerToPrevious[kmer] = v;
					}
				}
			}
		}
	}
}



uint nbReadsSupporting(string& kmer2, unordered_map<string, vector<string>>::const_iterator kmer){
	uint count(0), others(0);
	//~ cout << "func" << endl;
	for (string s : kmer->second){
		//~ cout << s << " " << kmer2 << endl;
		if (s == kmer2){
			//~ cout << "ok" << endl;
			++count;
		} else {
			++others;
		}
	}
	if (others > THRESHOLD3){
		//~ cout << count << endl;
		return count;
	} else { // means that there is only one kmer preceeding or following = no junction
		return 0;
	}
}



void computeExons(vector<string>& vectSequences, uint k,  unordered_map<string, vector<pair<uint, uint>>>& kmerToReadPosi,  unordered_map<string, vector<int>>& kmerToSignature, unordered_map<string, vector<string>>& kmerToPrevious, unordered_map<string, vector<string>>& kmerToNext){
	string kmer, kmer2, lastKmer("");
	
	for (uint r(0); r < vectSequences.size(); ++r){
		uint posi(0), posi2(1), previousDetectedPosi = 0;
		vector <uint> signatures;
		signatures.push_back(0); // signature of the first kmer
		while (posi < vectSequences[r].size() - k and posi2 < vectSequences[r].size() - k + 1){
			kmer = vectSequences[r].substr(posi, k);
			lastKmer = kmer;
			auto got(kmerToReadPosi.find(kmer));
			if (got != kmerToReadPosi.end()){
				kmer2 = vectSequences[r].substr(posi + 1, k);
				auto got2(kmerToReadPosi.find(kmer2));
				if (got2 != kmerToReadPosi.end()){ // both kmers errorless - compare signatures
					lastKmer = kmer2;
					vector<int> sign1 = kmerToSignature[kmer];
					vector<int> sign2 = kmerToSignature[kmer2];
					vector<int> subs;
					transform(sign1.begin(), sign1.end(), sign2.begin(), std::back_inserter(subs), [](int a, int b) { return abs(a-b); });
					uint sum(0);
					for (uint n : subs){
						sum+=n;
					}
					if (sum < THRESHOLD2){
						signatures.push_back(signatures.back());
					} else {
						signatures.push_back(signatures.back() + 1); // change the signature
						if (posi != vectSequences.size() - k and posi != previousDetectedPosi){
							unordered_map<string, vector<string>>::const_iterator gotn (kmerToNext.find(kmer));
							bool lookForPrev(false);
							if (gotn != kmerToNext.end()){
								uint nbReadsS(nbReadsSupporting(kmer2, gotn));
								if (nbReadsS > 0){ // there are several types of kmers following : junction
									cout << "junction at position " << posi + k -1<< endl;
									previousDetectedPosi = posi + k - 1;
								} else if (not posi == 0){
									lookForPrev = true;
								}
							} else  if (not posi < k){
								lookForPrev = true;
							}
							if (lookForPrev){
								unordered_map<string, vector<string>>::const_iterator gotp (kmerToPrevious.find(kmer2));
								if (gotp != kmerToPrevious.end()){
									string prevk(kmer);
									//~ string prevk(vectSequences[r].substr(posi2, k));
									uint nbReadsS(nbReadsSupporting(prevk, gotp));
									if (nbReadsS > 0){ // there are several types of kmers following : junction
										cout << "junction at position *" << posi << endl;
									}
								}
								lookForPrev = false;
							}
						}
					}
				} else { // second kmer with error : we propagate the signature
					signatures.push_back(signatures.back());
				}
			} else { // first kmer has an error
				auto got2(kmerToReadPosi.find(kmer));
				if (got2 != kmerToReadPosi.end()){
				// compare to the last errorless kmer's sign if it exists.
					if (not lastKmer.empty()){
						vector<int> sign1 = kmerToSignature[lastKmer];
						vector<int> sign2 = kmerToSignature[kmer2];
						vector<int> subs;
						transform(sign1.begin(), sign1.end(), sign2.begin(), std::back_inserter(subs), [](int a, int b) { return abs(a-b); });
						uint sum(0);
						for (uint n : subs){
							sum+=n;
						}
						if (sum < THRESHOLD2){
							signatures.push_back(signatures.back());
						} else {
							signatures.push_back(signatures.back() + 1); // change the signature
						}
					} else {
						signatures.push_back(signatures.back());
					}
				} else {
					signatures.push_back(signatures.back());
				}
			}
			++posi; ++posi2;
		}
		cout << "read " << r << endl;
		for (uint s : signatures){
			cout << s ;
		}
		cout << endl;
	}
	
}


void correctErrors(vector<string>& vectSequences, unordered_map<uint, vector<uint>>& readToErrorPosition, uint k, unordered_map<string, vector<pair<uint, uint>>>& kmerToReadPosi,  uint threshold=D){
	string kmer, kmer2, prevKmer;
	vector<string> altKmers(3);
	uint sizeSubSeq;
	uint nextTruePosi;
	bool okCorrect;
	string nextTrueKmer;
	string consensus;
	//~ bool first(false);
	uint pos(0);
	for (auto r(readToErrorPosition.begin()); r != readToErrorPosition.end(); ++r){
		for (uint posi1 : r->second){
			kmer = vectSequences[r->first].substr(posi1, k);
			auto gotPrev(kmerToReadPosi.end());
			if (posi1 != 0){
				prevKmer = vectSequences[r->first].substr(posi1 - 1, k);
				gotPrev = kmerToReadPosi.find(prevKmer);
			}
			if (gotPrev != kmerToReadPosi.end()){ // correct only if the precedent kmer is ok
				okCorrect = false;
				for (uint posi2(posi1 + 1); posi2 < posi1 + (k -1) + threshold + 1; ++posi2){
					kmer2 = vectSequences[r->first].substr(posi2, k);
					auto got = kmerToReadPosi.find(kmer2);
					if (got != kmerToReadPosi.end()){ // this kmer is not erroneous
						sizeSubSeq = posi2 - posi1 + k - 1;
						nextTrueKmer = kmer2;
						nextTruePosi = posi2 + k - 1;
						okCorrect = true;
						break;
					} // vérifier ici que ça marche avec une erreur trop grande
				}
				if (okCorrect){
					bool corrected(false);
					string consensus("");
					altKmers = generateAltKmer(kmer);
					for (string alt : altKmers){
						vector<string> recruited;
						vector<uint> vecSizes;
						auto got = kmerToReadPosi.find(alt);
						auto got2 = kmerToReadPosi.find(nextTrueKmer);
						if(got != kmerToReadPosi.end() and got2 != kmerToReadPosi.end() ){
							for(auto it = got->second.begin(); it != got->second.end(); ++it ){
								uint read = it->first;
								uint position1 = it->second;
								//todo : instead of taking the next true kmer, mutate the last false kmer 
								string recruitedSeq(returnSeqBetween2Kmers(k, read, position1, nextTrueKmer, vectSequences)); // get the sequence between the two flanking right kmers in reads that contain them
								if (not recruitedSeq.empty()){
									recruited.push_back(recruitedSeq);
									vecSizes.push_back(recruitedSeq.size());
								}
							}
							if (recruited.size() > 1){
								uint lenToKeep(sizeToKeep(vecSizes));
								if (lenToKeep != 0){
									vector<string> recruited2;
									for (string s : recruited){
										if (s.size() == lenToKeep){
											recruited2.push_back(s);
										}
									}
									consensus = makeConsensus(recruited2);
									if (consensus.size() > 0){
										corrected = true;
									}
								}
							} else if (not recruited.empty()) { // discutable
								consensus = recruited.back();
								corrected = true;
							}
							if (corrected){ // greedy
								vectSequences[r->first] = vectSequences[r->first].substr(0, posi1) + consensus + vectSequences[r->first].substr(nextTruePosi);
								corrected = false;
								okCorrect = false;
								break; //TODO: break if corrected = T
							}
						}
					}
				}
				pos = posi1;
			}
		}
	}
}


//~ void addNewReads(vector<string>& vectSequences, vector<uint>& correctedReads, vector<string>& newReads){
	//~ for (uint r(0); r < correctedReads.size(); ++r){
		//~ vectSequences[correctedReads[r]] = newReads[r]; 
	//~ }
//~ }

int main(int argc, char **argv){
	//~ vector<string> vectSequences({"AGGTTTAATCTG","AGGTTGGGGTAACTG","AGGTTTAACTG","AGGTTTACT","AGGTTTAACTG"}); //todo: kmers du debut et de la fin...
	string ref("TATTCTTTTCTAGATATATTTTTTTGTACACTGGTTTTTTTTCAACATTTGCGGACTCTGTTTTTATTCGATACCTTAGCTCTACGACACTCAATATCCAGCGGTAGTTTACCATTTGCCAGCAGGTTTTGGAATGGCTTGACAGGACCCAAAGTGTGGGATCTGTCTTGCCCAGTATCTGGATATTGTGGGTCATGTATCGAAATCATTTGCTGAGCAGAAGATGGATTTCTATGTTCCTCGCTGATCAAAACCTGTGGAAATGCCAGACTTGTCGCAAGAAATACAGTAGCACGAGATGTTTGCCAGAACTCCTGTGAGAGGACAGAATCTGGCGGAATAAATAAAACTGGAAACCCC");
	
	//~ vector<string> vectSequences({			"TATTCTTTTCTAGATATATTTTTTTGTACACTGGTTTTTTTTCAACATTTGCGGACTCTGTTTTTTTATTCGATACCTTAGCTCTACGACACTCAATATCCAGCGGTAGTTTACCATTTGCCAGCAGGTTTTGGAATGGCTTGACAGGACCCAAAGTGTGGGATCTGTCTTGCCCAGTATCTGGATATTGTGGGTCATGTATCGAAATCATTTGCTGAGCAGAAGATGGATTTCCCTATGTTCCTCGCCTGATCAAAACCTGTGGAAATGCCAGACTTGTCGCAAGAAAGTACAGTAGCACGAGATGTTTGCCAGAACTCCTGTGAGAGGACAGAATCTGGCGGAATAAATAAAACTGGAAACCCC",

	//~ "TATTCTTTTCTAGATATATTTTTTTGTACACTGGTTTTTTTTCAACATTTGCGGACTCTGTTTTTATTCGATAACCTTAGCTCTACGACACTCAATATCCAGCGGTAGTTTACCATTTGCCAGCAGGTTTTGGAATGGCTTGACAGGACCCAAAGTGTGGGATCTGTCTTGCCCAGTATCTGGATATTGTGGGTCATGTATCGAAATCATTTGCTGAGCAGAAGATGGATTTCTATGTTCCTCGCTGATCAAAACCTGTGGAAATGCCAGACTTGTCGCAAGAAATACAGTAGCACGAGATGTTTGCCAGAACTCCTGTGATGAGGACAGAATCTGGCGGAATAAATAAAACTGGAAACCCC",

	//~ "TATTCTTTTCTAGATATATTTTTTTGTACACTGGTTTTTTTGTCAACATTTGCGCACTCTGTTTCCTATTCGATACCTTAGCTCTACGACACTCAATATCCAGCGGTAGTTTACCATTTGCCAGCAGGTTTTGGGAATGGCTTGACAGGACCCAAAGTGTGGGATCTGTCTTGCCCAGTATCTGGATATTGTGGGTCATGTATCGAAATCATTTGCTGTGCAGAAGATGGATTTCTATGTTCCTCGCTGATCAAAACCTGTGGAAATGCCAGACTTGTCGCAAGAAATACAGTAGCACGAGATGTTTGCCAGAACTCCTGTGAGAGGACAGAATCTGGCGGAATAAATAAAACTGGAAACCCC",

	//~ "TATTCTTTTCTAGATATATTTTTTTGTACACTGGTTTTTTTTCAACATTTGCGGACTCTGTTTTTATTCGATACCTTAGCTCTACGACACTCAATATCCAGCGGTAGTTTACCATTTGCCAAGCAGGTTTTGGAATGGCTTGACAGGACCCAAAGTGTGGGATCTGTCTTGCCCAGTATCTGGAGCCTATTGTGGGTCATGTATCGAAATCATTTGCTGAGCAGAAGATGGATTTCTATGTTCCTCGCTGATCAAAACCTGTGGAAATGCCAGACTTGTCGCAAGAAATACAGTAGCACGAGATGTTTGCCAGAACTCCTGTGAGAGGACAGAATCTGGCGGAATAAATAAAACTGGAAACCCC",

	//~ "TATTCTTTTCTAGATATATTTTTTGTACACTGGTTTTTTTTCAACATTTGCGGACTCTGTTTTTATTCGATACCTTAGCTCTACGACACTCAATATCCAGCGGTAGTTTACCATTTGCCAGCAGGTTTTGGAATGGCTTGACAGGACCCAAAGTGTGGGATCTGTCTTGCCCAGTATCTGGATATTGTGGGTCATGTATCGAAATCATTTGCTGAGCAGAAGATGGATTTCTATGTTCCTCGCTGATCAAAACCTGTGGAAATGCCAGACTTGTCGCAAGAAATACAGTAGCACGAGATGTTTGCCAGAACTCCTGTGAGAGGACAGAATCTGGCGGAATAAATAAAACTGGAAACCCC"});
	

	vector<string> vectSequences({"TATTCTTTTCTAGAATATATTTTTTTGTACACTGGTTTTTTTTCAACATTTGCGGACTCTGTTTTTATTCGATGCAAGAAATACAGTAGCACGAGATGTTTGCCAGAACTCCTGTGAGAGGACAGAATCTGGCGGAATAAATAAAACTGGAAACCCC",


	"TATTCTTTTCTAGATATATTTTTTTGTACACTGGTTTTTTTTCAACATTTGCGGAACTCTGTTTTTATTCGATGCAAGAAATACAGTAGCACGAGATGTTTGCCAGAACTCCTGTGAGAGGACAGAATCTGGCGGAATAAATAAAACTGGAAACCCC",
	"TATTCTTTTCTAGATATATTTTTTTGTACACTGGTTTTTTTTCAACATTTGCGGACTCTGTTTTTATTCGATGCAAGAAATACAGTAGCACGAGATGTTTGCCAGAACTCCTGTGAGAGGACAGAATCTGGCGGAATAAATAAAACTGGAAACCCC",
	"TATTCTTTTCTAGATATATTTTTTTGTACACTGGTTTTTTTTCAACATTTGCGGACTCTGTTTTTATTCGATGCAAGAAATACAGTAGCACGAGATGTTTGCCAGAACTCCTGTGAGAGGACAGAATCTGGCGGAATAAATAAAACTGGAAACCCC",

	"TATTCTTTTCTAGATATATTTTTTTGTACACTGGTTTTTTTTCAACATTTGCGGACTCTGTTTTTATTCGATACCTTAGCTCTACGACACTCAATATCCAGCGGTAGTTTACCATTTGCCAGCAGGTTTTGGAATGGCTTGACAGGACCCAAAGTGTGGGATCTGTCTTGCCCAGTATCTGGATATTGTGGGTCATGTATCGAAATCATTTGCTGAGCAGAAGATGGATTTCTATGTTCCTCGCTGATCAAAACCTGTGGAAATGCCAGACTTGTCGCAAGAAATACAGTAGCACGAGATGTTTGCCAGAACTCCTGTGAGAGGACAGAATCTGGCGGAATAAATAAAACTGGAAACCCC",
	"TATTCTTTTCTAGATATATTTTTTTGTACACTGGTTTTTTTTCAACATTTGCGGACTCTGTTTTTATTCGATACCTTAGCTCTACGACACTCAATATCCAGCGGTAGTTTACCATTTGCCAGCAGGTTTTGGAATGGCTTGACAGGACCCAAAGTGTGGGATCTGTCTTGCCCAGTATCTGGATATTGTGGGTCATGTATCGAAATCATTTGCTGAGCAGAAGATGGATTTCTATGTTCCTCGCTGATCAAAACCTGTGGAAATGCCAGACTTGTCCAAGAAATACAGTAGCACGAGATGTTTGCCAGAACTCCTGTGAGAGGACAGAATCTGGCGGAATAAATAAAACTGGAAACCCC",

	"TATTCTTTTCTAGATATATTTTTTTGTACACTGGTTTTTTTTCAACATTTGCGGACTCTGTTTTTATTCGATACCTTAGCTCTACGACACTCAATATCCAGCGGTAGTTTACCATTTGCCAGCAGGTTTTGGAATGGCTTGACAGGACCCAAAGTGTGGGATCTGTCTTGCCCAGTATCTGGATATTGTGGGTCATGTATCGAAATCATTTGCTGAGCAGAAGATGGATTTCTATGTTCCTCGCTGATCAAAACCTGTGGAAATGCCAGACTTGTCGCAAGAAATACAGTAGCACGAGATGTTTGCCAGAACTCCTGTGAGAGGACAGAATCTGGCGGAATAAATAAAACTGGAAACCCC",

	"TATTCTTTTCTAGATATATTTTTTTGTACACTGGTTTTTTTTCAACATTTGCGGACTCTGTTTTTATTCGATACCTTAGCTCTACGACACTCAATATCCAGCGGTAGTTTACCATTTGCCAGCAGGTTTTGGAATGGCTTGACAGGACCCAAAGTGTGGGATCTGTCTTGCCCAGTATCTGGATATTGTGGGTCATGTATCGAAATCATTTGCTGAGCAGAAGATGGATTTCTATGTTCCTCGCTGATCAAAACCTGTGGAAATGCCAGACTTGTCGCAAGAAATACAGTAGCACGAGATGTTTGCCAGAACTCCTGTGAGAGGACAGAATCTGGCGGAATAAATAAAACTGGAAACCCC"});
	//pas de correction dans le stretch de T
	
	//todo: kmers du debut et de la fin...
	uint k(9);
	//~ for (uint i(2); i < 6; ++i){
		unordered_map<string, vector<int>> kmerToSignature;
		fillSignatures(vectSequences, k, kmerToSignature);
		unordered_map<string, vector<pair<uint, uint>>> kmerToReadPosi;
		unordered_map<uint, vector<uint>> readToErrorPosition;
		dispatchErrors(vectSequences, k, kmerToSignature, kmerToReadPosi, readToErrorPosition);
		unordered_map<string, vector<string>> kmerToPrevious, kmerToNext;
		followingKmers(vectSequences, k, kmerToReadPosi, kmerToPrevious, kmerToNext);
		computeExons(vectSequences, k, kmerToReadPosi,  kmerToSignature,  kmerToPrevious, kmerToNext);
		// todo : garder les kmers des junctions et les positions, faire un consensus de la position de jonction

		//~ correctErrors(vectSequences, readToErrorPosition, k, kmerToReadPosi,  i);
	//~ }
	// une seconde passe après la correction pour des signatures plus propres => description des exons
	for (string s : vectSequences){
		if (ref == s){
			cout << "ok" << endl;
		} else {
			cout << "pas ok" << endl;
		}
	}
	return 0;
}
