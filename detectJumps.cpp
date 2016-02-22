#include <fstream>
#include <cstring>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include "detectJumps.h"
#include "utils.h"
#include "consensus.h"


using namespace std;


uint64_t transformRegionInReadToHash(regionInRead region){
	hash<int> regionHash;
	return regionHash(region.read + region.firstWindow + region.lastWindow);
}


void detectJumps(const vector<pairOfIndexWindow>& vec, uint indexReadT, uint indexReadQ, unordered_map <regionInRead, vector<regionInRead>>& correspondance){
	uint indexT(0);
	uint indexQ(0);
	for (uint i(0); i<vec.size(); ++i){
		regionInRead regionTarget({indexT, vec[i].target, indexReadT});
		if (i != vec.size()-1 and absolute(int(vec[i].target)-int(vec[i+1].target))>2){
			if (correspondance.count(regionTarget)){
				correspondance[regionTarget].push_back({indexQ, vec[i].query, indexReadQ});
			} else {
				regionInRead r({indexQ, vec[i].query, indexReadQ});
				vector <regionInRead> v;
				v.push_back(r);
				correspondance[regionTarget] = v;
			}
			indexT = vec[i+1].target;
			indexQ = vec[i+1].query;
		}
		if (i != vec.size()-1 and absolute(int(vec[i].query)-int(vec[i+1].query))>2){
			if (correspondance.count(regionTarget)){
				correspondance[regionTarget].push_back({indexQ, vec[i].query, indexReadQ});
			} else {
				regionInRead r({indexQ, vec[i].query, indexReadQ});
				vector <regionInRead> v;
				v.push_back(r);
				correspondance[regionTarget] = v;
			}
			indexT = vec[i+1].target;
			indexQ = vec[i+1].query;
		}
	}
	regionInRead regionTarget({indexT, vec.back().target, indexReadT});
	if (correspondance.count(regionTarget)){
		correspondance[regionTarget].push_back({indexQ, vec.back().query, indexReadQ});
	} else {
		regionInRead r({indexQ, vec.back().query, indexReadQ});
		vector <regionInRead> v;
		v.push_back(r);
		correspondance[regionTarget] = v; 
	}
}


void consensusBetweenRegions(const unordered_map <regionInRead, vector<regionInRead>>& correspondance, vector<string>& readSet, uint w, uint k){
	for (auto iter = correspondance.begin(); iter != correspondance.end(); ++iter){
		if (iter->second.size()>1){
			vector <string> readRegionSeqs;
			//~ string targetSequence(getCanonical(readSet[iter->first.read]));
			string targetSequence(readSet[iter->first.read]);
			string targetRegion(getSequenceInConsecutiveWindows(targetSequence, w, k, iter->first.firstWindow, iter->first.lastWindow));
			for (uint i(0); i < iter->second.size(); ++i){
				string seq(getSequenceInConsecutiveWindows(readSet[iter->second[i].read], w, k, iter->second[i].firstWindow, iter->second[i].lastWindow));
				if (seq.size() == targetRegion.size()){
					readRegionSeqs.push_back(seq);
				}
			}
			vector <nucleotide> nucl;
			setColumnsOfNt(readRegionSeqs, targetRegion, nucl);
			string consensus(ntToString(nucl));
			//~ string newSeq;
			cout << "former " << targetSequence << endl;
			correctConsecutiveWindows(readSet[iter->first.read], consensus, w, k, iter->first.firstWindow, iter->first.lastWindow);
			cout << "new " << readSet[iter->first.read] << endl;
			cout << "SIZES" << targetSequence.size() << " " << readSet[iter->first.read].size() << endl;
		}
	}
}


