#include <fstream>
#include <cstring>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include "detectJumps.h"
#include "utils.h"



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
				//~ cout << "correspondance " << indexT << " " << vec[i].target << " : "<< indexQ << vec[i].query << endl; 
			} else {
				regionInRead r({indexQ, vec[i].query, indexReadQ});
				vector <regionInRead> v;
				v.push_back(r);
				correspondance[regionTarget] = v;
				//~ cout << "correspondance " << indexT << " " << vec[i].target << " : "<< indexQ << " " << vec[i].query << endl; 
			}
			//~ cout << "jump in target: " << vec[i].target << " to " << vec[i+1].target << endl;
			indexT = vec[i+1].target;
			indexQ = vec[i+1].query;
		} else if (i != vec.size()-1 and absolute(int(vec[i].query)-int(vec[i+1].query))>2){
			if (correspondance.count(regionTarget)){
				correspondance[regionTarget].push_back({indexQ, vec[i].query, indexReadQ});
			} else {
				regionInRead r({indexQ, vec[i].query, indexReadQ});
				vector <regionInRead> v;
				v.push_back(r);
				correspondance[regionTarget] = v;
			}
			//~ cout << "jump in query: " << vec[i].query <<  " to " << vec[i+1].query <<  endl;
			indexT = vec[i+1].target;
			indexQ = vec[i+1].query;
		}
	}
	regionInRead regionTarget({indexT, vec.back().target, indexReadT});
	//~ if (indexT == 0 and indexQ == 0){
		//~ cout << "no jump" << endl;
	//~ }
	if (correspondance.count(regionTarget)){
		correspondance[regionTarget].push_back({indexQ, vec.back().query, indexReadQ});
		//~ cout << "correspondance " << indexT << " " << vec.back().target << " : "<< indexQ << vec.back().query << endl; 
	} else {
		regionInRead r({indexQ, vec.back().query, indexReadQ});
		vector <regionInRead> v;
		v.push_back(r);
		correspondance[regionTarget] = v;
		//~ cout << "correspondance " << indexT << " " << vec.back().target << " : "<< indexQ << " " << vec.back().query << endl; 
	}
}
