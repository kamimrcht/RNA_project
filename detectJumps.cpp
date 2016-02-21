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
		//~ cout << vec[i].target << " " << vec[i].query << endl;
		regionInRead regionTarget({indexT, vec[i].target, indexReadT});
		if (i != vec.size()-1 and absolute(int(vec[i].target)-int(vec[i+1].target))>2){
			cout << "jump" << endl;
			//~ if (correspondance.count(regionTarget)){
				//~ correspondance[regionTarget].push_back({indexQ, vec[i].query, indexReadQ});
			//~ } else {
				//~ regionInRead r({indexQ, vec[i].query, indexReadQ});
				//~ vector <regionInRead> v;
				//~ v.push_back(r);
				//~ correspondance[regionTarget] = v;
			//~ }
			//~ cout << "jump in target" << indexT << "to" << vec[i+1].target << endl;
			indexT = vec[i+1].target;
			indexQ = vec[i+1].query;
		} else if (i != vec.size()-1 and absolute(int(vec[i].query)-int(vec[i+1].query))>2){
			cout << "jump" << endl;
			//~ if (correspondance.count(regionTarget)){
				//~ correspondance[regionTarget].push_back({indexQ, vec[i].query, indexReadQ});
			//~ } else {
				//~ regionInRead r({indexQ, vec[i].query, indexReadQ});
				//~ vector <regionInRead> v;
				//~ v.push_back(r);
				//~ correspondance[regionTarget] = v;
			//~ }
			//~ cout << "jump in query" << indexQ <<  "to" << vec[i+1].query <<  endl;
			indexT = vec[i+1].target;
			indexQ = vec[i+1].query;
		}
	}
	if (indexT==0){
		cout << "no jump" << endl;
		//~ regionInRead regionTarget({vec[0].target, vec.back().target, indexReadT});
		//~ if (correspondance.count(regionTarget)){
				//~ correspondance[regionTarget].push_back({vec[0].query, vec.back().query, indexReadQ});
		//~ } else {
				//~ regionInRead r({vec[0].query, vec.back().query, indexReadQ});
				//~ vector <regionInRead> v;
				//~ v.push_back(r);
				//~ correspondance[regionTarget] = v;
		//~ }
	}
}
