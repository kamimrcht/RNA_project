#include <fstream>
#include <cstring>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>


using namespace std;


struct pairOfIndexWindow{
	uint target;
	uint query;
};


struct regionInRead{
	uint firstWindow;
	uint lastWindow;
	uint read;
	 bool operator==(const regionInRead& a) const
	{
		return (firstWindow == a.firstWindow and lastWindow == a.lastWindow and read==a.read);
	}
};


struct comparePairOfIndexWindow{
    bool operator()(const pairOfIndexWindow& pair1, const pairOfIndexWindow& pair2){
        return pair1.query < pair2.query;
    }
};


uint64_t transformRegionInReadToHash(regionInRead region);


namespace std { template <> struct hash<regionInRead> {
	typedef regionInRead argument_type;
	typedef uint64_t result_type; uint64_t operator()(regionInRead key) const { return transformRegionInReadToHash(key); } };
	}

void detectJumps(const vector<pairOfIndexWindow>& vec, uint indexReadT, uint indexReadQ, unordered_map <regionInRead, vector<regionInRead>>& correspondance);
void consensusBetweenRegions(const unordered_map <regionInRead, vector<regionInRead>>& correspondance, vector<string>& readSet, uint w, uint k);
