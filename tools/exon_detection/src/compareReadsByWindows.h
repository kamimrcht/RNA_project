#include <fstream>
#include <cstring>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include "detectJumps.h"

#ifndef COMPARE_READS
#define COMPARE_READS

using namespace std;

struct window{
	uint index;
	uint read;

    bool operator==(const window& a) const
	{
		return (index == a.index && read == a.read);
	}
};


struct compareWindow{
    bool operator()(const window& win1, const window& win2){
        return win1.index <win2.index;
    }
};

uint64_t transformWindowToHash(window win);

namespace std { template <> struct hash<window> {
	typedef window argument_type;
	typedef uint64_t result_type; uint64_t operator()(window key) const { return transformWindowToHash(key); } }; }


vector<window> removeDuplicates(const vector<window>& vect);
void setKmersToWindows(uint indexRead, uint indexWindow, string kmer, const unordered_map<string, uint>& solidKmers, unordered_map <string, vector <window>>& kmerToWindows);
void getStartKmerEnum(int& position, uint& startForKmer, uint sizeSeq, uint k, uint w);
void getKmersinFromReadsInMap(uint k, const vector <string>& readSet, unordered_map <string, uint>& kmersFromFile);
void getKmersinWindowsFromReads(uint k, uint w, const vector <string>& readSet, const unordered_map<string, uint>& solidKmers, unordered_map <string, vector<window>>& kmerToWindows);
void getKmersinFromReadsInMap(uint k, const vector <string>& readSet, unordered_map <string, uint>& kmersFromFile);
uint getSimilarity(uint posiForKmer, uint& nbKmers, uint indexReadTarget, const string& readSequence, const unordered_map<string, uint> solidKmers, const unordered_map<string, vector<window>> kmerToWindows, unordered_map<window, double>& similarity, uint k, uint w);
void getSolidKmers(const unordered_map<string, uint>& kmersFromFile, unordered_map<string, uint>& solidKmers);
void getSimilarWindowsPairs(uint readSetSize, uint indexWindowTarget, unordered_map<window, double>& similarity, uint nbKmers, unordered_map <uint, vector <pairOfIndexWindow>>& pairsOfaRead, uint w);
void compareAndCorrectReadWindows(uint k, uint w, vector<string>& readSet,const unordered_map<string, uint> solidKmers, const unordered_map<string, vector<window>> kmerToWindows);

#endif

