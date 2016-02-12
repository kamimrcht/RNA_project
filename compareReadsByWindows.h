#include <fstream>
#include <cstring>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>

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
//~ string revComp(const string& seq);
//~ string getCanonical(const string& seq);
//~ string getKmer(const string& sequence, int posi, int k);
void setKmersToWindows(uint indexRead, uint indexWindow, string kmer, const unordered_map<string, uint>& solidKmers, unordered_map <string, vector <window>>& kmerToWindows);

void getKmersinFromReadsInMap(uint k, const vector <string>& readSet, unordered_map <string, uint>& kmersFromFile);
void getKmersinWindowsFromReads(uint k, uint w, const vector <string>& readSet, const unordered_map<string, uint>& solidKmers, unordered_map <string, vector<window>>& kmerToWindows);
void getKmersinFromReadsInMap(uint k, const vector <string>& readSet, unordered_map <string, uint>& kmersFromFile);

void getSolidKmers(const unordered_map<string, uint>& kmersFromFile, unordered_map<string, uint>& solidKmers);
void compareReadWindows(uint k, uint w, const vector<string>& readSet,const unordered_map<string, uint> solidKmers, const unordered_map<string, vector<window>> kmerToWindows);

