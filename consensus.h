#include <fstream>
#include <cstring>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>

using namespace std;

struct nucleotide{
	char n;
	double score;
	bool operator< (const nucleotide& nt2) const{
		return this->score < nt2.score;
	}
};


void setColumnsOfNt(const vector<string>& regionSet, string& targetRegion, vector <nucleotide>& nucl);
nucleotide sumForColumn(vector <char> columnNt);
nucleotide compareNtByScore(const nucleotide& nt1, const nucleotide& nt2);
string ntToString(const vector <nucleotide>& vec);
