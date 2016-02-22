#include <fstream>
#include <cstring>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include "consensus.h"

using namespace std;


// todo before funct call : check if size of regionSet items are the same then the size of targetRegion
void setColumnsOfNt(const vector<string>& regionSet, string& targetRegion, vector <nucleotide>& nucl){
	for (uint pos(0); pos < targetRegion.size(); ++pos){
		vector <char> columnNt;
		columnNt.push_back(targetRegion[pos]);
		for (uint index(0); index < regionSet.size(); ++index){
			columnNt.push_back(regionSet[index][pos]);
		}
		nucleotide consensus(sumForColumn(columnNt));
		nucl.push_back(consensus);
	}
}


nucleotide compareNtByScore(const nucleotide& nt1, const nucleotide& nt2){
	return  nt1.score >= nt2.score ? nt1 : nt2;
}


nucleotide sumForColumn(vector <char> columnNt){
	nucleotide A({'A',0});
	nucleotide C({'C',0});
	nucleotide T({'T',0});
	nucleotide G({'G',0});
	for (uint i(0); i < columnNt.size(); ++i){
		switch (columnNt[i]) {
			case 'A':
				++A.score;
				break;
			case 'T':
				++T.score;
				break;
			case 'C':
				++C.score;
				break;
			case 'G':
				++G.score;
				break;
		}
	}
	nucleotide consensus(max(A,max(T,max(C,G))));
	consensus.score /= (columnNt.size());
	return consensus;
}

string ntToString(const vector <nucleotide>& vec){
	string s;
	for(uint i(0); i < vec.size(); ++i){
		s += vec[i].n;
	}
	return s;
}
