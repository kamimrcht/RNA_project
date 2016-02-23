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


void correctConsecutiveWindows(string & readSequence, const string& replacementSeq, uint w, uint k, uint firstIndexWin, uint lastIndexWin){
	string newSeq;
	if (replacementSeq.size() > w){
		if (replacementSeq.size() == readSequence.size()){ //  case 1: the whole sequence is replaced
			//~ newSeq = replacementSeq;
			readSequence = replacementSeq;
			//~ cout << "case 1" << endl;
		} else if (firstIndexWin == 0){ // case 2: the beginning is replaced
			newSeq = replacementSeq;
			newSeq += readSequence.substr((lastIndexWin+1)*w);
			readSequence = newSeq;
			//~ cout << "case 2 " << endl;
		}
		else {
			uint span(0);
			bool lastWindowIncluded(false);
			if (w * lastIndexWin + w == readSequence.size()){
				lastWindowIncluded = true;
				span = w * firstIndexWin;
			} else {
				if (not (w * (lastIndexWin + 1) + (w + k -1)> readSequence.size())){
					lastWindowIncluded = true;
					if (firstIndexWin == lastIndexWin) {
						span = readSequence.size() - w - k + 2;
					} else {
						span =  w * firstIndexWin;
					}
				} else if (w * lastIndexWin + w > readSequence.size()){
					lastWindowIncluded = true;
					span = w * firstIndexWin;
				}
			}
			if (not lastWindowIncluded){ // case 3: some region in the middle is replaced
				//~ cout << "case 3" << endl;
				newSeq = readSequence.substr(0, firstIndexWin*w);
				newSeq += replacementSeq;
				newSeq += readSequence.substr((lastIndexWin+1)*w);
				readSequence = newSeq;
			} else { // case 4 : the end of the sequence is replaced
				//~ cout <<  " case 4" << endl;
				newSeq = readSequence.substr(0, span);
				newSeq += replacementSeq;
				readSequence = newSeq;
			}
		}
	}
}
