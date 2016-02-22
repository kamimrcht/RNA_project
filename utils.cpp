#include <fstream>
#include <cstring>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>

using namespace std;

string revComp(const string& seq){
	string revCompSeq = "";
	int pos = seq.size()-1;
	char nt;
	do{
		nt = seq[pos];
		switch (nt) {
			case 'A':
				revCompSeq += 'T';
				break;
			case 'T':
				revCompSeq += 'A';
				break;
			case 'C':
				revCompSeq += 'G';
				break;
			case 'G':
				revCompSeq += 'C';
				break;
		}
		--pos;
	} while (pos>=0);
	return revCompSeq;
}


string getCanonical(const string& seq){
	string revCompSeq = revComp(seq);
	return min(seq, revCompSeq); 
}


string getKmer(const string& sequence, int posi, int k){
	return getCanonical(sequence.substr(posi, k));
}


string getSequenceInWindow(const string & readSequence, uint w, uint k, uint indexWin){
	uint position(indexWin*w);
	string region;
	if (position + k - 1 < readSequence.size()){
		region = readSequence.substr(position, w);
	} else {
		uint posi= readSequence.size() - w - k + 1;
		region = readSequence.substr(posi);
	}
	return region;
}


string getSequenceInConsecutiveWindows(const string & readSequence, uint w, uint k, uint firstIndexWin, uint lastIndexWin){
	string region;
	for (uint index(firstIndexWin); index <= lastIndexWin; ++index){
		region += getSequenceInWindow(readSequence, w, k, index);
	}
	return region;
}

uint absolute(int a){
	return  a >= 0 ? a : -a;
}
