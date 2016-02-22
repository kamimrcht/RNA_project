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


void correctConsecutiveWindows(const string & readSequence, const string& replacementSeq, string& newSeq, uint w, uint k, uint firstIndexWin, uint lastIndexWin){
	string region;
	if (replacementSeq.size() == readSequence.size()){ //  case 1: the whole sequence is replaced
		newSeq = replacementSeq;
		cout << "case 1" << endl;
	} else if (firstIndexWin == 0){ // case 2: the beginning is replaced
		newSeq += replacementSeq;
		newSeq += readSequence.substr((lastIndexWin+1)*w);
		
		cout << "case 2 " << replacementSeq.size() << " " <<readSequence.size() << endl;
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
				cout << "span in case 4 "<<lastIndexWin *w<< endl;
				span = w * firstIndexWin;
			}
		}
		//~ cout << firstIndexWin << " " << lastIndexWin << endl;
		//~ cout << "w * lastIndexWin + w" << w * lastIndexWin + w << endl;
		//~ cout << "w * (lastIndexWin+1) + (w + k -1))" << w * (lastIndexWin +1) + (w + k -1) << endl;
		//~ cout << "readSequence.size()" << readSequence.size() << endl;
		//~ cout << "bool" << lastWindowIncluded << endl;
		if (not lastWindowIncluded){ // case 3: some region in the middle is replaced
			cout << "case 3" << endl;
			newSeq += readSequence.substr(0, firstIndexWin*w);
			newSeq += replacementSeq;
			newSeq += readSequence.substr((lastIndexWin+1)*w);
		} else { // case 4 : the end of the sequence is replaced
			cout <<  " case 4" << " " << span << endl;
			newSeq += readSequence.substr(0, span);
			newSeq += replacementSeq;
		}
	}
}
