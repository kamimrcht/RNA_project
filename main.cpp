#include <fstream>
#include <cstring>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include "compareReadsByWindows.h"
#include "utils.h"

using namespace std;


int main(int argc, char ** argv){
	if (argc < 4){
		cout << "command line: ./rnaLR reads.fasta k w" << endl;
	} else {
		string fileName = argv[1];
		uint k = stoi(argv[2]);
		uint w = stoi(argv[3]);
		ifstream readFile(fileName);
		vector <string> readSet;
		string sequence;
		while (not readFile.eof()){
            getline(readFile, sequence);
			getline(readFile, sequence);
			readSet.push_back(sequence);
		}
		unordered_map <string, uint> kmersFromFile; //  TODO: destroy
		unordered_map <string, uint> solidKmers;
		unordered_map <string, vector<window>> kmerToWindows;
		getKmersinFromReadsInMap(k, readSet, kmersFromFile);
		getSolidKmers(kmersFromFile, solidKmers);
		getKmersinWindowsFromReads(k, w, readSet, solidKmers, kmerToWindows);
		compareReadWindows(k, w, readSet, solidKmers, kmerToWindows);
	}
	return 0;
}

