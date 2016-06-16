#include <fstream>
#include <cstring>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include "compareReadsByWindows.h"
#include "correctionGraph.h"
#include "utils.h"

using namespace std;


int main(int argc, char ** argv){
	if (argc < 4){
		/* test graph */
		string fileName = argv[1];
		ifstream readFile(fileName);
		vector <string> readSet;
		string sequence;
		while (not readFile.eof()){
            getline(readFile, sequence);
			getline(readFile, sequence);
			if (not sequence.empty()){
				readSet.push_back(sequence);
			}
		}
		Graph graph(4);
		//~ string readRegion0("ACGTAGCATAGATTGA");
		//~ string readRegion1("ACGTAGCATAGATTGA");
		//~ string readRegion2("ACGTAGCATAGATTGA");
		string readRegion0("ACGTA");
		string readRegion1("ACGTA");
		string readRegion2("ACGTA");
		vector <string> vecReads({readRegion0, readRegion1, readRegion2});
		uint bestKmer(graph.createGraphFromSetofRegions(vecReads));
		//~ uint bestKmer(graph.createGraphFromSetofRegions(readSet));
		//~ cout << "bk "<< bestKmer << endl;
		vector <Node*> backbone;
		graph.getBackBone(bestKmer, backbone);
		cout << bestKmer << " size bb " << backbone.size() << endl;
		graph.getStartingNodes();
		graph.sequences2dot();
		system("dot -Tpng out.dot > output.png");
		//~ graph.duplicateNode(0);
		graph.graphClear();
		/* end test */
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
			readSet.push_back(getCanonical(sequence));
		}
		unordered_map <string, uint> kmersFromFile; //  TODO: destroy
		unordered_map <string, uint> solidKmers;
		unordered_map <string, vector<window>> kmerToWindows;
		getKmersinFromReadsInMap(k, readSet, kmersFromFile);
		getSolidKmers(kmersFromFile, solidKmers);
		getKmersinWindowsFromReads(k, w, readSet, solidKmers, kmerToWindows);
		compareAndCorrectReadWindows(k, w, readSet, solidKmers, kmerToWindows);
	}
	return 0;
}

