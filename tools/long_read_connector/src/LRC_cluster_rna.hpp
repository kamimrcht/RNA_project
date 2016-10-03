#ifndef _LRC_cluster_rna_HPP_
#define _LRC_cluster_rna_HPP_


#include <fstream>
#include <unordered_set>
#include <string>
#include <vector>
#include <iostream>

using namespace std;


class LRC_cluster_rna
{
    string inFile;
    string outFile;
//~ private:

    //~ string inFile;

public:
	
    // Constructor
	LRC_cluster_rna();
	LRC_cluster_rna(const string& fileIn, const string& fileOut);

    // Actual job done by the tool is here
    void split_src_lines(const string& line, char delim, vector<string>& vec);
    void make_clusters(ifstream& srcFile, vector<unordered_set<int>>& clusterVec);
    void write_result(ofstream& out, const vector<unordered_set<int>>& clusterVec);
    void execute();
};


#endif /* _LRC_cluster_rna_HPP_ */
