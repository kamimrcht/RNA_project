#ifndef _LRC_pseudocluster_rna_HPP_
#define _LRC_pseudocluster_rna_HPP_


#include <fstream>
#include <unordered_set>
#include <string>
#include <vector>
#include <iostream>

using namespace std;
//~ #include <gatb/gatb_core.hpp>
//~ #include "../../../thirdparty/IteratorKmerH5/IteratorKmerH5.hpp"
//~ #include "../../../thirdparty/quasi_dictionary/src/quasidictionary.h"
//~ #include "common.hpp"


//~ class LRC_cluster_rna : public Tool

class LRC_pseudocluster_rna
{
    string inFile;
    string outFile;
//~ private:

    //~ string inFile;

public:
	
    // Constructor
	LRC_pseudocluster_rna();
	LRC_pseudocluster_rna(const string& fileIn, const string& fileOut);

    // Actual job done by the tool is here
    void split_src_lines(const string& line, char delim, vector<string>& vec);
    void make_pseudo_clusters(ifstream& srcFile, vector<unordered_set<int>>& pseudoClusterVec, unordered_set<int>& treatedElements);
    void write_result(ofstream& out, const vector<unordered_set<int>>& pseudoClusterVec);
    void execute();
};


#endif /* _LRC_cluster_rna_HPP_ */
