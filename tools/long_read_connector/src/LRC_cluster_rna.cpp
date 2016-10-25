#include <fstream>
#include <unordered_set>
#include <string>
#include <vector>
#include <iostream>
#include "LRC_cluster_rna.hpp"
using namespace std;


//~ static const char* STR_IN_FILE = "-f";


LRC_cluster_rna::LRC_cluster_rna(){
	inFile = "long_read_connector_res.tmp";
	outFile = "long_read_connector_res.txt";
}


LRC_cluster_rna::LRC_cluster_rna(const string& in, const string& out){
	inFile = in;
	outFile = out;
}


//~ LRC_cluster_rna::LRC_cluster_rna()  : Tool ("SRC_cluster_rna"){
	//~ getParser()->push_back (new OptionOneParam (STR_IN_FILE, "input_file",    true));
//~ }



void LRC_cluster_rna::split_src_lines(const string& line, char delim, vector<string>& vec) {
    auto i = 0;
    auto pos = line.find(delim);
    while (pos != string::npos) {
	vec.push_back(line.substr(i, pos-i));
	i = ++pos;
	pos = line.find(delim, pos);
	if (pos == string::npos)
		vec.push_back(line.substr(i, line.length()));
    }
}


void LRC_cluster_rna::make_clusters(ifstream& srcFile, vector<unordered_set<int>>& clusterVec){
//~ void LRC_cluster_rna::make_pseudo_clusters(ifstream& srcFile, vector<unordered_set<int>>& pseudoClusterVec, unordered_set<int>& treatedElements){
    string line;
    getline(srcFile, line); // header
    while (not srcFile.eof() ){
	getline(srcFile, line);
	vector<string> firstSplit;
	vector<string> secondSplit;
	if (not line.empty()){
	    split_src_lines(line, ':', firstSplit);
	    int element(stoi(firstSplit[0]));
	    split_src_lines(firstSplit[1], ' ', secondSplit);
	    vector<int> group;
	    unordered_set<int> test;
	    test.insert(element);
	    for (int i(0); i < secondSplit.size(); ++i){
		if (not secondSplit[i].empty()){
		    test.insert(stoi(secondSplit[i]));
		}
	    }
	    clusterVec.push_back(test);
	}
    }
    bool inClust(true);
    uint c(0), cc(0);
    while (inClust){
	
	cc = 0;
	for (uint i(0); i < clusterVec.size() ; ++i){
	    
	    if (not clusterVec[i].empty()){
		for (uint j(i+1); j < clusterVec.size(); ++j){
		    c = 0;
		    if (not clusterVec[j].empty()){
			    for (auto elt(clusterVec[i].begin()); elt != clusterVec[i].end(); ++elt){
				if (clusterVec[j].count(*elt)){
				    clusterVec[i].insert(clusterVec[j].begin(), clusterVec[j].end());
				    clusterVec[j].clear();
				    inClust = true;
				    goto test;
				} else {
				    ++c;
				}
			    }
			    if (c == clusterVec[i].size()){
				++ cc;
			    }
		    } else {
			++cc;
		    }
		}
		if (cc == clusterVec.size() - 1){
		    inClust = false;
		}
		test:
		{
		    if (inClust){
			break;
		    }	
		}
	    }
	}
    }
}


void LRC_cluster_rna::write_result(ofstream& out, const vector<unordered_set<int>>& clusterVec){
	cout << clusterVec.size() << endl;
	for (uint i(0); i < clusterVec.size(); ++i){
	    if (not clusterVec[i].empty()){
		string toPrint;
		bool first(true);
		for (auto j(clusterVec[i].begin()); j != clusterVec[i].end(); ++j){
			if (first){
				toPrint += to_string(*j);
				first = false;
			} else {
				toPrint += " " + to_string(*j);
			}
		}
		out << toPrint << endl;
	    }
	}
}

void LRC_cluster_rna::execute(){
	ofstream outputFile(outFile);
	//~ string srcFileName = getInput()->getStr(STR_IN_FILE).c_str();
	ifstream srcFile(inFile);
	//~ ifstream srcFile(srcFileName);
	vector<unordered_set<int>> clusterVec;
	make_clusters(srcFile, clusterVec);
	write_result(outputFile, clusterVec);
}


