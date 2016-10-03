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
	    group.push_back(element);
	    vector <uint> toMerge;
	    for (int i(0); i < secondSplit.size(); ++i){
		if (not secondSplit[i].empty()){
		    group.push_back(stoi(secondSplit[i])); // create a new group with elements from one line of the file
		}
	    }
	    bool inClust(false);
	    for (uint elts(0); elts < group.size(); ++elts){
		for (uint clust(0); clust < clusterVec.size(); ++clust){
		    //~ vect <uint> toMerge;
		    if (clusterVec[clust].count(group[elts])){
			inClust = true;
			toMerge.push_back(clust);
		    }
		}
	    }
	    if (not inClust){ // create a new cluster
		unordered_set<int> cluster;
		for (int j(0); j < group.size(); ++j){
		    cluster.insert(group[j]);
		    clusterVec.push_back(cluster);
		}
	    } else {
		if (toMerge.size() <= 1){ // merge group in cluster
		    
		} else { // merge several clusters + group
		    for (uint i(1); i < toMerge.size() ; ++i){
			clusterVec[toMerge[0]].insert(clusterVec[toMerge[i]].begin(), clusterVec[toMerge[i]].end());
			clusterVec[toMerge[i]].clear();
			for (int j(0); j < group.size(); ++j){
			    clusterVec[toMerge[0]].insert(group[j]);
			}
		    }
		}
	    }
	    //~ if (treatedElements.count(element)){
		//~ for (int i(0); i < pseudoClusterVec.size(); ++i){
		    //~ if (pseudoClusterVec[i].count(element)){
			//~ for (int j(0); j < group.size(); ++j){
			    //~ treatedElements.insert(group[j]);
			    //~ pseudoClusterVec[i].insert(group[j]);
			//~ }
		    //~ }
		//~ }
	    //~ } else {
		//~ unordered_set<int> pseudoCluster;
		//~ pseudoCluster.insert(element);
		//~ treatedElements.insert(element);
		//~ for (int j(0); j < group.size(); ++j){
		    //~ treatedElements.insert(group[j]);
		    //~ pseudoCluster.insert(group[j]);
		//~ }
		//~ pseudoClusterVec.push_back(pseudoCluster);
	    //~ }
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


