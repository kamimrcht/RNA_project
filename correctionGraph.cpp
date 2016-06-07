#include <vector>
#include <string>
#include "correctionGraph.h"
#include <iostream>
#include <fstream>



using namespace std;
string nPrefix(uint n, const string& sequence){
	return sequence.substr(0, n);
}


string nSuffix(uint n,const string& sequence){
	return sequence.substr(sequence.size()-n, n);
}

// ----- Class Node ----- //

// Node constructor
Node::Node(uint indexI, const string& kmerI, uint position, uint readIndex){
	index = indexI;
	kmer = kmerI;
	readsAndPositions[readIndex].push_back(position);
}

// add a new position (and read) of a k-mer in the node
void Node::addPosition(uint position, uint readIndex){
	readsAndPositions[readIndex].push_back(position);
}

// add a new in-node to the node's neighbours (node )
void Node::addInNode(Node* nodeIn){
	inNodes.push_back(nodeIn);
}

// add a new out-node to the node's neighbours
void Node::addOutNode(Node* nodeOut){
	outNodes.push_back(nodeOut);
}
// ----- end Class Node ----- //




// ----- Class Graph ----- //


// graph constructor
Graph::Graph(uint kmerSize){
	k = kmerSize;
}

// store graph's starting node
void Graph::addFirstNode(Node* node){
	firstPositions.push_back(node);
}

//construct a graph with all fragments of reads corresponding to a region
void Graph::addKmersInGraph(const string& region, uint index){
	uint posi(0);
	uint indexNode(0);
	while (posi + k < region.size()){
		string kmer(region.substr(posi, k));
		if (kmersToNode.count(kmer)){ // complete node information
			(kmersToNode.at(kmer)).addPosition(posi, index);
		} else { // new kmer = new node
			Node newNode(indexNode, kmer, posi, index);
			kmersToNode.insert({kmer, newNode});
			string prefix(nPrefix(k, kmer));
			string suffix(nSuffix(k, kmer));
			//  look for nodes to connect with
			if (prefixes.count(suffix)){
				for (uint i(0); i<prefixes[suffix].size(); ++i){
					newNode.addOutNode(&prefixes[suffix][i]);
					prefixes[suffix][i].addInNode(&newNode);
				}
			}
			if (suffixes.count(prefix)){
				for (uint i(0); i<suffixes[prefix].size(); ++i){
					newNode.addInNode(&suffixes[prefix][i]);
					suffixes[prefix][i].addOutNode(&newNode);
				}
			}
			// add prefix and suffix in tables
			prefixes[prefix].push_back(newNode); //cree si clef n'existe pas, ajoute sinon
			suffixes[suffix].push_back(newNode);
			++indexNode;
		}
		posi += k;
	}
}


void Graph::sequences2dot(const Node& node){
    ofstream out("out.dot",ofstream::out);
    out<<"digraph ham {"<<endl;
    // title
    out << "labelloc=\"t\"" << endl ;
    out << "label = \"k="  << k << "\"" <<  endl;
    out << firstPositions[0]->index << "->" << endl;
    out<<"}"<<endl;
}
// ----- end Class Graph ----- //
 

