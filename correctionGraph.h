#include <vector>
#include <string>
#include <iostream>
#include <unordered_map>


#ifndef CORRECGPH
#define CORRECGPH

using namespace std;

typedef unsigned int uint;

string nPrefix(uint n, const string& sequence);
string nSuffix(uint n, const string& sequence);
// ----- Class Node ----- //

class Node {
 public:
	uint index;
	string kmer;
	unordered_map <uint, vector<uint>> readsAndPositions;
	vector <Node*> inNodes;
	vector <Node*> outNodes;
	vector <uint> readsIndex;
	//vector <pair<uint, uint>> readsStrands; //  0 for +, 1 for -

	Node(uint index, const string& kmer, uint position, uint readIndex);
	void addPosition(uint position, uint readIndex);
	void addInNode(Node* nodeIn);
	void addOutNode(Node* nodeOut);
	void duplicateNode();
};

// ----- end Class Node ----- //

// ----- Class Graph ----- //
class Graph {
 public:
	Graph(uint k);
	uint k;
	uint size;
	vector <Node*> firstPositions;
	unordered_map <string, Node> kmersToNode;
	unordered_map <string, vector<Node>> prefixes;
	unordered_map <string, vector<Node>> suffixes;

	
	void addFirstNode(Node* node);
	void addKmersInGraph(const string& region, uint index);
	void sequences2dot(const Node& node);
};
// ----- end Class Graph ----- //


 
#endif
