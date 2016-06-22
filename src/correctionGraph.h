#include <vector>
#include <string>
#include <iostream>
#include <unordered_map>


#ifndef CORRECGPH
#define CORRECGPH

using namespace std;

typedef unsigned int uint;


// ----- Class Node ----- //

class Node {
 public:
	uint index;
	string kmer;
	uint globalPosition;
	bool visited;
	unordered_map <uint, vector<uint>> readsAndPositions;
	vector <Node*> inNodes;
	vector <Node*> outNodes;
	vector <uint> readsIndex;
	//vector <pair<uint, uint>> readsStrands; //  0 for +, 1 for -

	Node(uint index, const string& kmer, uint position, uint readIndex);
	void addPosition(uint position, uint readIndex);
	//~ void addInNode(Node* nodeIn);
	//~ void addOutNode(Node* nodeOut);
	void duplicateNode();
	void callOutNodes(ofstream* out);
	void getGlobalPosition();
};

// ----- end Class Node ----- //

struct compareNode{
    bool operator()(Node* n1, Node* n2){
        return n1->globalPosition < n2->globalPosition;
    }
};

// ----- Class Graph ----- //
class Graph {
 public:
	Graph(uint k);
	uint k;
	uint size;
	vector <Node*> backbone;
	vector <Node*> firstPositions;
	unordered_map <string, Node*> kmersToNode;
	unordered_map <string, vector<Node*>> prefixes;
	unordered_map <string, vector<Node*>> suffixes;
	void graphClear();
	
	//~ void addFirstNode(Node* node);
	//~ void addKmersInGraph(const string& region, uint index);
	uint createGraphFromSetofRegions(const vector <string>& regionVec);
	void sequences2dot();
	void getStartingNodes();
	void duplicateNode(uint indexNode);
	//~ void callOutNodes(const Node& node);
	void getBackBone(uint bestKmer);
	//~ void traversalBetweenTwoNodes();
	void greedyTraversal();
	void greedyDFStoNextBBNode(vector <Node*>& traversal, vector <Node*>& finalTraversal, Node* nextBBNode);
	uint scoreTraversal(const vector<Node*> traversal);
};
// ----- end Class Graph ----- //

string nPrefix(uint n, const string& sequence);
string nSuffix(uint n, const string& sequence);
bool compareNodebyPosition(const Node* n1, const Node* n2);
 
#endif
