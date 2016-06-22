#include <vector>
#include <string>
#include "correctionGraph.h"
#include <iostream>
#include <fstream>
#include <algorithm>



using namespace std;

// ----- Class Node ----- //

// Node constructor
Node::Node(uint indexI, const string& kmerI, uint position, uint readIndex){
	index = indexI;
	kmer = kmerI;
	readsAndPositions[readIndex].push_back(position);
	visited = false;
}

// add a new position (and read) of a k-mer in the node
void Node::addPosition(uint position, uint readIndex){
	readsAndPositions[readIndex].push_back(position);
}




//~ // add a new in-node to the node's neighbours (node )
//~ void Node::addInNode(Node* nodeIn){
	//~ inNodes.push_back(nodeIn);
//~ }

//~ // add a new out-node to the node's neighbours
//~ void Node::addOutNode(Node* nodeOut){
	//~ outNodes.push_back(nodeOut);
//~ }

void Node::callOutNodes(ofstream* out){
	for (uint i(0); i < outNodes.size(); ++i){
		if (outNodes[i]->visited == false){
			*(out) << kmer  << readsAndPositions.size() << "->" << outNodes[i]->kmer << outNodes[i]->readsAndPositions.size() << endl;
			//~ cout << kmer << "->" << outNodes[i]->kmer << endl;
			outNodes[i]->callOutNodes(out);
			outNodes[i]->visited = true;
		} else {
			*(out) << kmer << readsAndPositions.size() << "->" << outNodes[i]->kmer << outNodes[i]->readsAndPositions.size()  << endl;
		}
	}
}

void Node::getGlobalPosition(){
	uint avg(0);
	uint index(0);
	bool repeated(false);
	for (auto i(readsAndPositions.begin()); i != readsAndPositions.end(); ++i){
		if (i->second.size()>1){
			repeated = true;
		}
		if (not repeated){
			avg += i->second[0];
			++ index;
		}
	}
	if (not repeated){
		globalPosition = avg/index;
	} else {
		globalPosition = 0;
	}
}
// ----- end Class Node ----- //




// ----- Class Graph ----- //


// graph constructor
Graph::Graph(uint kmerSize){
	k = kmerSize;
}



//construct a graph with all fragments of reads corresponding to a region
//~ void Graph::addKmersInGraph(const string& region, uint index){
uint  Graph::createGraphFromSetofRegions(const vector <string>& regionVec){
	uint kmerSupport(1);
	for (uint index(0); index < regionVec.size(); ++index){
		uint posi(0);
		uint indexNode(0);
		string region(regionVec[index]);
		//~ cout << kmerSupport;
		while (posi + k <= region.size()){
			string kmer(region.substr(posi, k));
			if (kmersToNode.count(kmer)){ // complete node information
				(kmersToNode.at(kmer))->addPosition(posi, index);
				//~ cout << kmersToNode.at(kmer)->readsAndPositions.size() << endl;
				if (kmersToNode.at(kmer)->readsAndPositions.size() > kmerSupport){
					kmerSupport = kmersToNode.at(kmer)->readsAndPositions.size();
				}
			} else { // new kmer = new node
				Node* newNode = new Node(indexNode, kmer, posi, index);
				string prefix(nPrefix(k-1, kmer));
				string suffix(nSuffix(k-1, kmer));
				//  look for nodes to connect with
				if (prefixes.count(suffix)){
					for (uint j(0); j < prefixes[suffix].size(); ++j){
						newNode->outNodes.push_back(prefixes[suffix][j]);
						prefixes[suffix][j]->inNodes.push_back(newNode);
						//~ out << newNode.kmer<< "->" << prefixes[suffix][j].kmer << endl;
					}
				}
				if (suffixes.count(prefix)){
					for (uint l(0); l < suffixes[prefix].size(); ++l){
						newNode->inNodes.push_back(suffixes[prefix][l]);
						suffixes[prefix][l]->outNodes.push_back(newNode);
						//~ out << suffixes[prefix][l].kmer << "->" << newNode.kmer << endl;
					}
				}
				// add prefix and suffix in tables
				prefixes[prefix].push_back(newNode); //cree si clef n'existe pas, ajoute sinon
				suffixes[suffix].push_back(newNode);
				kmersToNode.insert({kmer, newNode});
				++indexNode;
			}
			posi += 1;
		}
	}
	return kmerSupport;
	//~ out<<"}"<<endl;
}


void Graph::duplicateNode(uint indexNode){
	
	for (auto i(kmersToNode.begin()); i !=kmersToNode.end(); ++i){
		vector <vector<uint>> vecPositions;
		vector <vector<uint>> vecReadIndex;
		uint indexVec(0);
		//~ for (uint read(0); read<i->second->readsAndPositions.size() ; ++i){
			//~ cout << i->second->readsAndPositions[read][0];
			//~ if(i->second->readsAndPositions[read][1].size() == 1){ // if the kmer is duplicated inside the read itself we dont treat it
				//~ if (read == 0){
					//~ vecPositions[0].push_back(i->second->readsAndPositions[read][1][0]); //readsAndPositions[i] is vector [indexread, [positions in read]]
					//~ vecReadIndex[0].push_back(i->second->readsAndPositions[read][0]);
				//~ } else {
					//~ if (vecPositions[indexVec][vecPositions[indexVec].size()-1] + (uint)5 < i->second->readsAndPositions[read][0]){
						//~ vecReadIndex[indexVec].push_back(i->second->readsAndPositions[read][0]); // read index corresponds to position
					//~ } else {
						//~ ++indexVec;
						//~ vecPositions[indexVec].push_back(i->second->readsAndPositions[read][1][0]);
						//~ vecReadIndex[indexVec].push_back(i->second->readsAndPositions[read][0]);
					//~ }
				//~ }
			//~ }
		//~ }
		//~ if (vecPositions.size() > 1){
			//~ vector <uint> means;
			//~ for (uint val(0); val < vecPositions.size(); ++val){
				//~ uint m(0);
				//~ for (uint j(0); j < vecPositions[val].size(); ++j){
					//~ m += vecPositions[val][j];
				//~ }
				//~ m /= vecPositions[val].size();
				//~ means.push_back(m);
				//~ if (val>0){
					//~ Node* newNode = new Node(indexNode, i->first, m, vecReadIndex[val]); //todo: store readIndex
				//~ }
			//~ }
		//~ }
	}
}


void Graph::graphClear(){
	for (auto i(kmersToNode.begin()); i!=kmersToNode.end(); ++i){
		delete(i->second);
	}
}


void Graph::getStartingNodes(){
	for (auto i(kmersToNode.begin()); i !=kmersToNode.end(); ++i){
		if (i->second->inNodes.empty()){
			firstPositions.push_back(i->second);
		}
	}
}


void Graph::sequences2dot(){
    ofstream out("out.dot",ofstream::out);
    out<<"digraph G {"<<endl;
    // title
    out << "labelloc=\"t\"" << endl ;
    out << "label = \"k="  << k << "\"" <<  endl;
    for (uint i(0); i < firstPositions.size(); ++i){
		firstPositions[i]->callOutNodes(&out);
	}
    out<<"}"<<endl;
}


void Graph::getBackBone(uint bestKmer){
	for (auto i(kmersToNode.begin()); i!=kmersToNode.end(); ++i){
		if (i->second->readsAndPositions.size() == bestKmer){
			bool repeated(false);
			for (auto j(i->second->readsAndPositions.begin()); j!= i->second->readsAndPositions.end(); ++j){// check if the kmer is not repeated
				if (j->second.size() > 1){
					repeated = true;
					break;
				}
			}
			if (not repeated){
				i->second->getGlobalPosition();
				backbone.push_back(i->second);
			}
		}
	}
	sort(backbone.begin(), backbone.end(), compareNode());
	//~ for (uint i(0); i < backbone.size(); ++i){
		//~ cout << backbone[i]->kmer << endl;
	//~ }
}


//~ void Graph::traversalBetweenTwoNodes(){
	//~ vector<Node*> traversal;
	//~ for (uint i(0); i < backbone.size()-1; ++i){
		//~ cout << "bb " << backbone[i]->kmer << " " << backbone[i]->globalPosition << endl;
		//~ uint posi(backbone[i]->globalPosition);
		//~ uint nextPosi(backbone[i+1]->globalPosition);
		//~ uint mid(posi + nextPosi/2);
		//~ cout << mid << endl;
		//~ Node* currentNode(backbone[i]);
		//~ for (uint position(posi); position < nextPosi; ++position){
			//~ traversal.push_back(currentNode);
			//~ for (uint next(0); next < currentNode->outNodes.size(); ++next){
				//~ currentNode->outNodes[next]->getGlobalPosition();
				//~ if (currentNode->outNodes[next] == backbone[i+1]){
					//~ currentNode = currentNode->outNodes[next];
					//~ break;
				//~ }
				//~ if (currentNode->outNodes[next]->globalPosition <= mid){
					//~ cout << currentNode->outNodes[next]->kmer << endl;
					//~ currentNode = currentNode->outNodes[next];
				//~ }
			//~ }
		//~ }
	//~ }
	//~ cout << "*******" << endl;
	//~ for (uint i(0); i < traversal.size(); ++i){
		//~ cout << traversal[i]->kmer << endl;
	//~ }
//~ }

void Graph::greedyTraversal(){
	for (uint index1(0),index2(1); index1< backbone.size()-1 and index2 < backbone.size() ; ++index1, ++index2){
		vector<Node*> traversal;
		traversal.push_back(backbone[index1]);
		cout << "*****" << backbone[index1]->kmer <<  endl;
		vector <Node*> finalTraversal;
		greedyDFStoNextBBNode(traversal, finalTraversal, backbone[index2]);
		//~ cout << finalTraversal.size() << endl;
		for (uint i(0); i < finalTraversal.size(); ++i){
			//~ cout << "p" << endl;
			cout << finalTraversal[i]->kmer << endl;
		}
	}
	
}


void Graph::greedyDFStoNextBBNode(vector <Node*>& traversal, vector <Node*>& finalTraversal, Node* nextBBNode){
	//~ uint finalScore(0);
	Node* currentNode(traversal[traversal.size()-1]);
	//~ vector <Node*> returnedTraversal;
	//~ cout << finalScore << endl;
	if (not currentNode->outNodes.empty()){
		for (uint i(0); i < currentNode->outNodes.size(); ++i){
			//~ vector<Node*> trav={currentNode};
			//~ cout << currentNode->kmer << endl;
			Node* nextNode(currentNode->outNodes[i]);
			if (nextNode == nextBBNode){ // we reached the next node of the backbone
				//~ break;
				uint score(scoreTraversal(traversal));
				uint bestScore(0);
				if (not finalTraversal.empty()){
					bestScore = scoreTraversal(finalTraversal);
				}
				//~ cout << score << " "<< finalScore << endl;
				//~ cout << score << endl;
				if (score >= bestScore){
					//~ cout << "ok" << endl;
					//~ finalScore = score;
					//~ cout << score << " " << bestScore << endl;
					finalTraversal = traversal;
				}
				traversal = {traversal[0]};
			} else {
				nextNode->getGlobalPosition();
				if (nextNode->globalPosition < nextBBNode->globalPosition and nextNode->globalPosition != 0){ // we have not gone over the position of the next bb node yet
					traversal.push_back(currentNode->outNodes[i]);
					greedyDFStoNextBBNode(traversal, finalTraversal, nextBBNode);
					//~ break;
				} else {
					traversal = {}; // if we have not reached
				}
			}
		}
	}
}


uint Graph::scoreTraversal(const vector<Node*> traversal){
	uint score(0);
	for (uint i(0); i<traversal.size(); ++i){
		//~ cout << "." << traversal[i]->kmer << endl;
		score += traversal[i]->readsAndPositions.size();
	}
	return score;
}

// ----- end Class Graph ----- //



string nPrefix(uint n, const string& sequence){
	return sequence.substr(0, n);
}


string nSuffix(uint n,const string& sequence){
	return sequence.substr(sequence.size()-n, n);
}

//~ bool compareNodebyPosition(const Node* n1, const Node* n2){
    //~ return n1->globalPosition < n2->globalPosition;
//~ }


 

