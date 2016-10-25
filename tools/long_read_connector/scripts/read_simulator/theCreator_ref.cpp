#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iterator>
#include <ctime>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
#include <iterator>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <algorithm>
#include <chrono>
#include <map>
#include <set>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>



using namespace std;



char randomNucleotide(){
	switch (rand() % 4){
		case 0:
			return 'A';
		case 1:
			return 'C';
		case 2:
			return 'G';
		case 3:
			return 'T';
	}
	return 'A';
}



string randomSequence(const uint length){
	string result(length, 'A');
	for(uint i(0); i < length; ++i){
		result[i] = randomNucleotide();
	}	
	return result;
}



//~ string mutateSequence(const string& referenceSequence, uint maxMutRate=6, vector <double> ratioMutation={0.06,0.73,0.21}){
string mutateSequence(const string& referenceSequence, uint maxMutRate=6, vector <double> ratioMutation={0.06,0.73,0.21}){
	string result;
	result.reserve(5 * referenceSequence.size());
	for(uint i(0); i < referenceSequence.size(); ++i){
		uint mutRate(maxMutRate);
		//~ uint mutRate(rand() % maxMutRate);
		double substitutionRate(mutRate * ratioMutation[0]);
		double insertionRate(mutRate * ratioMutation[1]);
		double deletionRate(mutRate * ratioMutation[2]);
		uint dice(rand() % 100);
		if(dice<substitutionRate){
			//SUBSTITUTION
			char newNucleotide(randomNucleotide());
			while(newNucleotide == referenceSequence[i]){
				newNucleotide = randomNucleotide();
			}
			result.push_back(newNucleotide);
			continue;
		}
		if(dice < deletionRate+substitutionRate){
			//DELETION
			
			continue;
		}
		if(dice < deletionRate + substitutionRate + insertionRate){
			//INSERTION
			char newNucleotide(randomNucleotide());
			result.push_back(newNucleotide);
			--i;
			continue;
		}
		//NO ERROR
		result.push_back(referenceSequence[i]);
		
	}
	return result;
}



//~ vector<string> generateAlternativeTranscriptReferences(uint transcriptNumber=3, uint totalExonNumber=15, uint exonNumber=6, uint sizeExons=100){
vector<vector<string>> generateAlternativeTranscriptReferences(ifstream& refFile, uint referenceNumber, uint transcriptNumber=3, uint totalExonNumber=14, uint exonNumber=8, uint sizeExons=100){
	string sequence;
	vector<vector<string>> result;
	uint nbRef(0);
	uint position(0);
	string transcript;
	unordered_set<uint> selectedExons;
	while (not refFile.eof() and nbRef < referenceNumber){
        getline(refFile, sequence);
		getline(refFile, sequence);
		uint lengthExons;
		if (sequence.size()/sizeExons < totalExonNumber){
			lengthExons = sequence.size()/(totalExonNumber+1);
		} else {
			lengthExons = sizeExons;
		}
		if (lengthExons > 30){
			position = 0;
			vector<string> exonList;
			uint i(0);
			while (i < totalExonNumber and position < sequence.size()){ // creation of exons from genomic sequence
				string exon(sequence.substr(position, lengthExons));
				exonList.push_back(exon);
				position += lengthExons;
				++i;
			}
			++nbRef;
			vector<string> transcriptList;
			for(uint i(0); i < transcriptNumber; ++i){
				uint dice1(rand() % (totalExonNumber - 3));
				uint transcriptExonNumber(dice1 + 3);
				transcript = "";
				selectedExons = {};
				while(selectedExons.size() != transcriptExonNumber){
					selectedExons.insert(rand() % exonList.size());					
				}
				for(uint ii(0); ii < exonList.size(); ++ii){
					if(selectedExons.count(ii) == 1){
						transcript += exonList[ii];
					}
				}
				transcriptList.push_back(transcript);
			}
			result.push_back(transcriptList);
		}
	}
	return result;
}


unordered_map<uint, uint> geneExpressionChunks(uint referencesNumber, uint numberReads){
	uint nbMidExpressed(referencesNumber * 0.6); // code 1
	uint nbHighlyExpressed(referencesNumber * 0.1+1); // code 0

	
	
	unordered_map<uint, uint> result;
	cout << "high:" << nbHighlyExpressed << " medium:" << nbMidExpressed << " low:" <<  referencesNumber - nbMidExpressed - nbHighlyExpressed << endl;
	unordered_set<uint> refs;
	for (uint i(0); i < referencesNumber; ++i){
		refs.insert(i);
	}
	uint i(0);
	for (auto j (refs.begin()); j != refs.end(), i < nbHighlyExpressed;){
			result.insert({*j, 0});
			j = refs.erase(j);
			++i;
	}
	i = 0;
	for (auto j (refs.begin()); j != refs.end(), i < nbMidExpressed;){
			result.insert({*j, 1});
			j = refs.erase(j);
			++i;
	}
	cout << "check size " << refs.size();
	for (auto j (refs.begin()); j != refs.end(); ++j){
			
			result.insert({*j, 2});
	}
	return result;
}



void generateReads(uint numberReads, ifstream& inRef, uint referencesNumber=25, const string& outFileName="simulatedReads.fa", const string& outRefFileName="RefFile"){
	ofstream out(outFileName);
	ofstream outRef(outRefFileName);
	vector<vector<string>> referenceList(generateAlternativeTranscriptReferences(inRef, referencesNumber));

	int highly(numberReads*0.5);
	int regular(numberReads*0.4);
	int rare(numberReads - highly - regular);
	
	
	unordered_map<uint, uint> geneExpression(geneExpressionChunks(referencesNumber, numberReads));

	cout << "nb reads:" << numberReads << " nb rare reads:" << rare << " nb regular reads:"<< regular << " nb high reads:"<< highly << endl;
	
	for(uint i(0);i < referenceList.size(); ++i){
		string expr;
		switch (geneExpression[i]){
		case 0:
			expr = "highExpression";
			break;
		case 1:
			expr = "regularExpression";
			break;
		case 2:
			expr = "shallowExpression";
			break;
		}
		for(uint ii(0); ii<referenceList[i].size(); ++ii){
			outRef << ">referenceNumber:" << i << " alternativeNumber" << ii << " " <<  expr <<endl;
			outRef << referenceList[i][ii] << endl;
		}
	}
	
	string refRead,realRead;

	int hi(0), re(0), ra(0);
	while (hi < highly or re < regular or ra < rare){
		uint dice1(rand() % referencesNumber);
		uint dice2(rand() % (referenceList[dice1].size() ));
		int expr(geneExpression[dice1]);
		if (expr == 0 and hi < highly){
			++hi;
			string expression = "highExpression";
			refRead = referenceList[dice1][dice2];
			realRead = mutateSequence(refRead);
			out << ">referenceNumber:" << dice1 << " alternativeNumber:" << dice2  << " " << expression << endl;
			out << realRead << endl;
		}
		if (expr == 1 and re < regular){
			++re;
			string expression = "regularExpression";
			refRead = referenceList[dice1][dice2];
			realRead = mutateSequence(refRead);
			out << ">referenceNumber:" << dice1 << " alternativeNumber:" << dice2  << " " << expression << endl;
			out << realRead << endl;
		}
		if (expr == 2 and ra < rare){
			++ra;
			string expression = "lowExpression";
			refRead = referenceList[dice1][dice2];
			realRead = mutateSequence(refRead);
			out << ">referenceNumber:" << dice1 << " alternativeNumber:" << dice2  << " " << expression << endl;
			out << realRead << endl;
		}
	}
}



int main(int argc, char ** argv){
	if (argc > 1){
		string refName = argv[1];
		ifstream refFile(refName);
		srand (time(NULL));
		auto startChrono = chrono::system_clock::now();
		generateReads(100000, refFile);
		auto end = chrono::system_clock::now(); auto waitedFor = end - startChrono;
		cout << "Time  in ms : " << (chrono::duration_cast<chrono::milliseconds>(waitedFor).count()) << endl;
	}
	return 0;
}
