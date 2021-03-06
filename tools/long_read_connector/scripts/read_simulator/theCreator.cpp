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


void insertion(uint rate, string& result){
	uint dice(rand() % 100);
	if(dice < rate){
		char newNucleotide(randomNucleotide());
		result.push_back(newNucleotide);
		insertion(rate, result);
	}
}


string mutateSequence(const string& referenceSequence, uint maxMutRate, vector <double> ratioMutation={0.06,0.73,0.21}){
//~ string mutateSequence(const string& referenceSequence, uint maxMutRate=3, vector <double> ratioMutation={0.06,0.73,0.21}){
	string result;
	result.reserve(5 * referenceSequence.size());
	for(uint i(0); i < referenceSequence.size(); ++i){
		uint mutRate(maxMutRate);
		//~ uint mutRate(rand() % maxMutRate);
		double substitutionRate(mutRate * ratioMutation[0]);
		double insertionRate(mutRate * ratioMutation[1]);
		double deletionRate(mutRate * ratioMutation[2]);
		uint dice(rand() % 100);


		if (dice <substitutionRate ){
			//SUBSTITUTION
			char newNucleotide(randomNucleotide());
			while(newNucleotide == referenceSequence[i]){
				newNucleotide = randomNucleotide();
			}
			result.push_back(newNucleotide);
			continue;
		} else if(dice < deletionRate+substitutionRate){
			//DELETION
			uint dice2(rand() % 100);
			while (dice2 < deletionRate+substitutionRate){ // deletions larger than 1
				++i;
				dice2 = rand() % 100;
			}
			continue;
		} else if (dice < deletionRate + substitutionRate + insertionRate){
			//INSERTION
			char newNucleotide(randomNucleotide());
			result.push_back(referenceSequence[i]);
			result.push_back(newNucleotide);
			//~ --i;
			insertion(deletionRate + substitutionRate + insertionRate, result); // larger than 1 insertions
			
			continue;
		} else {
		//NO ERROR
			result.push_back(referenceSequence[i]);
		}

	}
	return result;
}



vector<string> generateAlternativeTranscriptReferences(uint transcriptNumber=3, uint totalExonNumber=15, uint exonNumber=12, uint sizeExons=200){

	vector<string> result;
	vector<string> exonList;
	for(uint i(0); i < totalExonNumber; ++i){
		exonList.push_back(randomSequence(sizeExons));
	}
	string transcript;
	transcript.reserve(exonNumber*sizeExons);
	unordered_set<uint> selectedExons;
	uint dice1, transcriptExonNumber;
	for(uint i(0); i < transcriptNumber; ++i){
		if (transcriptNumber > 1){
			dice1 = (rand() % (totalExonNumber-3));
			transcriptExonNumber = (dice1 + 3);
		} else {
			transcriptExonNumber = 0;
		}
		transcript = "";
		selectedExons = {};
		if (transcriptExonNumber > 0){
			while(selectedExons.size() != transcriptExonNumber){
				selectedExons.insert(rand() % totalExonNumber);
			}
		} else {
			selectedExons.insert(0);
		}
		for(uint ii(0); ii < totalExonNumber; ++ii){
			if(selectedExons.count(ii) == 1){
				transcript += exonList[ii];
			}
			
		}
		result.push_back(transcript);
	}
	return result;
}




void generateReads(uint numberReads, uint mutRate, uint referencesNumber, const string& outFileName="simulatedReads.fa", const string& refFileName="RefFile"){
	ofstream out(outFileName);
	ofstream outRef(refFileName);
	vector<vector<string>> referenceList;
	for(uint i(0);i < referencesNumber; ++i){
		referenceList.push_back(generateAlternativeTranscriptReferences());
		for(uint ii(0); ii<referenceList[i].size(); ++ii){
			outRef << ">referenceNumber:" << i << " alternativeNumber" << ii << " length"<< referenceList[i][ii].size() << endl;
			outRef << referenceList[i][ii] << endl;
		}
	}

	string refRead,realRead;
	for(uint i(0); i < numberReads; ++i){
		uint dice1(rand() % referencesNumber);
		uint dice2(rand() % referenceList[dice1].size());
		refRead = referenceList[dice1][dice2];
		realRead = mutateSequence(refRead, mutRate);
		out << ">" << i + 1 << "_referenceNumber:" << dice1 << " alternativeNumber" << dice2 << " length" << realRead.size() << endl;
		out << realRead << endl;
	}
}



int main(int argc, char ** argv){
	srand (time(NULL));
	if (argc > 3){
		auto startChrono = chrono::system_clock::now();
		uint nbReads(stoi(argv[1]));
		uint mut(stoi(argv[2]));
		uint genes(stoi(argv[3]));
		generateReads(nbReads, mut, genes);
		auto end = chrono::system_clock::now(); auto waitedFor = end - startChrono;
		cout << "Time  in ms : " << (chrono::duration_cast<chrono::milliseconds>(waitedFor).count()) << endl;
	}
	return 0;
}
