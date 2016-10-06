#include <LRC_linker_rna.hpp>
#include "LRC_cluster_rna.hpp"
#include "unordered_set"
#include "mutex"

using namespace std;


/********************************************************************************/


// We define some constant strings for names of command line parameters
static const char* STR_URI_BANK_INPUT = "-bank";
//~ static const char* STR_URI_QUERY_INPUT = "-query";
static const char* STR_FINGERPRINT = "-fingerprint_size";
static const char* STR_GAMMA = "-gamma";
static const char* STR_THRESHOLD = "-kmer_threshold";
static const char* STR_WINDOW = "-window_size";
static const char* STR_OUT_FILE = "-out";
static const char* STR_CORE = "-core";
static const char* STR_SMALLK = "-small_k";
static const char* STR_NBSMALLK = "-nb_small_k";
static const char* STR_SZWINDOWS = "-size_windows";



uint countRm(0), prevCount(0);
mutex mutex1,mutex2,mutex3,mutex4, mutex5;


LRC_linker_rna::LRC_linker_rna()  : Tool ("SRC_linker_rna"){
	// We add some custom arguments for command line interface
	getParser()->push_back (new OptionOneParam (STR_URI_GRAPH, "graph input",   true));
	getParser()->push_back (new OptionOneParam (STR_URI_BANK_INPUT, "bank input",    true));
	getParser()->push_back (new OptionOneParam (STR_OUT_FILE, "output_file",    true));
	getParser()->push_back (new OptionOneParam (STR_THRESHOLD, "Minimal percentage of shared kmer in a region for considering 2 reads in a same group.",    false, "75"));
	getParser()->push_back (new OptionOneParam (STR_WINDOW, "Size of a region (putative exon).",    false, "80"));
	getParser()->push_back (new OptionOneParam (STR_GAMMA, "gamma value",    false, "2"));
	getParser()->push_back (new OptionOneParam (STR_FINGERPRINT, "fingerprint size",    false, "32"));
	getParser()->push_back (new OptionOneParam (STR_CORE, "Number of thread",    false, "1"));


	getParser()->push_back (new OptionOneParam (STR_SMALLK, "small k value",    false, "9"));
	getParser()->push_back (new OptionOneParam (STR_NBSMALLK, "nb of small k mers to look for",    false, "120"));
	getParser()->push_back (new OptionOneParam (STR_SZWINDOWS, "size of windows",    false, "15"));
	//~ nbRead=0;
}



void LRC_linker_rna::create_quasi_dictionary(int fingerprint_size, int nbCores){
	const int display = getInput()->getInt (STR_VERBOSE);
	// We get a handle on the HDF5 storage object.
	// Note that we use an auto pointer since the StorageFactory dynamically allocates an instance
	auto_ptr<Storage> storage (StorageFactory(STORAGE_HDF5).load (getInput()->getStr(STR_URI_GRAPH)));
	// We get the group for dsk
	Group& dskGroup = storage->getGroup("dsk");
	kmer_size = atoi(dskGroup.getProperty("kmer_size").c_str());
	// We get the solid kmers collection 1) from the 'dsk' group  2) from the 'solid' collection
	Partition<Kmer<>::Count>& solidKmers = dskGroup.getPartition<Kmer<>::Count> ("solid");
	nbSolidKmers = solidKmers.getNbItems();
	if(nbSolidKmers==0){
		cout<<"No solid kmers in bank -- exit"<<endl;
		exit(0);
	}
	IteratorKmerH5Wrapper iteratorOnKmers (solidKmers.iterator());
	quasiDico = quasidictionaryVectorKeyGeneric<IteratorKmerH5Wrapper, u_int32_t> (nbSolidKmers, iteratorOnKmers, fingerprint_size, nbCores, gamma_value);
}



struct FunctorIndexer{
	quasidictionaryVectorKeyGeneric <IteratorKmerH5Wrapper, u_int32_t > &quasiDico;
	int kmer_size;
	ISynchronizer* synchro;
	vector<string>& vec;

	FunctorIndexer(quasidictionaryVectorKeyGeneric <IteratorKmerH5Wrapper, u_int32_t >& quasiDico, int kmer_size, vector<string>& vec, ISynchronizer* synchro)  :  quasiDico(quasiDico), kmer_size(kmer_size), synchro(synchro), vec(vec) {
	}

	void operator() (Sequence& seq){
		if(not valid_sequence(seq,kmer_size)){return;}
		//~ nbRead++;
		Kmer<KMER_SPAN(1)>::ModelCanonical model (kmer_size);
		Kmer<KMER_SPAN(1)>::ModelCanonical::Iterator itKmer (model);
		itKmer.setData (seq.getData());
		u_int32_t read_id = static_cast<u_int32_t>(seq.getIndex()+1);
		for (itKmer.first(); !itKmer.isDone(); itKmer.next()){
			// Adding the read id to the list of ids associated to this kmer.note that the kmer may not exist in the dictionary if it was under the solidity threshold.in this case, nothing is done
			quasiDico.set_value((itKmer)->value().getVal(), read_id);
		}
		synchro->lock();
		uint indexWanted(seq.getIndex() + 1 ); // same index than the qdict, starting at 1
		if (indexWanted>= vec.size()){
		    vec.resize(indexWanted + 2);
		    vec[indexWanted] = seq.toString();
		} else {
		    vec[indexWanted] = seq.toString();
		}
		synchro->unlock();
	}
};



uint nuc2int(char n){
    if(n=='A'){
	return 0;
    }else if(n=='C'){
	return 1;
    }else if(n=='G'){
	return 2;
    }else{
	return 3;
    }
}



uint64_t string2int(const string& str){
    uint64_t res(0);
    for(uint i(0);i<str.size();++i){
	res<<=2;
	res+=nuc2int(str[i]);
    }
    return res;
}



void LRC_linker_rna::fill_quasi_dictionary(const int nbCores, const string& bankName, vector <string>& v){
	bool exists;
	IBank* bank = Bank::open (bankName);
	cout<<"Index "<<kmer_size<<"-mers from bank "<<getInput()->getStr(STR_URI_BANK_INPUT)<<endl;
	LOCAL (bank);
	ProgressIterator<Sequence> itSeq (*bank);
	Dispatcher dispatcher (nbCores, 10000);
	ISynchronizer* synchro = System::thread().newSynchronizer();
	dispatcher.iterate (itSeq, FunctorIndexer(quasiDico, kmer_size, v, synchro));
	//~ cout << "SIZE" <<v.size() << endl;
}



class FunctorQueryMatchingRegions{
public:
	ISynchronizer* synchro;
	FILE* outFile;
	int kmer_size;
	quasidictionaryVectorKeyGeneric <IteratorKmerH5Wrapper, u_int32_t>* quasiDico;
	std::unordered_map<uint64_t, vector<uint32_t>> reads_sharing_kmer_2_positions;  // store the position where a k-mer is seen in a read that can be potentially recruited
	std::unordered_map<uint64_t, vector<readGrouped>> read_group;  // for a read, get all reads sharing at least a window
	uint threshold;
	uint size_window;
	vector<uint32_t> associated_read_ids;
	Kmer<KMER_SPAN(1)>::ModelCanonical model;
	Kmer<KMER_SPAN(1)>::ModelCanonical::Iterator* itKmer;
	vector<string>* vecReads;
	unordered_map<uint64_t, vector<int>>* readRedundancy;
	int smallKsize;
	int nbSmallK;
	int nbWindows;


	FunctorQueryMatchingRegions(const FunctorQueryMatchingRegions& lol) // used by the dispatcher
	{
		size_window=lol.size_window;
		synchro=lol.synchro;
		outFile=lol.outFile;
		kmer_size=lol.kmer_size;
		quasiDico=lol.quasiDico;
		threshold=lol.threshold;
		associated_read_ids=lol.associated_read_ids;
		reads_sharing_kmer_2_positions = lol.reads_sharing_kmer_2_positions;
		//~ reads_sharing_kmer_2_positions.max_load_factor ( 0.5);
		//~ reads_sharing_kmer_2_positions.reserve ( 10000);
		read_group =  lol.read_group;
		model=lol.model;
		itKmer = new Kmer<KMER_SPAN(1)>::ModelCanonical::Iterator (model);
		vecReads = lol.vecReads;
		readRedundancy = lol.readRedundancy;
		smallKsize = lol.smallKsize;
		nbSmallK = lol.nbSmallK;
		nbWindows = lol.nbWindows;
	}


	FunctorQueryMatchingRegions(ISynchronizer* synchro, FILE* outFile,  const int kmer_size,  quasidictionaryVectorKeyGeneric <IteratorKmerH5Wrapper, u_int32_t >* quasiDico, const int threshold, const uint size_window, std::unordered_map<uint64_t, vector<uint>>& reads_sharing_kmer_2_positions, std::unordered_map<uint64_t, vector<readGrouped>> read_group, vector<string>* vecReads, unordered_map<uint64_t, vector<int>>* readRedundancy, int smallKsize, int nbSmallK, int nbWindows)
	: synchro(synchro), outFile(outFile), kmer_size(kmer_size), quasiDico(quasiDico), threshold(threshold),  size_window(size_window), reads_sharing_kmer_2_positions(reads_sharing_kmer_2_positions), read_group(read_group), vecReads(vecReads), readRedundancy(readRedundancy), smallKsize(smallKsize), nbSmallK(nbSmallK), nbWindows(nbWindows){
		model=Kmer<KMER_SPAN(1)>::ModelCanonical (kmer_size);
	}


	FunctorQueryMatchingRegions(){}
	void operator() (Sequence& seq){
		//~ return;
	    if (not valid_sequence(seq, kmer_size)){return;} // query sequence
		uint64_t seqIndex(seq.getIndex() + 1);
		if (not (*vecReads)[seqIndex].empty()){
		    bool exists;
		    associated_read_ids = {}; // list of the ids of reads from the bank where a kmer occurs
		    reads_sharing_kmer_2_positions = {};  // store the position where a k-mer is seen in a read that can be potentially recruited
		    itKmer->setData(seq.getData());
		    uint i(0); // position on the read
		    for (itKmer->first(); !itKmer->isDone(); itKmer->next()){
				quasiDico->get_value((*itKmer)->value().getVal(), exists, associated_read_ids); // warning: a same read id can be stored several times - indexes of reads sharing kmers with query are stored in associated_read_ids
				if(!exists) {++i; continue;}
				for (uint r(0); r < associated_read_ids.size(); ++r){ // map: associated_read to vector of kmers
					if (reads_sharing_kmer_2_positions.count(associated_read_ids[r])){
						if (i > reads_sharing_kmer_2_positions[associated_read_ids[r]][reads_sharing_kmer_2_positions[associated_read_ids[r]].size() - 1]) { /*because of the warning above */
							reads_sharing_kmer_2_positions[associated_read_ids[r]].push_back(i);
						}
					} else {
						if (associated_read_ids[r] != seq.getIndex() + 1){  // we dont want to store the info about a read similar to itself
							reads_sharing_kmer_2_positions[associated_read_ids[r]].push_back(i);
						}
					}
				}
				++i;
		    }
		    uint lenSeq((uint)(*vecReads)[seqIndex].size());
		    vector<uint> toErase;
		    bool found(false);
		    for ( auto it = reads_sharing_kmer_2_positions.begin(); it != reads_sharing_kmer_2_positions.end(); ++it ){
			found = false;
			int element((int)it->first);
			if(it->second.size() < 10){//TODO 3 is arbitrary threshold
			    toErase.push_back(it->first);
			}
			if (it->second.size() > (lenSeq - kmer_size + 1)*0.05 and  (*vecReads)[element].size()*0.95 <= lenSeq and  lenSeq <= (*vecReads)[element].size()*1.05 ){ // todo bad threshold
			    found = true;
			    for (int toEraseFromQD(0); toEraseFromQD < (int)(*vecReads)[element].size() - kmer_size + 1; ++toEraseFromQD){
				    string kmer((*vecReads)[element].substr(toEraseFromQD, kmer_size));
				    uint64_t kmerInt(string2int(kmer));
				    uint32_t elem(it->first);
				    mutex1.lock();
				    quasiDico->remove(kmerInt, elem, countRm);  
				    mutex1.unlock();
				}
				mutex2.lock();
				(*vecReads)[element] = "";
				mutex2.unlock();
			}
			if (found){
				bool confirm(false);
				mutex5.lock();
				if (read_group.count(seqIndex)){
				    read_group[seqIndex].push_back({element, confirm});
				} else {
				    readGrouped rg({element, confirm});
				    vector <readGrouped> v({rg});
				    read_group[seqIndex] = {v};
				}
				mutex5.unlock();
			}
			
		    }
		    for(uint locali(0); locali < toErase.size(); ++locali){
			reads_sharing_kmer_2_positions.erase(locali);
		    }
		    vector <uint> matrix(uint(lenSeq) * reads_sharing_kmer_2_positions.size(), 0);
		    uint indexJMatrix(0);
		    for (auto r(reads_sharing_kmer_2_positions.begin()); r != reads_sharing_kmer_2_positions.end(); ++r){
			found = false;
			int element((int)r->first);
			if (not (*vecReads)[element].empty()){
			    unordered_set <string> kmers_read_shared;
			    for (uint v(0); v < (*vecReads)[element].size() - smallKsize + 1; ++v){
				string kmer((*vecReads)[element].substr(v, smallKsize));
				kmers_read_shared.insert(kmer);
			    }
			    for (uint vv(0); vv < lenSeq - smallKsize + 1; ++vv){
				string smallKmer((*vecReads)[seqIndex].substr(vv, smallKsize)); // get all small k-mers from the query sequence
				if (kmers_read_shared.count(smallKmer)){
				    matrix[vv + lenSeq * indexJMatrix] = 1;
				}
			    }
			    uint start(0);
			    uint count(0);
			    uint startKmerPosi(0);
			    uint endKmerPosi(0);
			    uint nbWindowsHit(0);
			    for (uint w(0); w < lenSeq - smallKsize + 1; ++w){
				
				if (w < nbWindows){
				    if (matrix[w + lenSeq * indexJMatrix] == 1){
					endKmerPosi = w;
					++count;
				    }
			       
				} else {
				    start = w - nbWindows + 1;
				    endKmerPosi = w;
				    startKmerPosi = start;
				    if (matrix[start - 1 + lenSeq * indexJMatrix] == 1 and count > 0){
					--count;
				    }
				    if (matrix[w  + lenSeq * indexJMatrix] == 1){
					++count;
				    }
				}
				if (count > nbSmallK){
				    found = true;
				    ++nbWindowsHit;
				    break;
				}
			    }
			    if (found){
				    bool confirm(false);
				    mutex4.lock();
				    if (read_group.count(seqIndex)){
					read_group[seqIndex].push_back({element, confirm});
				    } else {
					readGrouped rg({element, confirm});
					vector <readGrouped> v({rg});
					read_group[seqIndex] = {v};
				    }
				    mutex4.unlock();
			    }
			    ++indexJMatrix; // one j per read recruited
			}
		    }
		    string toPrint;
		    bool read_id_printed = false; // Print (and sync file) only if the read is similar to something.
		    if (not read_group[seqIndex].empty()){
			    if (not read_id_printed){
			    read_id_printed = true;
			    toPrint = to_string(seqIndex) + ":";
			    }
			    for (uint i(0); i < read_group[seqIndex].size(); ++i){
				    toPrint += to_string(read_group[seqIndex][i].index) + " ";
			    }
			    if (read_id_printed){
			    mutex3.lock();
			    toPrint += "\n";
			    fwrite(toPrint.c_str(), sizeof(char), toPrint.size(), outFile);
			    mutex3.unlock();
			    }
		    }
		    if (countRm > prevCount){
			    prevCount = countRm;
		    }
	    }
	}
};



void LRC_linker_rna::parse_query_sequences(int threshold, uint size_window, const int nbCores, const string& bankName, vector <string>* vecReads, unordered_map<uint64_t, vector<int>>* readRedundancy, int smallK, int nbSmallK, int nbWindows){
    IBank* bank = Bank::open(bankName);
    cout<<"Query "<<kmer_size<<"-mers from bank "<< bankName <<endl;
    FILE * outFile;
    outFile = fopen ("long_read_connector_res.tmp", "wb");
    //~ outFile = fopen (getInput()->getStr(STR_OUT_FILE).c_str(), "wb");
    string message("#query_read_id [target_read_id-kmer_span (k="+to_string(kmer_size)+")-kmer_span query percentage]* or U (unvalid read, containing not only ACGT characters or low complexity read)\n");
    fwrite((message).c_str(), sizeof(char), message.size(), outFile);
    LOCAL (bank);
    ProgressIterator<Sequence> itSeq (*bank);
    ISynchronizer* synchro = System::thread().newSynchronizer();
    Dispatcher dispatcher (nbCores, 1000);
    std::unordered_map<uint64_t, vector<uint>> reads_sharing_kmer_2_positions;
    std::unordered_map<uint64_t, vector<readGrouped>> read_group;
    dispatcher.iterate(itSeq, FunctorQueryMatchingRegions(synchro, outFile, kmer_size, &quasiDico, threshold, size_window, reads_sharing_kmer_2_positions, read_group, vecReads, readRedundancy, smallK, nbSmallK, nbWindows));
    fclose (outFile);
    delete synchro;
}



void LRC_linker_rna::execute(){
	int nbCores = getInput()->getInt(STR_CORE);
	int fingerprint_size = getInput()->getInt(STR_FINGERPRINT);
	gamma_value = getInput()->getInt(STR_GAMMA);

	//~ int small_k =  9;
	//~ int nb_small_k = 100;
	//~ int size_windows = 41;
	int small_k = getInput()->getInt(STR_SMALLK);
	int nb_small_k = getInput()->getInt(STR_NBSMALLK);
	int size_windows = getInput()->getInt(STR_SZWINDOWS);
	// IMPORTANT NOTE:
	// Actually, during the filling of the dictionary values, one may fall on non solid non indexed kmers
	// that are quasi dictionary false positives (ven with a non null fingerprint. This means that one nevers knows in advance how much
	// values are gonna be stored for all kmers. This is why I currently us a vector<u_int32_t> for storing read ids associated to a kmer.

	// We need a non null finger print because of non solid non indexed kmers
	//	if (getInput()->getStr(STR_URI_BANK_INPUT).compare(getInput()->getStr(STR_URI_QUERY_INPUT))==0)
	//		fingerprint_size=0;
	vector <string> readsVector;
	cout<<"fingerprint = "<<fingerprint_size<<endl;
	create_quasi_dictionary(fingerprint_size, nbCores);
	string bankName(getInput()->getStr(STR_URI_BANK_INPUT));
	fill_quasi_dictionary(nbCores, bankName, readsVector);
	int threshold = getInput()->getInt(STR_THRESHOLD);
	uint size_window =  getInput()->getInt(STR_WINDOW);
	unordered_map <uint64_t, vector<int>> readRedundancy;
	parse_query_sequences(threshold, size_window, nbCores, bankName, &readsVector, &readRedundancy, small_k, nb_small_k, size_windows);
	LRC_cluster_rna clust("long_read_connector_res.tmp", getInput()->getStr(STR_OUT_FILE).c_str());
	clust.execute();


	getInfo()->add (1, &LibraryInfo::getInfo());
	getInfo()->add (1, "input");
	getInfo()->add (2, "Sequences bank",  "%s",  getInput()->getStr(STR_URI_BANK_INPUT).c_str());
	//~ getInfo()->add (2, "Query bank",  "%s",  getInput()->getStr(STR_URI_BANK_INPUT).c_str());
	getInfo()->add (2, "Kmer size",  "%d",  kmer_size);
	getInfo()->add (2, "Fingerprint size",  "%d",  fingerprint_size);
	getInfo()->add (2, "gamma",  "%d",  gamma_value);
	getInfo()->add (2, "Minimal kmer span percentage",  "%d",  threshold);
	getInfo()->add (1, "output");
	getInfo()->add (2, "Results written in",  "%s",  getInput()->getStr(STR_OUT_FILE).c_str());

	getInfo()->add (2, "Small k size",  "%d",  small_k);
	getInfo()->add (2, "Number of small k-mers",  "%d",  nb_small_k);
	getInfo()->add (2, "Size of windows",  "%d",  size_windows);

}
