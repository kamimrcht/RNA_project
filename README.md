| **Linux** | **Mac OSX** |
|-----------|-------------|
[![Build Status](https://ci.inria.fr/gatb-core/view/RConnector/job/tool-rconnector-build-debian7-64bits-gcc-4.7/badge/icon)](https://ci.inria.fr/gatb-core/view/RConnector/job/tool-rconnector-build-debian7-64bits-gcc-4.7/) | [![Build Status](https://ci.inria.fr/gatb-core/view/RConnector/job/tool-rconnector-build-macos-10.9.5-gcc-4.2.1/badge/icon)](https://ci.inria.fr/gatb-core/view/RConnector/job/tool-rconnector-build-macos-10.9.5-gcc-4.2.1/)

[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)


# What is Long  Read Connector RNA (LRC_RNA)?
Long read connector RNA  enables to group reads coming from the same gene family in a long reads set *S*. For each read from *S* it provides:
 * For each read of *S*, a  list of reads from *S* that share enough *k*-mers organized in exons
 
# Getting the latest source code

## Requirements

CMake 2.6+; see http://www.cmake.org/cmake/resources/software.html

c++ compiler; compilation was tested with gcc and g++ version>=4.5 (Linux) and clang version>=4.1 (Mac OSX).

## Instructions


    # get a local copy of source code
    git clone --recursive https://github.com/kamimrcht/RNA_project
    
    # compile the code an run a simple test on your computer
    cd RNA_project
    sh INSTALL

# Getting a binary stable release


# Quick start

##Minimal call
Run a simple test looking for reads from data/c1.fasta.gz that can be grouped altogether. Kmers indexed from data/c1.fasta are those occurring at least 2 times. 

	 ./long_read_connector.sh -f data/c1.fasta
     

## Options
	 -p prefix. All out files will start with this prefix. Default="short_read_connector_res"
	 -g: with this option, if a file of solid kmer exists with same prefix name and same k value, then it is re-used and not re-computed.
	 -k value. Set the length of used kmers. Must fit the compiled value. Default=21
	 -F value. Fingerprint size. Size of the key associated to each indexed value, limiting false positives. Default=12
	 -G value. gamma value. MPHF expert users parameter - Default=2
	 -a: kmer abundance min (kmer from bank seen less than this value are not indexed). Default=2
	 -s: minimal number of kmer shared by two reads to be considered as similar. Default=3
	 -t: number of thread used. Default=0


## Output Format
### Short reads counter
Command:
	         ./long_read_connector.sh -f data/c1.fasta

Two first lines of the output file: 

	 #query_read_id mean median min max number of shared 31mers with banq read set data/c2.fasta.gz
	 1 3 2
The first line is the file header. 
The second line can be decomposed as:
   * 1: id of the query read (starts at 1)
   * 3 2: numbers separated by " " indicate the id of other reads from the data-set linked to the query read 
   
## Input read sets
LRC_linker supports fasta file format.
   
#Contact

Contact: Camille Marchet: camille.marchet@irisa.fr
