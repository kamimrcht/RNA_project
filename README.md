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



# Quick start

##Minimal call
Run a simple test looking for reads from data/simulatedReads.fa that can be grouped altogether. Kmers indexed from data/c1.fasta are those occurring at least 2 times. 

	 ./long_read_connector.sh -f data/simulatedReads.fa
     

## Options
	 -p prefix. All out files will start with this prefix. Default="short_read_connector_res"
	 -g: with this option, if a file of solid kmer exists with same prefix name and same k value, then it is re-used and not re-computed.
	 -k value. Set the length of used kmers. Must fit the compiled value. Default=21
	 -F value. Fingerprint size. Size of the key associated to each indexed value, limiting false positives. Default=12
	 -G value. gamma value. MPHF expert users parameter - Default=2
	 -a: kmer abundance min (kmer from bank seen less than this value are not indexed). Default=2
	 -s: minimal number of kmer shared by two reads to be considered as similar. Default=3
	 -t: number of thread used. Default= as much as available


## Output Format
### Short reads counter
Command:
	         ./long_read_connector.sh -f data/simulatedReads.fa

Two first lines of the output file: 
		head -2 long_read_connector_res.txt 
	
561 577 189 906 558 649 358 261 974 964 188 673 479 913 816 355 731 133 327 77 842 560 366 57 251 47 517 3 391 952 758 247 635 750 820 105 207 817 328 93 187 972 27 638 628 971 486 1 277 763 569 475 5 178 836 916 795 810 325 169 555 707 660 49 801 611 417 793 600 741 536 924 505 818 430 54 336 239 297 677 580 156 350
736 957 52 723 145 687 37 425 46 51 342 87 960 171 190 567 191 896 605 291 362 832 738 845 371 131 497 2 311 893 996 705 535 341 41 623 135 426 717 427 136 912 100 946 809 236 624 203 140 819 255 800 744 405 890 770 576 91 867 256 254 545 686 880 298 249
   
	* on each line several reads are grouped because they share enough k-mers in a sufficiently long region. A first group consists of reads from number 561, 577 ... to read 350. The second group contains reads 736... 249.  The numbers correspond to the line of the read in the input file (starts at 1).

## Input read sets
LRC_linker supports fasta file format.
   
#Contact

Contact: Camille Marchet: camille.marchet@irisa.fr
