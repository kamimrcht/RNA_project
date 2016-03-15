#!/usr/bin/env python
import os, sys, random, subprocess
from Bio import SeqIO

gtfFile = sys.argv[1]
fastaFile = sys.argv[2]
outFile = open("isoseq.fa", 'w')
transcriptsToSequence = dict()

previousGene = None
with open(gtfFile) as infile:
	for line in infile:
		string = line.rstrip().split("\t")
		boolExon = True if string[2] == 'exon' else False
		if boolExon:
			chromosome = string[0]
			startPosition = int(string[3])
			endPosition = int(string[4])
			geneInfo = string[8].split(";")
			posGeneName = geneInfo[0].find('gene_id')
			posTranscriptName = geneInfo[1].find('transcript_id')
			if posGeneName != -1:
				geneName = geneInfo[0][len('gene_id') + 1:][1:-1]
				transcriptName = geneInfo[1][len('transcript_id') + 1:][2:-1]
			if previousGene is not None:
				if geneName in transcriptsToSequence.keys():
					if transcriptName in transcriptsToSequence[geneName].keys():
					#  another exon from the transcript: add it to what will be sequenced
						transcriptsToSequence[geneName][transcriptName]["start"].append(startPosition)
						transcriptsToSequence[geneName][transcriptName]["end"].append(endPosition)
					else:
					#  probability to sequence an alternative transcript from the gene
						probaToSeqAlt = random.uniform(0,1)
						if probaToSeqAlt <= 1:
							transcriptsToSequence[geneName][transcriptName] = {"start": [startPosition], "end": [endPosition]}
				else:
					probaToSeqGene = random.uniform(0,1)
					if probaToSeqGene < 0.4:
						previousGene = geneName
						transcriptsToSequence.setdefault(geneName, {transcriptName: {"start": [startPosition], "end": [endPosition]}})
			else:
				probaToSeqGene = random.uniform(0,1)
				if probaToSeqGene < 0.4:
					transcriptsToSequence.setdefault(geneName, {transcriptName: {"start": [startPosition], "end": [endPosition]}})
					previousGene = geneName
					
fasta = SeqIO.read(fastaFile, "fasta")
params = {"file": ""}
FNULL = open(os.devnull, 'w')
for gene in transcriptsToSequence.iterkeys():
	if len(transcriptsToSequence[gene]) > 1:  # alternative isoforms
		for transcript in transcriptsToSequence[gene].iterkeys():
			sequence = ""
			for start, end in zip(transcriptsToSequence[gene][transcript]["start"], transcriptsToSequence[gene][transcript]["end"]):
				#~ print start, end
				sequence += fasta.seq[start:end].upper() #  not sure if end or end+1 #  ? un ficher par seq pour pbsim ? // lancer pbsim dans le script python sur chaque sequence puis effacer le fasta de ref, noter tous les transcripts/genes/positions exons dans un fichier
			outputTemp = open(transcript + ".fa", 'w')
			outputTemp.seek(0)
			outputTemp.write(">gene_" + gene + "|" + str(len(transcriptsToSequence[gene])) + "|transcript_" + transcript + "\n")
			outputTemp.write(str(sequence))
			params["file"] = transcript
			if len(sequence) < 2000:
				params["type"] = "ccs"
				bashCommand = ["/home/marchet/bin/pbsim-1.0.3-Linux-amd64-bin/Linux-amd64/bin/pbsim", "%(file)s.fa" %params, "--prefix", "%(file)s" %params, "--depth", "10",  "--length-mean", "1500",  "--length-max", "2000", "--model_qc", "/home/marchet/bin/pbsim-1.0.3-Linux-amd64-bin/data/model_qc_ccs", "--data-type", "CCS"]
			else:
				params["type"] = "clr"
				bashCommand = ["/home/marchet/bin/pbsim-1.0.3-Linux-amd64-bin/Linux-amd64/bin/pbsim", "%(file)s.fa" %params, "--prefix", "%(file)s" %params, "--depth", "10",  "--length-mean", "5000", "--model_qc", "/home/marchet/bin/pbsim-1.0.3-Linux-amd64-bin/data/model_qc_clr"]
			proc = subprocess.Popen(bashCommand, stdout=FNULL, stderr=subprocess.PIPE)
			outputErr = proc.stderr.read()
			if not "ERROR" in outputErr:
				fastq = open("%(file)s_0001.fastq" %params, 'r')
				linesFq = fastq.readlines()
				for i in range(len(linesFq)):
					if i%4 == 0:
						line = ">gene_" + gene + "|" + str(len(transcriptsToSequence[gene])) + "|transcript_" + transcript + "|"  + "%(type)s_" %params + linesFq[i][1:] 
						outFile.write(line)
					if i%4==1:
						outFile.write(linesFq[i])
			remove = "rm %(file)s.fa" % params
			os.system(remove)
removeRef =  "rm *_0001.ref" % params
removeMaf = "rm *_0001.maf" % params
os.system(removeRef)
os.system(removeMaf)
			
