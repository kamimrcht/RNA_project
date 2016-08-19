#!/usr/bin/env python
import os, sys, random, subprocess
from Bio import SeqIO
from math import exp

#  determine transcripts' expression
def zipfLaw(rank,mostAbund):
	return mostAbund * rank**(-0.7) * exp(rank/10000.0 * (rank/10000.0)**2)

def main():
	if len(sys.argv) == 4:
		gtfFile = sys.argv[1]
		fastaFile = sys.argv[2]
		depth =  int(sys.argv[3])
		readSize = 1000.0
		outFile = open("isoseq.fa", 'w')
		
		transcriptomeLength = 0
		dictTranscripts = dict()
		#~ transcriptsToSequence = dict()
		listTranscript = list()

		#  1. read gtf and simulate expression
		with open(gtfFile) as infile:
			for line in infile:
				string = line.rstrip().split("\t")
				boolExon = True if string[2] == 'exon' else False
				if boolExon:  # for an exon line:
					chromosome = string[0]
					startPosition = int(string[3])
					endPosition = int(string[4])
					geneInfo = string[8].split(";")
					posGeneName = geneInfo[0].find('gene_id')
					posTranscriptName = geneInfo[1].find('transcript_id')
					transcriptomeLength += endPosition - startPosition
					if posGeneName != -1:
						geneName = geneInfo[0][len('gene_id') + 1:][1:-1]
						transcriptName = geneInfo[1][len('transcript_id') + 1:][2:-1]
					if transcriptName in dictTranscripts.keys(): #  add new exon of the transcript
						dictTranscripts[transcriptName]["start"].append(startPosition)
						dictTranscripts[transcriptName]["end"].append(endPosition)
					else: #  add new transcript
						dictTranscripts[transcriptName] = dict()
						dictTranscripts[transcriptName]["start"] = [startPosition]
						dictTranscripts[transcriptName]["end"] = [endPosition]
						dictTranscripts[transcriptName]["gene"] = geneName
						listTranscript.append(transcriptName)
		#~ print transcriptomeLength
		nbMol = transcriptomeLength/readSize * depth
		random.shuffle(listTranscript) #  attribute a rank of expression to each transcript
		for i in range(len(listTranscript)):
			dictTranscripts[listTranscript[i]]["expr"] = zipfLaw(i+1, 1000) * nbMol  # get expression level for each transcript
			#~ print int(round(zipfLaw(i+1, 1000),0))

		#  2. extract sequences of transcripts from ref fasta and simulate sequencing with PBSIM
		fasta = SeqIO.read(fastaFile, "fasta")		
		params = {"file": ""}
		for transcript in dictTranscripts.iterkeys():
			#~ if len(transcriptsToSequence[gene]) > 1:  # alternative isoforms
				#~ for transcript in transcriptsToSequence[gene].iterkeys():
			sequence = ""
			for start, end in zip(dictTranscripts[transcript]["start"], dictTranscripts[transcript]["end"]):
				#  concatenate all exons of a transcript #  not sure if end or end+1
				h= str(fasta.seq[start:end]).upper()
				sequence += h
			outputTemp = open(transcript + ".fa", 'w')
			outputTemp.seek(0)
			outputTemp.write(">gene_" + dictTranscripts[transcript]["gene"] + "|transcript_" + transcript + "\n")
			outputTemp.write(sequence)
			outputTemp.close()
			params["file"] = transcript
			if len(sequence)>100: #   required by PBSIM
				if len(sequence) < 4000:
					print "ccs"
					params["type"] = "ccs"
					params["depth"] = dictTranscripts[transcript]["expr"]
					FNULL = open('/dev/null', 'w')
					bashCommand = ["/home/marchet/bin/pbsim-1.0.3-Linux-amd64-bin/Linux-amd64/bin/pbsim", "%(file)s.fa" %params,  "--model_qc", "/home/marchet/bin/pbsim-1.0.3-Linux-amd64-bin/data/model_qc_clr"]
					proc = subprocess.Popen(bashCommand, stdout=FNULL, stderr=subprocess.PIPE)
					outputErr = proc.stderr.read()
					if not "ERROR" in outputErr:
						print "noerror"
					else:
						print outputErr

		#~ FNULL = open('/dev/null', 'w')
		#~ bashCommand = ["/home/marchet/bin/pbsim-1.0.3-Linux-amd64-bin/Linux-amd64/bin/pbsim", "/home/marchet/RNA_project/B.fasta",  "--model_qc", "/home/marchet/bin/pbsim-1.0.3-Linux-amd64-bin/data/model_qc_clr"]
		#~ proc = subprocess.Popen(bashCommand, stdout=FNULL, stderr=subprocess.PIPE)
		#~ outputErr = proc.stderr.read()
		#~ if not "ERROR" in outputErr:
			#~ print "noerror"
		#~ else:
			#~ print outputErr

if __name__ == '__main__':
    main()
