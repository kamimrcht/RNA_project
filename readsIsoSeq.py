#!/usr/bin/env python
import os, sys, random, subprocess
from Bio import SeqIO
from math import exp

#  determine transcripts' expression
def zipfLaw(rank, mostAbund):
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
		#~ nbMol = transcriptomeLength/(readSize * depth)
		nbMol = len(dictTranscripts) * 10
		#~ print nbMol
		random.shuffle(listTranscript) #  attribute a rank of expression to each transcript
		for i in range(len(listTranscript)):
			dictTranscripts[listTranscript[i]]["expr"] = zipfLaw(i+1, 100) * nbMol  # get expression level for each transcript
			print "expr", round(zipfLaw(i+1, nbMol/(depth*1.0)) * nbMol, 1) * nbMol
			print "len T", transcriptomeLength
			#~ print int(round(zipfLaw(i+1, 1000),0))

		#  2. extract sequences of transcripts from ref fasta and simulate sequencing with PBSIM
		fasta = SeqIO.read(fastaFile, "fasta")

		params = {"file": ""}
		FNULL = open('/dev/null', 'w')
		for transcript in dictTranscripts.iterkeys():
			#~ if len(transcriptsToSequence[gene]) > 1:  # alternative isoforms
				#~ for transcript in transcriptsToSequence[gene].iterkeys():
			sequence = ""
			for start, end in zip(dictTranscripts[transcript]["start"], dictTranscripts[transcript]["end"]):
				#  concatenate all exons of a transcript #  not sure if end or end+1
				h= str(fasta.seq[start:end]).upper()
				sequence += h
			outputTemp = open(transcript + ".fa", 'w')
			print transcript + ".fa"
			outputTemp.seek(0)
			outputTemp.write(">gene_" + dictTranscripts[transcript]["gene"] + "|transcript_" + transcript + "\n")
			outputTemp.write(sequence)
			outputTemp.close()
			params["file"] = transcript
			if len(sequence)>100: #   required by PBSIM
				if len(sequence) < 4000:
					params["type"] = "ccs"
					print dictTranscripts[transcript]["expr"]
					params["depth"] =  dictTranscripts[transcript]["expr"]
					#~ print "depth", transcriptomeLength/(1500.0 * dictTranscripts[transcript]["expr"])
					# depth = transcriptome length / (readSize * expr)
					bashCommand = ["/home/marchet/bin/pbsim-1.0.3-Linux-amd64-bin/Linux-amd64/bin/pbsim", "%(file)s.fasta" %params, "--prefix", "%(file)s" %params, "--depth", "%(depth)s" %params,  "--length-mean", "1500",  "--length-max", "4000", "--model_qc", "/home/marchet/bin/pbsim-1.0.3-Linux-amd64-bin/data/model_qc_ccs", "--data-type", "CCS"]
					#~ bashCommand = ["/home/marchet/bin/pbsim-1.0.3-Linux-amd64-bin/Linux-amd64/bin/pbsim", "%(file)s.fa" %params, "--prefix", "%(file)s" %params, "--depth", "10",  "--length-mean", "1500",  "--length-max", "4000", "--model_qc", "/home/marchet/bin/pbsim-1.0.3-Linux-amd64-bin/data/model_qc_ccs", "--data-type", "CCS"]
					#~ bashCommand = ["/home/marchet/bin/pbsim-1.0.3-Linux-amd64-bin/Linux-amd64/bin/pbsim", "/home/marchet/RNA_project/BC040516.fa",  "--model_qc", "/home/marchet/bin/pbsim-1.0.3-Linux-amd64-bin/data/model_qc_clr"]
					proc = subprocess.Popen(bashCommand, stdout=FNULL, stderr=subprocess.PIPE)
					outputErr = proc.stderr.read()
				else:
					params["type"] = "clr"
					#~ params["depth"] = dictTranscripts[transcript]["expr"]
					params["depth"] =  transcriptomeLength/(5000.0 * dictTranscripts[transcript]["expr"])
					bashCommand = ["/home/marchet/bin/pbsim-1.0.3-Linux-amd64-bin/Linux-amd64/bin/pbsim", "%(file)s.fasta" %params, "--prefix", "%(file)s" %params, "--depth", "%(depth)s" %params,  "--length-mean", "5000", "--model_qc", "/home/marchet/bin/pbsim-1.0.3-Linux-amd64-bin/data/model_qc_clr"]
					proc = subprocess.Popen(bashCommand, stdout=FNULL, stderr=subprocess.PIPE)
					outputErr = proc.stderr.read()
				
				if not "ERROR" in outputErr:
					#~ print "noerror"
					fastq = open("%(file)s_0001.fastq" %params, 'r')
					linesFq = fastq.readlines()
					for i in range(len(linesFq)):
						if i%4 == 0:
							line = ">gene_" + dictTranscripts[transcript]["gene"] +"|transcript_" + transcript + linesFq[i][1:] 
							outFile.write(line)
						if i%4==1:
							outFile.write(linesFq[i])
				else:
					print outputErr
				#~ remove = "rm %(file)s.fa" % params
				#~ os.system(remove)
			else:
				print "too short"
			
			removeRef =  "rm *_0001.ref" % params
			removeMaf = "rm *_0001.maf" % params
			os.system(removeRef)
			os.system(removeMaf)
			
			
	else:
		print "./readsIsoSeq.py file.gtf file.fasta depth"
#~ def main():
	#~ gtfFile = sys.argv[1]
	#~ fastaFile = sys.argv[2]
	#~ outFile = open("isoseq.fa", 'w')
	#~ transcriptsToSequence = dict()

	#~ previousGene = None
	#~ with open(gtfFile) as infile:
		#~ for line in infile:
			#~ string = line.rstrip().split("\t")
			#~ boolExon = True if string[2] == 'exon' else False
			#~ if boolExon:
				#~ chromosome = string[0]
				#~ startPosition = int(string[3])
				#~ endPosition = int(string[4])
				#~ geneInfo = string[8].split(";")
				#~ posGeneName = geneInfo[0].find('gene_id')
				#~ posTranscriptName = geneInfo[1].find('transcript_id')
				#~ if posGeneName != -1:
					#~ geneName = geneInfo[0][len('gene_id') + 1:][1:-1]
					#~ transcriptName = geneInfo[1][len('transcript_id') + 1:][2:-1]
				#~ if previousGene is not None:
					#~ if geneName in transcriptsToSequence.keys():
						#~ if transcriptName in transcriptsToSequence[geneName].keys():
						#~ #  another exon from the transcript: add it to what will be sequenced
							#~ transcriptsToSequence[geneName][transcriptName]["start"].append(startPosition)
							#~ transcriptsToSequence[geneName][transcriptName]["end"].append(endPosition)
						#~ else:
						#~ #  probability to sequence an alternative transcript from the gene
							#~ probaToSeqAlt = random.uniform(0,1)
							#~ if probaToSeqAlt <= 1:
								#~ transcriptsToSequence[geneName][transcriptName] = {"start": [startPosition], "end": [endPosition]}
					#~ else:
						#~ probaToSeqGene = random.uniform(0,1)
						#~ if probaToSeqGene < 0.4:
							#~ previousGene = geneName
							#~ transcriptsToSequence.setdefault(geneName, {transcriptName: {"start": [startPosition], "end": [endPosition]}})
				#~ else:
					#~ probaToSeqGene = random.uniform(0,1)
					#~ if probaToSeqGene < 0.4:
						#~ transcriptsToSequence.setdefault(geneName, {transcriptName: {"start": [startPosition], "end": [endPosition]}})
						#~ previousGene = geneName
						
	#~ fasta = SeqIO.read(fastaFile, "fasta")
	#~ params = {"file": ""}
	#~ FNULL = open(os.devnull, 'w')
	#~ for gene in transcriptsToSequence.iterkeys():
		#~ if len(transcriptsToSequence[gene]) > 1:  # alternative isoforms
			#~ for transcript in transcriptsToSequence[gene].iterkeys():
				#~ sequence = ""
				#~ for start, end in zip(transcriptsToSequence[gene][transcript]["start"], transcriptsToSequence[gene][transcript]["end"]):
					#~ sequence += fasta.seq[start:end].upper() #  not sure if end or end+1 #  ? un ficher par seq pour pbsim ? // lancer pbsim dans le script python sur chaque sequence puis effacer le fasta de ref, noter tous les transcripts/genes/positions exons dans un fichier
				#~ outputTemp = open(transcript + ".fa", 'w')
				#~ outputTemp.seek(0)
				#~ outputTemp.write(">gene_" + gene + "|" + str(len(transcriptsToSequence[gene])) + "|transcript_" + transcript + "\n")
				#~ outputTemp.write(str(sequence))
				#~ params["file"] = transcript
				#~ if len(sequence) < 2000:
					#~ params["type"] = "ccs"
					#~ bashCommand = ["/home/marchet/bin/pbsim-1.0.3-Linux-amd64-bin/Linux-amd64/bin/pbsim", "%(file)s.fa" %params, "--prefix", "%(file)s" %params, "--depth", "10",  "--length-mean", "1500",  "--length-max", "2000", "--model_qc", "/home/marchet/bin/pbsim-1.0.3-Linux-amd64-bin/data/model_qc_ccs", "--data-type", "CCS"]
				#~ else:
					#~ params["type"] = "clr"
					#~ bashCommand = ["/home/marchet/bin/pbsim-1.0.3-Linux-amd64-bin/Linux-amd64/bin/pbsim", "%(file)s.fa" %params, "--prefix", "%(file)s" %params, "--depth", "10",  "--length-mean", "5000", "--model_qc", "/home/marchet/bin/pbsim-1.0.3-Linux-amd64-bin/data/model_qc_clr"]
				#~ proc = subprocess.Popen(bashCommand, stdout=FNULL, stderr=subprocess.PIPE)
				#~ outputErr = proc.stderr.read()
				#~ if not "ERROR" in outputErr:
					#~ fastq = open("%(file)s_0001.fastq" %params, 'r')
					#~ linesFq = fastq.readlines()
					#~ for i in range(len(linesFq)):
						#~ if i%4 == 0:
							#~ line = ">gene_" + gene + "|" + str(len(transcriptsToSequence[gene])) + "|transcript_" + transcript + "|"  + "%(type)s_" %params + linesFq[i][1:] 
							#~ outFile.write(line)
						#~ if i%4==1:
							#~ outFile.write(linesFq[i])
				#~ remove = "rm %(file)s.fa" % params
				#~ os.system(remove)
	#~ removeRef =  "rm *_0001.ref" % params
	#~ removeMaf = "rm *_0001.maf" % params
	#~ os.system(removeRef)
	#~ os.system(removeMaf)
			
if __name__ == '__main__':
    main()
