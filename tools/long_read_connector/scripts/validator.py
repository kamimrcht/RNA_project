#!/usr/bin/env python
import sys

def intersect(a, b):
	#~ print len(set(a))
	#~ print len(set(b))
	#~ print len(set(a) & set(b))
	#~ print len(set(a) | set(b))
	#~ aa = set(a)
	#~ bb = set(b)
	#~ print bb.issubset(aa), len(list(set(a) & set(b)))
	return list(set(a) & set(b))

simulatedFile = open(sys.argv[1], "r")
lines = simulatedFile.readlines()
l = 0
readNumber = 0
refToReads = dict()
readToRef = []
while (l < len(lines)):
	ref = int(lines[l].split(" ")[0].split(":")[1])
	if ref in refToReads.keys():
		refToReads[ref].append(readNumber)
	else:
		refToReads[ref] = [readNumber]
	readToRef.append([readNumber, ref])
	readNumber += 1
	l += 2  # skip sequence

readToReads = dict()
for couple in readToRef:
	readToReads[couple[0]] = refToReads[couple[1]]
	
	#~ toPrint = str(couple[0]) + ":" + " ".join(str(x) for x in refToReads[couple[1]])
	#~ print toPrint

resultFile = open(sys.argv[2], "r")
res = resultFile.readlines()

recall = 0
specificity = 0
for i in res:
	if "#" in i:
		continue
	refRead =  int(i.split(":")[0]) - 1
	reads = [int(x)-1 for x in i.split(":")[1].split(" ")[:-1]]
	inter = intersect(reads, readToReads[refRead])
	recall += float( len(inter))/len(readToReads[refRead])
	specificity += float(len(inter) - len(reads)) / len(inter)
	#~ print recall, specificity
	#~ print len(readToReads[refRead]), len(reads)
print "****************************************************************************"
recall = float(recall)/len(res)
specificity = 1 + float(specificity)/len(res)
print recall, specificity
	#~ print refRead, len(readToReads[refRead]), len(reads), len(inter), len(reads)-len(inter), len(readToReads[refRead]) - len(reads)
	
	
	
