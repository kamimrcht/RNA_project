#!/usr/bin/env python
import sys

def readSRCFile(fileName, clusters):
    srcfile = open(fileName, "r")
    for line in srcfile.readlines():
        if line[0] == '#': #header
            continue
        line = line.rstrip()
        reads = line.split(' ')
        query_read_id = int(reads[0]) - 1
        for read in reads:
            if query_read_id in clusters.keys():
                clusters[query_read_id].add(int(read)-1)
            else:
                clusters[query_read_id] = set([int(read)-1])
        

            
def readSimulFile(fileName, clustersSimu,  types):
    simulFile = open(fileName, "r")
    readNumber = 0
    for line in simulFile.readlines():
        if line[0] == '>': #header >referenceNumber:1 alternativeNumber1
            ref = line.split(' ')[0].split(':')[1]
            if ref in clustersSimu.iterkeys():
                clustersSimu[ref].add(readNumber)
            else:
                clustersSimu[ref] = set([readNumber])
                types[ref] = line.split(' ')[2][:-1]
            readNumber += 1

def validateResults(clustersLRC, clustersSimu,  outfile, types):
    rec = {"highExpression":[], "lowExpression":[], "regularExpression":[]}
    pre = {"highExpression":[], "lowExpression":[], "regularExpression":[]}
    GLOBALR = []
    GLOBALP = []
    for keyRefSimu in clustersSimu.iterkeys():
        readsFound = set()
        recalls = []
        precisions = []
        typeExpr = types[keyRefSimu]
        #~ print typeExpr
        for keyLRC in clustersLRC.iterkeys() :
            if keyLRC in clustersSimu[keyRefSimu]:
                intersec = len(clustersLRC[keyLRC].intersection(clustersSimu[keyRefSimu]))
                recall = round(intersec * 1.0/len(clustersSimu[keyRefSimu]) * 100, 2)
                precision = round(intersec * 1.0 / len(clustersLRC[keyLRC]) * 100, 2)
                recalls.append(recall)
                if intersec != 0:
                    precisions.append(precision)
                readsFound.update(clustersLRC[keyLRC].intersection(clustersSimu[keyRefSimu]))

        if len(readsFound) != 0:
            RECALLS = round(len(readsFound) * 1.0 /len(clustersSimu[keyRefSimu]) * 100, 2)
        else:
            RECALLS = 0
        if len(precisions) != 0:
            PREC = round(sum(precisions) * 1.0/len(precisions), 2)
            
        else:
            PREC= "NA"
        
        #~ print "gene ", typeExpr, keyRefSimu,RECALLS, PREC
        GLOBALR.append(RECALLS)
        rec[typeExpr].append(RECALLS)
        if PREC != "NA":
            GLOBALP.append(PREC)
            pre[typeExpr].append(PREC)
    for key in rec.iterkeys():
        if len(rec[key]) != 0:
            if len(pre[key]) != 0:
                print key, round(sum(rec[key]) * 1.0/len(rec[key]), 2),  round(sum(pre[key]) * 1.0/len(pre[key]), 2)
            else:
                print key, round(sum(rec[key]) * 1.0/len(rec[key]), 2), 0
        else:
            print key, 0, "NA"
    print "global", round(sum(GLOBALR) * 1.0/len(GLOBALR), 2), round(sum(GLOBALP) * 1.0/len(GLOBALP), 2)
        
        




if len(sys.argv) < 3 or len(sys.argv) > 3:
     print("COMMAND")
     print(sys.argv[0]," <SRC_linker output file> <simulation file>")
     sys.exit(1)

SRC_output = sys.argv[1]
simul = sys.argv[2]
clustersLRC = dict()
clustersSimu = dict()
types = dict()
readSRCFile(SRC_output, clustersLRC)
readSimulFile(simul, clustersSimu, types)
globalValues = open("/home/marchet/RNA_project/expe/global_results.txt", 'w')
validateResults(clustersLRC, clustersSimu, globalValues, types)


