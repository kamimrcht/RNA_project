#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  kmers_occurrences.py
#  
#  Copyright 2016 camille marchet <camille.marchet@irisa.fr>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  


def main(args):
    k = 11
    while k < 12:
        indexRead = -1
        kmersMap = dict()
        readsToKmers = dict()
        with open(sys.argv[1]) as infile:
            for line in infile:
                sequence = line.rstrip()
                if not ">" in sequence:
                    indexRead += 1
                    kmersSet = set()  # to keep only kmers unique in the sequence
                    for i in range(len(sequence) - k + 1):
                        kmer = sequence[i:i+k]
                        if kmer not in kmersSet:
                            kmersSet.add(kmer)
                            readsToKmers[indexRead] = kmersSet
                            if kmer in kmersMap.keys():
                                kmersMap[kmer] += 1
                            else:
                                kmersMap[kmer] = 1
        print "***** ", k, " *****"
        # stats
        #~ occurrence = []
        #~ for nb in kmersMap.itervalues():
            #~ occurrence.append(nb)
        
        #~ occurr = sorted(occurrence)
        #~ c = 0
        #~ values = []
        #~ for val in occurr:
            #~ if val > c:
                #~ values.append([val, occurr.count(val)])
                #~ c = val
        #~ sumV = 0
        #~ for v in values:
            #~ sumV += v[1]
        #~ acc = sumV - values[0][1]
        #~ print min(occurrence), max(occurrence), float(sum(occurrence))/len(occurrence), float(acc)/len(occurrence)

        # 1 v 1 read comparison
        listCommon = []
        checked = 0
        listIndex = range(0, indexRead)
        for i in range(50):
            random.shuffle(listIndex)
            index = listIndex[0]
            listIndex.remove(index)
        
            for read in readsToKmers.iterkeys():
            #~ for readComp in readsToKmers.iterkeys():
                i = len(readsToKmers[index].intersection(readsToKmers[read]))  # nb of kmers in common in two reads
                listCommon.append(i)
        print min(listCommon), max(listCommon), float(sum(listCommon))/len(listCommon)
        
        k += 2
    return 0

if __name__ == '__main__':
    import sys
    from numpy import random
    sys.exit(main(sys.argv))
