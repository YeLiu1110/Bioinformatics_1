import time
import os
import sys
import pandas as pd
import math
import random
import numpy as np

def StringComposition(k, Text):
    kmerlist = [Text[i:i+k] for i in range(len(Text)-k+1)]
    return sorted(kmerlist)

def ReconstructGenome(PatternList):
    '''
    String Spelled by a Genome Path Problem
    '''
    Pattern = list(PatternList[0])
    for i in PatternList[1:]:
        Pattern.append(i[-1])
    return ''.join(Pattern)

def Overlap(Patterns):
    OverlapDic = {i:[] for i in Patterns}
    for i in Patterns:
        for j in Patterns:
            if i[1:] == j[:-1]:
                OverlapDic[i].append(j)
    OverlapDic1 = {}
    for i in OverlapDic.keys():
        if OverlapDic[i] != []:
            OverlapDic1.update({i:sorted(OverlapDic[i])})
    return OverlapDic1

def sew(Patterns):
    Pattern2 = []
    Pattern1 = Patterns
    while Pattern1:
        for i in Pattern1:
            Pattern2.append(i)
            Pattern1.remove(i)

            for le in Pattern1:
                if le[:-1] == ini[1:]:
                    Pattern2.append(le)
                    ini = le
            import pdb;pdb.set_trace()
    return Pattern2

def PathGraph(Text,k):
    PathDic = {}
    for i in range(len(Text)-k+1):
        if Text[i:i+k-1] in PathDic.keys():
            PathDic[Text[i:i+k-1]].append(Text[i+1:i+k])
        else:
            PathDic.update({Text[i:i+k-1]:[Text[i+1:i+k]]})
    PathDic1 = {}
    sortedkeylist = sorted(PathDic.keys())
    for i in sortedkeylist:
        PathDic1.update({i:sorted(PathDic[i])})
    return PathDic1

def DeBruijn(Patterns):
    DeDic = {}
    for i in Patterns:
        if i[:-1] in DeDic.keys():
            DeDic[i[:-1]].append(i[1:])
        else:
            DeDic.update({i[:-1]:[i[1:]]})
    DeDic1 = {}
    SortedKeyList = sorted(DeDic.keys())
    for i in SortedKeyList:
        DeDic1.update({i:sorted(DeDic[i])})
    return DeDic1

starttime = time.time()
'''
k = 8
t = 5
Dna = [
'CGCCCCTCTCGGGGGTGTTCAGTAACCGGCCA',
'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
'AATCCACCAGCTCCACGTGCAATGTTGGCCTA'
]
Patterns = [
    '0000','1000','0100','1101',
    '0001','0011','1010','1011',
    '0010','0101','1100','0111',
    '0100','1001','1110','1111'
]

fil = open(sys.argv[1],'r')
linlist = []
for lins in fil:
    linlist.append(lins.strip('\r\n'))
Patterns = linlist

Patterns = [
    '000','001','010','100','011','101','110','111'
]
'''
Patterns = [
'ATGCG',
'GCATG',
'CATGC',
'AGGCA',
'GGCAT'
]
job = Overlap(Patterns)

print ('job done by {} s.'.format(time.time()-starttime))
'''
cwd = os.getcwd()

txtFile1 = cwd + '/results.txt'

with open(txtFile1,'w') as f:
    for keys,val in job.items():
        f.write(('{} -> {}\n').format(keys,','.join(val)))
f.close()
'''
import pdb;pdb.set_trace()
for i in job:   print (i)