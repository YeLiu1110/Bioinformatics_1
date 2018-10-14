import time
import os
import sys
import pandas as pd
import math
import random
import numpy as np

def EulerianCycle(graphdic):
    def Cycle(graphdic,initialCycle, newDic, First = False):
        randomStart = initialCycle[0]
        while initialCycle.count(initialCycle[0]) != 2:
            stri = graphdic[randomStart][0]
            initialCycle.append(stri)
            if First:
                if len(graphdic[randomStart]) == 1:
                    newDic.update({randomStart: graphdic[randomStart]})
                    graphdic.pop(randomStart)

                else:
                    newDic.update({randomStart: [graphdic[randomStart][0]]})
                    del (graphdic[randomStart][0])
            else:

                if (len(graphdic[randomStart]) == 1) & (randomStart in newDic):
                    newDic[randomStart].extend(graphdic[randomStart])
                    graphdic.pop(randomStart)
                elif (len(graphdic[randomStart]) == 1) & ((randomStart in newDic) == False):
                    newDic.update({randomStart:graphdic[randomStart]})
                    graphdic.pop(randomStart)
                else:

                    newDic[randomStart].append(graphdic[randomStart][0]) if (randomStart in newDic) else newDic.update({randomStart:[graphdic[randomStart][0]]})
                    del (graphdic[randomStart][0])
            randomStart = stri

        return graphdic, initialCycle, newDic
    def AddingList(TotList,AddList, position):
        for i in range(len(AddList)-1):
            TotList.insert(position + i, AddList[i])
        return TotList
    keysList = list(graphdic.keys())
    randomStart = keysList[random.randint(0,len(keysList)-1)]

    initialCycle = [randomStart]
    popDic = {}
    graphdic, initialCycle, popDic = Cycle(graphdic, initialCycle, popDic, True)

    while graphdic:
        for i in set(initialCycle):
            if i in graphdic:
                position = initialCycle.index(i)
                newCycle = [i]

                graphdic, newCycle, popDic = Cycle(graphdic, newCycle, popDic)
                initialCycle = AddingList(initialCycle,newCycle,position)
            else:
                pass
    return initialCycle

def EulerianPath(graphdic):
    def FindingStart(graphdic):
        keysList = list(graphdic.keys())
        valuelist = []
        for i in graphdic.values():
            valuelist.extend(i)
        for i in keysList:
            if len(graphdic[i]) >  valuelist.count(i):
                Start = i
        for i in valuelist:
            if ((i in keysList)==False) or (len(graphdic[i]) <  valuelist.count(i)):
                End = i
        import pdb;pdb.set_trace()
        return Start, End

    def Cycle(graphdic,initialCycle, newDic, End, First = False ):
        Start = initialCycle[0]
        while (initialCycle[-1]!=End) and (initialCycle.count(initialCycle[0]) != 2):
            stri = graphdic[Start][0]
            initialCycle.append(stri)
            if First:
                if len(graphdic[Start]) == 1:
                    newDic.update({Start: graphdic[Start]})
                    graphdic.pop(Start)
                else:
                    newDic.update({Start: [graphdic[Start][0]]})
                    del (graphdic[Start][0])
            else:
                if (len(graphdic[Start]) == 1) & (Start in newDic):
                    newDic[Start].extend(graphdic[Start])
                    graphdic.pop(Start)
                elif (len(graphdic[Start]) == 1) & ((Start in newDic) == False):
                    newDic.update({Start:graphdic[Start]})
                    graphdic.pop(Start)
                else:
                    newDic[Start].append(graphdic[Start][0]) if (Start in newDic) else newDic.update({Start:[graphdic[Start][0]]})
                    del (graphdic[Start][0])
            Start = stri
        return graphdic, initialCycle, newDic
    def AddingList(TotList,AddList, position):
        for i in range(len(AddList)-1):
            TotList.insert(position + i, AddList[i])
        return TotList
    Start, End = FindingStart(graphdic)
    initialCycle = [Start]
    popDic = {}
    graphdic, initialCycle, popDic = Cycle(graphdic, initialCycle, popDic, End, True)
    while graphdic:
        for i in set(initialCycle):
            if i in graphdic:
                position = initialCycle.index(i)
                newCycle = [i]
                graphdic, newCycle, popDic = Cycle(graphdic, newCycle, popDic, End)
                initialCycle = AddingList(initialCycle,newCycle,position)
            else:
                pass
    return initialCycle

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

def PairedDeBruijn(GappedPatterns):
    DeDic = {}
    for i in GappedPatterns:
        Prefix = i.split('|')[0][:-1] + '|' + i.split('|')[1][:-1]
        Suffix = i.split('|')[0][1:] + '|' + i.split('|')[1][1:]
        if Prefix in DeDic:
            DeDic[Prefix].append(Suffix)
        else:
            DeDic.update({Prefix:[Suffix]})
    DeDic1 = {}
    SortedKeyList = sorted(DeDic.keys())
    for i in SortedKeyList:
        DeDic1.update({i:sorted(DeDic[i])})
    return DeDic1

def ReconstructGenome(PatternList):
    '''
    String Spelled by a Genome Path Problem
    '''
    Pattern = list(PatternList[0])
    for i in PatternList[1:]:
        Pattern.append(i[-1])
    return ''.join(Pattern)

def StringSpelledByGappedPatterns(GappedPatterns, k, d):
    FirstSequence = [i.split('|')[0] for i in GappedPatterns]
    SeconSequence = [i.split('|')[1] for i in GappedPatterns]
    PrefixSequence = ReconstructGenome(FirstSequence)
    SuffixSequence = ReconstructGenome(SeconSequence)
    for i in range(k+d,len(PrefixSequence)):
        if PrefixSequence[i] != SuffixSequence[i-k-d]:
            return "there is no string spelled by the gapped patterns"
    return PrefixSequence + SuffixSequence[-(k+d):]

def MaximalNonBranchingPaths(graph):
    def TestSameLoop(list1,list2):
        if len(list1) != len(list2):
            return False
        else:
            for i in range(len(list1)):
                if list1[i] not in list2:
                    return False
                else:
                    list3 = list2 + list2[:list2.index(list1[i])]
                    ind = []
                    ind.append(i - list3.index(list1[i]))
            if ind.count(ind[0]) == len(ind):
                return True
            else:
                return False
    Paths = []
    keysList = list(graph.keys())
    valueList = []
    usedkeysList = []
    for i in list(graph.values()):
        valueList.extend(i)
    for i in graph:
        if (len(graph[i]) != 1) or (valueList.count(i) != 1):
            if len(graph[i]) > 0:
                usedkeysList.append(i)
                for edge in graph[i]:
                    NonBranchingPath = [i,edge]
                    while (valueList.count(edge) == 1) and (edge in keysList) and (len(graph[edge]) == 1):
                        NonBranchingPath.extend(graph[edge])
                        usedkeysList.append(edge)
                        edge = graph[edge][0]
                    Paths.append(NonBranchingPath)
    Cycle = []
    notusedkeysList = set(keysList)-set(usedkeysList)
    for i in notusedkeysList:
        if (len(graph[i]) == 1) and (valueList.count(i) == 1):
            edge = graph[i][0]
            cycle = [i,edge]
            while (valueList.count(edge) == 1) and (len(graph[edge]) == 1) and (cycle[-1] != cycle[0]):
                cycle.extend(graph[edge])
                edge = graph[edge][0]
                if edge not in keysList:
                    break
            if len(Cycle) == 0:
                Cycle.append(cycle)
            else:
                TrueList = []
                for cy in Cycle:
                    TrueList.append(TestSameLoop(cycle,cy))
                if TrueList.count(False) == len(TrueList):
                    Cycle.append(cycle)
    for i in Cycle:
        Paths.append(i)
    return Paths

starttime = time.time()
'''
fil = open(sys.argv[1],'r')
linlist = []
for lins in fil:
    linlist.append(lins.strip('\r\n'))
Patterns = linlist[1:]
FirstPatterns,SecondPatterns = [],[]
for i in Patterns:
    FirstPatterns.append(i.split('|')[0])
    SecondPatterns.append(i.split('|')[1])
k,d = linlist[0].split(' ')

GraphDic = {}
for i in linlist:
    ke, value = i.split('->')
    valuelist = value.strip(' ').split(',')
    GraphDic.update({ke.strip(' '):valuelist})
print ('Dic done by {} s.'.format(time.time()-starttime))
#import pdb;pdb.set_trace()

GraphDic = {
'0' : ['2'],
     '1' : ['3'],
     '2' : ['1'],
     '3' : ['0','4'],
     '6' : ['3','7'],
     '7' : ['8'],
     '8' : ['9'],
     '9' : ['6']
}

Patterns = [
    'CTTA',
     'ACCA',
     'TACC',
     'GGCT',
     'GCTT',
     'TTAC'
]

k,d = 4,2
Patterns = [
'GAGA|TTGA',
'TCGT|GATG',
'CGTG|ATGT',
'TGGT|TGAG',
'GTGA|TGTT',
'GTGG|GTGA',
'TGAG|GTTG',
'GGTC|GAGA',
'GTCG|AGAT'
]

Patterns = [
'ATG',
'ATG',
'TGT',
'TGG',
'CAT',
'GGA',
'GAT',
'AGA',


]
'''
Patterns = [
'AAAT',
'AATG',
'ACCC',
'ACGC',
'ATAC',
'ATCA',
'ATGC',
'CAAA',
'CACC',
'CATA',
'CATC',
'CCAG',
'CCCA',
'CGCT',
'CTCA',
'GCAT',
'GCTC',
'TACG',
'TCAC',
'TCAT',
'TGCA'
]

Patterns = [
'ACC|ATA',
'ACT|ATT',
'ATA|TGA',
'ATT|TGA',
'CAC|GAT',
'CCG|TAC',
'CGA|ACT',
'CTG|AGC',
'CTG|TTC',
'GAA|CTT',
'GAT|CTG',
'GAT|CTG',
'TAC|GAT',
'TCT|AAG',
'TGA|GCT',
'TGA|TCT',
'TTC|GAA']


GraphDic = PairedDeBruijn(Patterns)

print ('Dic done by {} s.'.format(time.time()-starttime))
import pdb;pdb.set_trace()
job1 = EulerianPath(GraphDic)
finaljob = StringSpelledByGappedPatterns(job1,3,1)
print ('EulerianPath done by {} s.'.format(time.time()-starttime))
#import pdb;pdb.set_trace()

print ('job done by {} s.'.format(time.time()-starttime))

cwd = os.getcwd()

txtFile1 = cwd + '/results.txt'

with open(txtFile1,'w') as f:


    f.write(finaljob)



f.close()

import pdb;pdb.set_trace()
for i in job:   print (i)