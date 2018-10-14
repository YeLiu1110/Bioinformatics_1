import time
import os
import sys
import pandas as pd
import math

def NumberToPattern(numbe, divvalue):
    stalst, indlst = ['A', 'C', 'G', 'T'], [0, 1, 2, 3]
    numbe = int(numbe)
    mid = []
    while True:
        if numbe == 0: break
        numbe, rem = divmod(numbe, 4)
        mid.append(rem)
    if len(mid) < divvalue:
        x=0
        appendlst = []
        for x in range(divvalue - len(mid)):
            appendlst.append(0)
            x=x+1
    #import pdb;pdb.set_trace()
        midtot = mid + appendlst
    else:
        midtot = mid
    numstr = ''.join([stalst[indlst.index(x)] for x in midtot[::-1]])
    return numstr

def HammingDistance(Seq1, Seq2):
    if len(Seq1) != len(Seq2):
        #import pdb;pdb.set_trace()
        return
    else:
        mismatchcount = 0
        for i in range(0,len(Seq1)):
            if Seq1[i] != Seq2[i]:
                mismatchcount =mismatchcount + 1
        return mismatchcount

def Neighbors(Pattern,d):

    def Suffix(Pattern):
        return Pattern[1:]
    if d == 0:
        return [Pattern]
    if len(Pattern) == 1:
        return ['A','T','C','G']
    Neighborhood = set()

    #import pdb;pdb.set_trace()
    SuffixNeighbors = Neighbors(Suffix(Pattern),d)
    for str in SuffixNeighbors:
        #import pdb;pdb.set_trace()
        if HammingDistance(Suffix(Pattern),str) < int(d):
            for x in ['A','T','C','G']:
                Neighborhood.add('%s%s'%(x,str))
            #import pdb;pdb.set_trace()

        else:
            Neighborhood.add('%s%s' % (Pattern[len(Pattern)-len(str)-1], str))
    return Neighborhood

def MotifEnumeration(Dna, k, d):
    '''
    Dna: a list of dna strings
    k: integer
    d: integer
    '''
    def KmerInList(kmer,list):
        return int(kmer in list)
    Patterns = set()
    dic = {}
    for num in range(len(Dna)):
        name = 'String' + str(num)
        KmerSet = set()
        Kmer1Set = set()
        for i in range(len(Dna[num])-k+1):
            KmerSet.add(Dna[num][i:i+k])
            for kmer in KmerSet:
                Kmer1Set = Kmer1Set | set(Neighbors(kmer,d))
        dic.update({name:Kmer1Set})

    for kmer in dic['String0']:
        tot = 1
        for i in range(1,len(Dna)):
            name = 'String' + str(i)
            tot = tot + KmerInList(kmer, dic[name])
        if tot == len(Dna):
            Patterns.add(kmer)

    return Patterns
def ProfileForming(Motifs):
    '''
    Motifs: a list of sequence string
    '''
    def ProbabilityDistribution(numlist):
        return sum([ i for i in numlist])

    MatriceDF = pd.DataFrame([list(i) for i in Motifs])
    BaseDF = pd.DataFrame([[0]*len(Motifs[0])]*4,index = ['A','C','G','T'])
    ScoreDF = MatriceDF.append(BaseDF)
    consensus = []
    for i in ScoreDF.columns.tolist():
        CountDic = {j: ScoreDF[i][0:len(Motifs)].tolist().count(j) for j in ['A', 'T', 'C', 'G']}

        consensus.append(max(CountDic, key=CountDic.get))
    consensusStrig = ''.join(consensus)

    for i in ScoreDF.columns.tolist():
        for j in range(len(Motifs)):
            if ScoreDF[i][j] == 'A':
                ScoreDF[i]['A'] = ScoreDF[i]['A'] + 1
            elif ScoreDF[i][j] == 'C':
                ScoreDF[i]['C'] = ScoreDF[i]['C'] + 1
            elif ScoreDF[i][j] == 'G':
                ScoreDF[i]['G'] = ScoreDF[i]['G'] + 1
            elif ScoreDF[i][j] == 'T':
                ScoreDF[i]['T'] = ScoreDF[i]['T'] + 1

    ScoreDF.loc['A':'T'] = ScoreDF.loc['A':'T'].apply((lambda x: x/len(Motifs)))
    #import pdb;pdb.set_trace()
    #HammingDistanceList = [HammingDistance(''.join(ScoreDF.loc[i].tolist()), consensusStrig) for i in range(len(Motifs))]


    ScoreWantDF = ScoreDF.loc['A':'T']
    ScoreDic = {i: ScoreWantDF.loc[i].tolist() for i in ScoreWantDF.index.tolist()}
    #import pdb;pdb.set_trace()
    return ScoreDic
def ProfileForming1(Motifs):
    '''
    Motifs: a list of sequence string
    '''
    def ProbabilityDistribution(numlist):
        return sum([ i for i in numlist])

    MatriceDF = pd.DataFrame([list(i) for i in Motifs])
    BaseDF = pd.DataFrame([[0]*len(Motifs[0])]*4,index = ['A','C','G','T'])
    ScoreDF = MatriceDF.append(BaseDF)
    consensus = []
    for i in ScoreDF.columns.tolist():
        CountDic = {j: ScoreDF[i][0:len(Motifs)].tolist().count(j) for j in ['A', 'T', 'C', 'G']}

        consensus.append(max(CountDic, key=CountDic.get))
    consensusStrig = ''.join(consensus)

    for i in ScoreDF.columns.tolist():
        for j in range(len(Motifs)):
            if ScoreDF[i][j] == 'A':
                ScoreDF[i]['A'] = ScoreDF[i]['A'] + 1
            elif ScoreDF[i][j] == 'C':
                ScoreDF[i]['C'] = ScoreDF[i]['C'] + 1
            elif ScoreDF[i][j] == 'G':
                ScoreDF[i]['G'] = ScoreDF[i]['G'] + 1
            elif ScoreDF[i][j] == 'T':
                ScoreDF[i]['T'] = ScoreDF[i]['T'] + 1

    ScoreDF.loc['A':'T'] = ScoreDF.loc['A':'T'].apply((lambda x: (x+1)/(len(Motifs)+4)))
    #import pdb;pdb.set_trace()
    #HammingDistanceList = [HammingDistance(''.join(ScoreDF.loc[i].tolist()), consensusStrig) for i in range(len(Motifs))]


    ScoreWantDF = ScoreDF.loc['A':'T']
    ScoreDic = {i: ScoreWantDF.loc[i].tolist() for i in ScoreWantDF.index.tolist()}
    #import pdb;pdb.set_trace()
    return ScoreDic

def minHammingDistance(Pattern,Text):
    HDlist = []
    #import pdb;pdb.set_trace()
    if len(Text) == len(Pattern):
        HDlist.append(HammingDistance(Pattern,Text))
    else:
        for i in range(len(Text)-len(Pattern)+1):
            HDlist.append(HammingDistance(Pattern,Text[i:i+len(Pattern)]))
    #import pdb;pdb.set_trace()
    ind = HDlist.index(min(HDlist))
    return min(HDlist)
def sumofminHammingDistance(Pattern,TextList):
    Distance = []
    for Text in TextList:
        Distance.append(minHammingDistance(Pattern,Text))
    return sum(Distance)

def MedianString(Dna, k):
    distance = float("inf")
    kmerslist = [NumberToPattern(i,k) for i in range(4**k)]
    MedianList = []
    for kmer in kmerslist:
        MedianList.append(sumofminHammingDistance(kmer,Dna))
        if distance > sumofminHammingDistance(kmer,Dna):
            distance = sumofminHammingDistance(kmer,Dna)
            Median = kmer
    #import pdb;pdb.set_trace()

    MedianNeedList = []
    for ind,num in enumerate(MedianList):
        if num == min(MedianList):

            MedianNeedList.append(kmerslist[ind])
    import pdb;pdb.set_trace()
    return Median

def Probability(Pattern,profile):
    '''
    profile: a dictionary containing 'A','C','G','T' and its probability
    '''
    score = 1
    for i in range(len(Pattern)):
        score = score * float(profile[Pattern[i]][i])
    #import pdb;pdb.set_trace()
    return score

def ProfileMostProbable(Text, k, matrix):
    kmerlist = []
    for i in range(len(Text)-k):
        kmerlist.append(Text[i:i+k])
    #import pdb;pdb.set_trace()
    maxscore = 0.0
    kmermost = Text[0:k]
    for kmer in kmerlist:
        prob = Probability(kmer,matrix)
        if maxscore < prob:
            maxscore = prob
            kmermost = kmer

    return kmermost

def Score(Motifs):
    MatriceDF = pd.DataFrame([list(i) for i in Motifs])
    Scorelist = []
    for i in MatriceDF.columns.tolist():
        maxCount = max([MatriceDF[i].tolist().count(j) for j in ['A', 'C', 'G', 'T'] ])
        Scorelist.append(len(Motifs)-maxCount)
    return sum(Scorelist)

def GreedyMotifSearch(Dna, k, t):

    BestMotifs = [Dna[i][0:k] for i in range(len(Dna))]
    #import pdb;pdb.set_trace()
    for i in range(len(Dna[0])-k+1):
        motif = [Dna[0][i:i+k]]
        for j in range(1,t):
            profile = ProfileForming1(motif)
            motif.append(ProfileMostProbable(Dna[j],k, profile))
        if Score(motif) < Score(BestMotifs):
            BestMotifs = motif
    #import pdb;pdb.set_trace()

    return BestMotifs

def DistanceBetweenPatternAndStrings(Pattern,Dna):

    k = len(Pattern)
    distance = 0
    for text in Dna:
        hammingdistance = float("inf")
        for i in range(len(text)-k+1):
            kmer = text[i:i+k]
            #import pdb;pdb.set_trace()
            eachdistance = HammingDistance(Pattern,kmer)
            if hammingdistance > eachdistance:
                hammingdistance = eachdistance
        distance = distance + hammingdistance
    return distance

starttime = time.time()

fil = open(sys.argv[1],'r')
linlist = []
for lins in fil:
    linlist.append(lins.strip('\r\n'))
Pattern = linlist[0]
Dna = linlist[1].split(' ')
'''
for ind,bas in enumerate(['A','C','G','T']):
    #import pdb;pdb.set_trace()
    kmermatrice.update({bas:linlist[ind+2].split(' ')})
#Dna = ['AAAAA','AAAAA','AACAA']

Motifs = [
"TCGGGGGTTTTT",
"CCGGTGACTTAC",
"ACGGGGATTTTC",
"TTGGGGACTTTT",
"AAGGGGACTTCC",
"TTGGGGACTTCC",
"TCGGGGATTCAT",
"TCGGGGATTCCT",
"TAGGGGAACTAC",
"TCGGGTATAACC"
]

DNA = ['TGATGATAACGTGACGGGACTCAGCGGCGATGAAGGATGAGT','CAGCGACAGACAATTTCAATAATATCCGCGGTAAGCGGCGTA','TGCAGAGGTTGGTAACGCCGGCGACTCGGAGAGCTTTTCGCT','TTTGTCATGAACTCAGATACCATAGAGCACCGGCGAGACTCA','ACTGGGACTTCACATTAGGTTGAACCGCGAGCCAGGTGGGTG','TTGCGGACGGGATACTCAATAACTAAGGTAGTTCAGCTGCGA','TGGGAGGACACACATTTTCTTACCTCTTCCCAGCGAGATGGC','GAAAAAACCTATAAAGTCCACTCTTTGCGGCGGCGAGCCATA','CCACGTCCGTTACTCCGTCGCCGTCAGCGATAATGGGATGAG','CCAAAGCTGCGAAATAACCATACTCTGCTCAGGAGCCCGATG']


kmermatrice = {
'A':[0.7, 0.2, 0.1, 0.5, 0.4, 0.3, 0.2, 0.1],
'C':[0.2, 0.2, 0.5, 0.4, 0.2, 0.3, 0.1, 0.6],
'G':[0.1, 0.3, 0.2, 0.1, 0.2, 0.1, 0.4, 0.2],
'T':[0.0, 0.3, 0.2, 0.0, 0.2, 0.3, 0.3, 0.1]}
profile = {

    'A': [0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.9, 0.1, 0.1, 0.1, 0.3, 0.0],
    'C': [0.1, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.1, 0.2, 0.4, 0.6],
    'G': [0.0, 0.0, 1.0, 1.0, 0.9, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0],
    'T': [0.7, 0.2, 0.0, 0.0, 0.1, 0.1, 0.0, 0.5, 0.8, 0.7, 0.3, 0.4]
}

Dna = [
'GCACATCATTATCGATAACGATTCATTGCCAGGCGGCCGC',
'TCATCGAATAACTGACACCTGCTCTGGCTCATCCGACCGC',
'TCGGCGGTATAGCCAGATAGTGCCAATAATTTCCTAAGCG',
'GTGGTGAAGTGTGGGTTATGGGGAAAGGCAGACTGAGTCG',
'GACGGCAACTACGGTTACAACGCAGCAAGAATATTAACCG',
'TCTGTTGTTGCTAACACCGTTAAGCGACGGCAACTAGGCG',
'GCCAACGTAGGCGCGGCTTGGCATCTCGGTGTGTGAAGCG',
'AAAGGCGCATCTTACTCTTTTCGCTTTCAAAAAAAAATTG']
'''

Dna = [
    'CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCC',
    'GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTC',
    'GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG'
]
job = MedianString(Dna, 7)

#job, job1 = minHammingDistance('GATTCTCA', 'GCAAAGACGCTGACCAA')
#print (' '.join(job))
print ('job done by {} s.'.format(time.time()-starttime))
import pdb;pdb.set_trace()