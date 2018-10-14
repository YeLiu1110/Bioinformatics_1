import time
import os
import sys
import pandas as pd
import math
import random
import numpy as np

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
    #import pdb;pdb.set_trace()
    def ProbabilityDistribution(numlist):
        return sum([ i for i in numlist])

    MatriceDF = pd.DataFrame([list(i) for i in Motifs])
    BaseDF = pd.DataFrame([[0]*len(Motifs[0])]*4,index = ['A','C','G','T'])
    ScoreDF = MatriceDF.append(BaseDF)

    for i in ['A','C','G','T']:
        for j in ScoreDF.columns.tolist():
            #import pdb;pdb.set_trace()
            ScoreDF[j][i] = ScoreDF[j][:(int(t)-1)].tolist().count(i)
    '''
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
    '''
    ScoreDF.loc['A':'T'] = ScoreDF.loc['A':'T'].apply((lambda x: (x+1)/(len(Motifs)+4)))
    #import pdb;pdb.set_trace()
    #HammingDistanceList = [HammingDistance(''.join(ScoreDF.loc[i].tolist()), consensusStrig) for i in range(len(Motifs))]


    ScoreWantDF = ScoreDF.loc['A':'T']
    ScoreDic = {i: ScoreWantDF.loc[i].tolist() for i in ScoreWantDF.index.tolist()}

    return ScoreDic

def ProfileForming2(Motifs):

    ArrayOri = np.array([list(i) for i in Motifs])
    ArrayCount = np.transpose(ArrayOri)
    #import pdb;pdb.set_trace()


    ScoreDic = {}
    CountList = []
    for i in range(len(ArrayCount)):
        CountList.append([ArrayCount[i].tolist().count(j) for j in ['A','C','G','T']])
    for i in range(len(CountList)):
        Eachsum = sum(CountList[i])
        for j in range(len(CountList[i])):
            CountList[i][j] = (CountList[i][j] + 1)/(Eachsum +4)
    ArrayScore = np.transpose(np.array(CountList))

    for ind, i in enumerate(['A','C','G','T']):
        ScoreDic.update({i:ArrayScore[ind].tolist()})

    return ScoreDic

def Probability(Pattern,profile):
    '''
    profile: a dictionary containing 'A','C','G','T' and its probability
    '''
    score = 1
    for i in range(len(Pattern)):
        score = score * float(profile[Pattern[i]][i])
    #import pdb;pdb.set_trace()
    return score

def Score(Motifs):
    MatriceDF = pd.DataFrame([list(i) for i in Motifs])
    Scorelist = []
    for i in MatriceDF.columns.tolist():
        maxCount = max([MatriceDF[i].tolist().count(j) for j in ['A', 'C', 'G', 'T'] ])
        Scorelist.append(len(Motifs)-maxCount)
    return sum(Scorelist)

def ProfileMostProbable(Text, k, matrix):
    kmerlist = []
    for i in range(len(Text)-k+1):
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

def Motifs(profile, Dna, k):
    Motifslist = []
    for i in range(len(Dna)):
        #import pdb;pdb.set_trace()
        Motifslist.append(ProfileMostProbable(Dna[i], k, profile))
    return Motifslist
def RandomizedMotifSearch(Dna, k, t):
    motifs = []
    for i in range(t):
        RandomNum = random.randint(0,len(Dna[i])-k-1)
        motifs.append(Dna[i][RandomNum:RandomNum+k])
    BestMotifs = motifs
    var = 1
    while var == 1:

        Profile = ProfileForming2(motifs)
        #import pdb;pdb.set_trace()
        motifs = Motifs(Profile, Dna, k)
        if Score(motifs) < Score(BestMotifs):
            BestMotifs = motifs
            #import pdb;pdb.set_trace()
        else:
            return BestMotifs, Score(BestMotifs)

def GibbsSampler(Dna, k, t, N):
    def ProfileRandomlyGenaratedKmer(profile, text):

        proList = []
        for ind in range(len(text) - k + 1):
            proList.append(Probability(text[ind:ind+k],profile))
        return proList

    def weight_choice(text, weightList,k):

        rnd = random.random() * sum(weightList)
        for i, w in enumerate(weightList):
            rnd -= w
            if rnd < 0:
                return text[i:i+k]


    motifs = []
    for i in range(t):
        RandomNum = random.randint(0, len(Dna[i]) - k - 1)
        motifs.append(Dna[i][RandomNum:RandomNum + k])
    BestMotifs = motifs
    BestScore = Score(BestMotifs)
    for j in range(N):
        i = random.randint(0, t-1)
        #import pdb;pdb.set_trace()
        del(motifs[i])
        Profile = ProfileForming2( motifs)

        proList = ProfileRandomlyGenaratedKmer(Profile, Dna[i])

        motifs.insert(i, weight_choice(Dna[i],proList,k))
        #import pdb;pdb.set_trace()
        score1 = Score(motifs)
        if score1 < BestScore:
            BestMotifs = motifs
            BestScore = score1
        j += 1
    return BestMotifs, BestScore

def OuterLoopForRandomized(Dna, k, t, tim):
    k = int(k)
    t = int(t)
    minscore = float("inf")
    OuterBestMotif = []
    for i in range(tim):
        BestMotifs, score = RandomizedMotifSearch(Dna, k, t)

        if score < minscore:
            OuterBestMotif = BestMotifs
            minscore = score
        #import pdb;pdb.set_trace()
    return OuterBestMotif, minscore

def OuterLoopForGibbsSampler(Dna, k, t, N, Round):
    minscore = float("inf")
    OuterBestMotif = []
    for i in range(Round):
        BestMotifs, score = GibbsSampler(Dna, k, t, N)

        if score < minscore:
            OuterBestMotif = BestMotifs
            minscore = score
            #print ('job done by {} s.'.format(time.time() - starttime))
            #import pdb;pdb.set_trace()
    return OuterBestMotif, minscore

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
'''
fil = open(sys.argv[1],'r')
linlist = []
for lins in fil:
    linlist.append(lins.strip('\r\n'))
k, t = linlist[0].split(' ')

Dna = linlist[1:]

job, score = OuterLoopForRandomized(Dna, int(k), int(t), 1000)

print ('job done by {} s.'.format(time.time()-starttime))
import pdb;pdb.set_trace()
for i in job:   print (i)