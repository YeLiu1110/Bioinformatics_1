import time
import os
import sys
from numpy import *

def PatternToNumber1(DNAsta):
    stalst, indlst = ['A','C','G','T'], [0,1,2,3]
    DNAstaindlst = []
    for ind in DNAsta:
        DNAstaindlst.append(indlst[stalst.index(ind)])
    #import pdb;pdb.set_trace()
    DNAstaindlst = [str(x) for x in DNAstaindlst]
    strnum = ''.join(DNAstaindlst)
    return int(strnum,4) #base on 4
def Reversecomplement(DNAstring):
    DNAreverse = DNAstring[::-1]
    DNArecolst = []
    for ind in DNAreverse:
        if ind == 'G':
            DNArecolst.append('C')
        elif ind == 'C':
            DNArecolst.append('G')
        elif ind == 'A':
            DNArecolst.append('T')
        elif ind == 'T':
            DNArecolst.append('A')
        elif ind == 'g':
            DNArecolst.append('c')
        elif ind == 'c':
            DNArecolst.append('g')
        elif ind == 'a':
            DNArecolst.append('t')
        elif ind == 't':
            DNArecolst.append('a')
        else:
            DNArecolst.append(ind)
    return ''.join(DNArecolst)
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

def Skew(Genome):
    skew = [0]
    minSkew = [0]
    minSkewindex = []
    for i in range(1, len(Genome) + 1):
        # outlst.append(str(Skew(Genome,i)))
        if Genome[i - 1] == 'G':
            skew.append(skew[i - 1] + 1)
        elif Genome[i - 1] == 'C':
            skew.append(skew[i - 1] - 1)
        else:
            skew.append(skew[i - 1])
        if skew[i] <= min(minSkew):
            if skew[i] < min(minSkew):
                minSkew = [skew[i]]
                minSkewindex = [str(i)]
            else:
                minSkew.append(skew[i])
                minSkewindex.append(str(i))
        else:
            pass
    return minSkew, minSkewindex
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

def FindMismatch(String, Genome, d):
    '''
    Input: Strings Pattern and Text along with an integer d.
    Output: All starting positions where Pattern appears as a substring of Text with at most d mismatches.
    :param String: string
    :param Genome:
    :param d: integer
    :return: All starting positions where Pattern appears as a substring of Text with at most d mismatches.
    '''
    MismatchStart = []
    for i in range(0,(len(Genome)-len(String))+1):
        #import pdb;pdb.set_trace()
        if int(HammingDistance(String, Genome[i:i+len(String)])) <= int(d) :
            MismatchStart.append(str(i))
        else:
            pass
    return MismatchStart
def ApproximatePatternCount(bline, aline, d):
    count = 0
    ind = 0
    for ind in range(len(aline)-len(bline)+1):
        #import pdb;pdb.set_trace()
        if HammingDistance(aline[ind: ind + len(bline)], bline) <= int(d):
            count = count + 1
        else:
            pass
        ind = ind + 1
    return count
def ImmediateNeighbors(Patternlist):
    Neighborhood = []
    for pat in Patternlist:
        #import pdb;pdb.set_trace()
        for i in range(1,len(pat)):
            for x in set(['A','T','C','G']) - set(pat[i]):
                Neighbor = pat[1:].replace(pat[i],x)
                Neighborhood.append(Neighbor)
    #import pdb;pdb.set_trace()
    return Neighborhood

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
def IterativeNeighbors(Pattern, d):
    Neighborhood = set(Neighbors(Pattern,d))
    for j in range(1,d):
        for eah in Neighborhood:
            Neighborhood.add(eah)
            set(Neighborhood)
    return Neighborhood

def FrequentWordsWithMismatches(Genome, k, d):
    '''
    It first generates all neighbors (with up to d mismatches) for all k-mers in Text and combines them all into an array NeighborhoodArray. Note that a k-mer Pattern appears Countd(Text, Pattern) times in this array. The only thing left is to sort this array, count how many times each k-mer appears in the sorted array, and then return the k-mers that occur the maximum number of times.
    '''
    FrequentPatterns = set()
    Neighborhoods = []


    for i in range(len(Genome)-k+1):
        Neighborhoods.append(Neighbors(Genome[i:i+k],d))
    NeighborhoodArray = array(Neighborhoods)
    #import pdb;pdb.set_trace()
    IndexList = []
    Count = []
    for i in range(len(Neighborhoods)):
        Pattern = NeighborhoodArray[i]
        #import pdb;pdb.set_trace()
        for j in Pattern:
            IndexList.append(PatternToNumber1(j))
            Count.append(1)
    #import pdb;pdb.set_trace()
    SortedIndex = sort(IndexList)
    for i in range(len(Neighborhoods)*len(Neighborhoods[0])):
        if SortedIndex[i] == SortedIndex[i + 1]:
            Count[i + 1] = Count[i] + 1
    maxCount = max(Count)
    #import pdb;pdb.set_trace()
    for i in range(len(Neighborhoods)*len(Neighborhoods[0])):
        if Count[i] == maxCount:
            Pattern = NumberToPattern(SortedIndex[i], k)
            FrequentPatterns.add(Pattern)
    #import pdb;pdb.set_trace()
    return FrequentPatterns

def FrequentWordsWithMismatches2(Genome, k, d):
    '''
    Not Working!!!
    It first generates all neighbors (with up to d mismatches) for all k-mers in Text and combines them all into an array NeighborhoodArray. Note that a k-mer Pattern appears Countd(Text, Pattern) times in this array. The only thing left is to sort this array, count how many times each k-mer appears in the sorted array, and then return the k-mers that occur the maximum number of times.
    '''
    FrequentPatterns = set()
    Neighborhoods = []


    for i in range(len(Genome)-k):
        Neighborhoods.extend(list(Neighbors(Genome[i:i+k],d)))

    CountDic = {i:0 for i in Neighborhoods}
    SortedIndex = sort(Neighborhoods)
    for i in range(len(Neighborhoods)-1):
        if SortedIndex[i] == SortedIndex[i + 1]:
            CountDic[SortedIndex[i]] = CountDic[SortedIndex[i]] + 1
    KmerRevComDic = {}
    #import pdb;pdb.set_trace()
    for kmer in CountDic.keys():
        if Reversecomplement(kmer) in CountDic.keys():
        # import pdb;pdb.set_trace()
            KmerRevComDic.update({kmer: CountDic[kmer] + CountDic[Reversecomplement(kmer)]})
        else:
            KmerRevComDic.update({kmer: CountDic[kmer]})
    maxCountRe = max(KmerRevComDic.values())
    for kmer in KmerRevComDic.keys():
        if KmerRevComDic[kmer] == maxCountRe:
            FrequentPatterns.add(kmer)
    import pdb;pdb.set_trace()
    return FrequentPatterns

def FrequentWordsWithMismatches1(Genome, k, d):
    FrequentPattern = set()
    numDic = {i: 0 for i in range(4 ** k)}
    kmerDic = {}
    for num in numDic.keys():
        kmer = NumberToPattern(num, k)
        kmerDic[kmer] = ApproximatePatternCount(kmer, Genome, d)
    #import pdb;pdb.set_trace()
    KmerRevComDic = {}
    for kmer in kmerDic.keys():
        # import pdb;pdb.set_trace()
        KmerRevComDic.update({kmer: kmerDic[kmer] + kmerDic[Reversecomplement(kmer)]})

    maxCountRe = max(KmerRevComDic.values())
    for kmer in KmerRevComDic.keys():
        if KmerRevComDic[kmer] == maxCountRe:
            FrequentPattern.add(kmer)
    return FrequentPattern

def ComputingFrequenciesWithMismatches(Genome,k,d):
    FrequentPattern = set()
    FrequencyDic = {i:0 for i in range(4**k)}
    for i in range(len(Genome)-k):
        Pattern = Genome[i:i+k]
        Neighborhood = Neighbors(Pattern, d)
        for kmer in Neighborhood:
            j = PatternToNumber1(kmer)
            FrequencyDic[j] = FrequencyDic[j] + 1
    maxCount = max(FrequencyDic.values())
    for frekmer in FrequencyDic.keys():
        if FrequencyDic[frekmer] == maxCount:
            FrequentPattern.add(NumberToPattern(frekmer,k))

    return FrequentPattern
def ComputingFrequenciesWithMismatchesAndCom(Genome,k,d):
    FrequentPattern = set()

    KmerDic = {NumberToPattern(i,k):0 for i in range(4**k)}
    for i in range(len(Genome)-k):
        Pattern = Genome[i:i+k]
        Neighborhood = Neighbors(Pattern, d)
        for kmer in Neighborhood:
            KmerDic[kmer] = KmerDic[kmer] + 1

    KmerRevComDic = {}
    for kmer in KmerDic.keys():
        #import pdb;pdb.set_trace()
        KmerRevComDic.update({kmer:KmerDic[kmer] + KmerDic[Reversecomplement(kmer)]})

    maxCountRe = max(KmerRevComDic.values())
    for kmer in KmerRevComDic.keys():
        if KmerRevComDic[kmer] == maxCountRe:
            FrequentPattern.add(kmer)
    #import pdb;pdb.set_trace()
    return FrequentPattern

starttime = time.time()

fil = open(sys.argv[1],'r').readlines()
Stringpattern , Genome, d = fil[0].strip('\n'), fil[1].strip('\n'), fil[2].strip('\n')
patterncount = ApproximatePatternCount(Stringpattern, Genome, d)

fil = open(sys.argv[1],'r').readlines()
linlist = []
for lins in fil[1:]:
    linlist.append(lins.strip('\r\n'))
WholeGenome = ''.join(linlist).strip('\n')
#import pdb;pdb.set_trace()
SkewPosition = Skew(WholeGenome)[1]
Genome = WholeGenome[int(SkewPosition[0])-500:int(SkewPosition[-1])+500]
#import pdb;pdb.set_trace()
k = 9
d = 1

job = ComputingFrequenciesWithMismatchesAndCom(Genome,k,d)


job = Neighbors('ACGT',3)
print (job)
print ('job done by {} s.'.format(time.time()-starttime))
import pdb;pdb.set_trace()