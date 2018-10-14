import sys
import time
def PatternCount(aline, bline):
    count = 0
    ind = 0
    for ind in range(len(aline)):
        if aline[ind: ind + len(bline)] == bline:
            count = count + 1
            ind = ind + 1
    return count

def Frequentkmer(aline,kvalue, dic={}):
    count = 0
    ind = 0
    for ind in range(len(aline)-kvalue):
        shortstr = aline[ind:ind+kvalue]
        kmervalue = findtheB(aline, shortstr)
        dic[shortstr] = kmervalue
    return dic
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
def findplace(pattern, sequence):
    ind = 0
    indlst = []
    for ind in range(len(sequence)):
        if sequence[ind: ind + len(pattern)] == pattern:
            indlst.append(str(ind))
            ind = ind + 1
    #import pdb;pdb.set_trace()
    return ' '.join(indlst)
def PatternToNumber1(DNAsta):
    stalst, indlst = ['A','C','G','T'], [0,1,2,3]
    DNAstaindlst = []
    for ind in DNAsta:
        DNAstaindlst.append(indlst[stalst.index(ind)])
    #import pdb;pdb.set_trace()
    DNAstaindlst = [str(x) for x in DNAstaindlst]
    strnum = ''.join(DNAstaindlst)
    return int(strnum,4)
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

def ComputingFrequencies(Text, k):

    FrequencyArray = {i:0 for i in range(4**k)}
    for i in range(len(Text)-k+1):
        Pattern = Text[i:i+k]
        j = PatternToNumber2(Pattern)
        FrequencyArray[j] = FrequencyArray[j] + 1
    #import pdb;pdb.set_trace()
    '''
    numlst = [i for i in range(4**k)]
    keylst = [NumberToPattern(ind,k) for ind in numlst]
    for kii in list(set(keylst) - set(dic.keys())):
        dic[kii] = 0
    valuelst = [str(dic[kii]) for kii in keylst]
    totlst = [numlst,
              keylst,
              valuelst
            ]
    '''
    return FrequencyArray
def ComputingFrequenciesECOLI(Text,k,dic):
    for i in range(len(Text)-k):
        Pattern = Text[i:i+k]
        dic[Pattern] = dic[Pattern] + 1
    return dic
def PatternToNumber2(Pattern):
    stalst, indlst = ['A', 'C', 'G', 'T'], [0, 1, 2, 3]
    if len(Pattern) == 0:
        return 0
    symbol = Pattern[-1]
    Prefix = Pattern[0:-1]
    return 4 * PatternToNumber2(Prefix) + indlst[stalst.index(symbol)]
def LocateFre(Text, kvalue ,fre):
    totallist, kmerdic = ComputingFrequencies(Text, kvalue)
    return ' '.join(str(x) for x in kmerdic.keys() if kmerdic[x] == fre)
def ClumpFind(Genome, k, t, L):
    FrequentPatterns = []
    Clump = {i:0 for i in range(4**k)}
    for i in range(len(Genome)-L):
        Text = Genome[i:i+L]
        FrequencyArray = ComputingFrequencies(Text,k)
        for ind in range(4**k):
            if FrequencyArray[ind]>=t:
                Clump[ind] = 1
        for i in range(4**k):
            if Clump[i] == 1:
                Pattern = NumberToPattern(i,k)
                FrequentPatterns.append(Pattern)
        return FrequentPatterns
def BetterClumpFine(Genome,k,t,L):
    FrequentPatterns = []
    Clump = {i: 0 for i in range(4 ** k)}
    Text = Genome[0:L]
    FrequencyArray = {i: 0 for i in range(4 ** k)}
    FrequencyArray.update(ComputingFrequencies(Text,k))
    for i in range(4**k):
        if FrequencyArray[i]>=t:
            Clump[i]=1
    for i in range(len(Genome)-L):
        FirstPattern = Genome[i-1,i+k-1]
        index = PatternToNumber2(FirstPattern)
        FrequencyArray[index] = FrequencyArray[index] - 1
        LastPattern = Genome[i+L-k,i+L]
        index = PatternToNumber2(LastPattern)
        FrequencyArray[index] = FrequencyArray[index] - 1
        if FrequencyArray(index)>=t:
            Clump[index] = 1
        for i in range(4**k):
            if Clump(i) == 1:
                Pattern = NumberToPattern(i,k)
                FrequentPatterns.append(Pattern)
        return FrequentPatterns
def BetterClumpFind1(Genome,k,t,L,timestart):
    FrequentPatterns = []
    Clump = {i: 0 for i in range(4 ** k)}
    Text = Genome[0:L]
    FrequencyArray = {i: 0 for i in range(4 ** k)}
    FrequencyArray.update(ComputingFrequencies(Text, k))
    print('Time for first dictionary is {}'.format(time.time()-timestart))
    #import pdb;pdb.set_trace()
    for i in range(4 ** k):
        if FrequencyArray[i] >= t:
            Clump[i] = 1
    #import pdb;pdb.set_trace()
    for i in range(len(Genome) - L):
        FirstPattern = Genome[i - 1: i + k - 1]
        index = PatternToNumber2(FirstPattern)
        FrequencyArray[index] = FrequencyArray[index] - 1
        LastPattern = Genome[i + L - k: i + L]
        index = PatternToNumber2(LastPattern)
        FrequencyArray[index] = FrequencyArray[index] + 1
        if FrequencyArray[index] >= t:
            Clump[index] = 1
        #print i
    for ind in range(4 ** k):
        if Clump[ind] == 1:
            Pattern = NumberToPattern(ind, k)
            FrequentPatterns.append(Pattern)
    #import pdb;pdb.set_trace()
    return FrequentPatterns
def BetterClumpFindForEcoli(Genome,k,t,L,timestart):
    FrequentPatterns = []
    keyslst = set([Genome[i:i+9] for i in range(len(Genome)-9)])
    Zerolst = len(keyslst) * [0]
    #import pdb;pdb.set_trace()
    FrequencyArray = dict(zip(keyslst,Zerolst))
    print ('Time for first dictionary is {}'.format(time.time() - timestart))
    Text = Genome[0:L]
    FrequencyArray.update(ComputingFrequenciesECOLI(Text, k,FrequencyArray))
    print ('Time for update dictionary is {}'.format(time.time() - timestart))
    #import pdb;pdb.set_trace()
    # import pdb;pdb.set_trace()
    for ind in FrequencyArray.keys():
        if FrequencyArray[ind]>=t:
            FrequentPatterns.append(ind)
    #import pdb;pdb.set_trace()
    for i in range(len(Genome) - L-1):
        FirstPattern = Genome[i : i + k ]
        FrequencyArray[FirstPattern] = FrequencyArray[FirstPattern] - 1
        LastPattern = Genome[i + L - k+1: i + L+1]
        FrequencyArray[LastPattern] = FrequencyArray[LastPattern] + 1
        if FrequencyArray[LastPattern] >= t:
            FrequentPatterns.append(LastPattern)
        i += 1
        #print i
    # import pdb;pdb.set_trace()
    return FrequentPatterns

#Text1 = open(sys.argv[1],'r').readline()
Text1 = 'CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA'
kvalue = 9
limit = 3
Lvalue =500
#import pdb;pdb.set_trace()
timestart =time.time()
needlst = BetterClumpFindForEcoli(Text1,kvalue,limit,Lvalue,timestart)
print (len(set(needlst)))
print('{} seconds needed'.format(time.time()-timestart))

import pdb;pdb.set_trace()