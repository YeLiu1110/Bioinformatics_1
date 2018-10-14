import time
import os
import sys
import pandas as pd
import math
import random
import numpy as np

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

def Translate(RNAstring,Codondic):
    Proteinstring = []
    i = 0
    while i in range(len(RNAstring)-len(RNAstring)%3):
        Proteinstring.append(Codondic[RNAstring[i:i+3]])
        i = i +3
    return ''.join(Proteinstring)

def CodonDic():
    cwd = os.getcwd()
    DicFil = cwd + '/RNA_codon_table_1.txt'
    Codondic = {}
    fil = open(DicFil,'r')
    for lins in fil:
        Codondic.update({lins.strip('\r\n').split(' ')[0]:lins.strip('\r\n').split(' ')[1]})
    for ke,va in Codondic.items():
        if Codondic[ke] == '':
            Codondic[ke] = '*'
    return Codondic

def MassDic():
    cwd = os.getcwd()
    DicFil = cwd + '/integer_mass_table.txt'
    Masdic = {}
    fil = open(DicFil, 'r')
    for lins in fil:
        Masdic.update({lins.strip('\r\n').split(' ')[0]: int(lins.strip('\r\n').split(' ')[1])})
    return Masdic

global Massdic
Codondic = CodonDic()
Massdic = MassDic()

def PeptideEncoding(DNAstring, Peptide, Codondic):
    def DNAtoRNA(DNAstring):
        return DNAstring.replace('T','U')
    Enstring = []
    PepLen = len(Peptide)
    DNATranslateString = [Translate(DNAtoRNA(DNAstring)[i:],Codondic) for i in range(3)]
    ReverDNAstring = Reversecomplement(DNAstring)
    DNAreverTransString = [Translate(DNAtoRNA(ReverDNAstring)[i:],Codondic) for i in range(3)]
    for tim in range(3):

        for i in range(len(DNATranslateString[tim])):
            if DNATranslateString[tim][i:i+PepLen] == Peptide:
                Enstring.append(DNAstring[i*3+tim:(i+PepLen)*3+tim])
        #import pdb;pdb.set_trace()

        for pi in range(len(DNAreverTransString[tim])):
            if DNAreverTransString[tim][pi:pi+PepLen] == Peptide:
                Enstring.append(Reversecomplement(ReverDNAstring[(pi*3+tim):(pi+PepLen)*3+tim]))
    #import pdb;pdb.set_trace()

    return Enstring

def SubpeptideCount(CyclicPeptideLength):
    le = CyclicPeptideLength
    return le*(le-1)

def LinearSpectrum(Peptide):
    PrefixMass = [0]
    for i in range(len(Peptide)):
        PrefixMass.append(PrefixMass[i] + Massdic[Peptide[i]])
    Linear = [Massdic[Peptide[i]] for i in range(len(Peptide))]
    for i in range(len(Peptide)):
        for num in range(i+1,len(Peptide)):
            Linear.append(PrefixMass[num+1] - PrefixMass[i])
    Linear.insert(0,0)
    return sorted(Linear)

def LinearSpectrum1(Peptide):
    PrefixMass = [0]
    for i in range(len(Peptide)):
        PrefixMass.append(PrefixMass[i] + Peptide[i])
    Linear = [i for i in Peptide]
    for i in range(len(Peptide)):
        for num in range(i+1,len(Peptide)):
            #import pdb;pdb.set_trace()
            Linear.append(PrefixMass[num+1] - PrefixMass[i])
    Linear.insert(0,0)
    return sorted(Linear)

def Cyclospectrum(Peptide):
    PrefixMass = [0]
    for i in range(len(Peptide)):
        PrefixMass.append(PrefixMass[i] + Massdic[Peptide[i]])
    peptideMass = PrefixMass[-1]
    Cyclic = [Massdic[Peptide[i]] for i in range(len(Peptide))]
    for i in range(len(Peptide)):
        for num in range(i + 1, len(Peptide)):
            Cyclic.append(PrefixMass[num+1] - PrefixMass[i])
            if (i > 0) and (num < len(Peptide)):
                Cyclic.append(peptideMass - PrefixMass[num] + PrefixMass[i])
    Cyclic.insert(0, 0)
    return sorted(Cyclic)

def CyclopeptideSequencing(Spectrum):
    def Consistent(shortlist,longlist):
        shortset = set(shortlist)
        for i in shortset:
            if i not in longlist:
                return False
            #elif shortlist.count(i) > longlist(i):
                #return False
        return True
    def Mass(Peptide):
        #import pdb;pdb.set_trace()
        return sum([Massdic[i] for i in Peptide])
    output = []
    AClist = list(Massdic.keys())
    Peptides = [[i] for i in AClist if Massdic[i] in Spectrum]
    PeptideRemove = []
    #import pdb;pdb.set_trace()
    round = 1

    while Peptides:

        lenth = len(Peptides)
        for i in range(lenth):
            for j in AClist:
                Peptides.append(Peptides[i]+ [j])
        #import pdb;pdb.set_trace()
        for Pep in Peptides:
            if len(Pep) == round:
                PeptideRemove.append(Pep)
            elif Mass(Pep) >= Spectrum[-1]:
                if Cyclospectrum(Pep) == Spectrum:
                    output.append(Pep)
                PeptideRemove.append(Pep)
            elif Consistent(LinearSpectrum(Pep),Spectrum) == False:
                PeptideRemove.append(Pep)
        #import pdb;pdb.set_trace()
        for i in PeptideRemove:
            Peptides.remove(i)
        PeptideRemove = []
        round = round + 1
    #import pdb;pdb.set_trace()
    outputnum = []
    for i in output:
        outputnum.append([str(Massdic[j]) for j in i])
    outputnumset = []
    for i in outputnum:
        if i not in outputnumset:
            outputnumset.append(i)
    return output

def CountingPeptides(n):
    '''
    Counting Peptides with Given Mass Problem: Compute the number of peptides of given mass.
     Input: An integer m.
     Output: The number of linear peptides having integer mass m.
    Exercise Break: Solve the Counting Peptides with Given Mass Problem. Recall that we assume that peptides are formed
    from the following 18 amino acid masses:
    G	A	S	P	V	T	C	I/L	N	D	K/Q	E	M	H	F	R	Y	W
    57	71	87	97	99	101	103	113	114	115	128	129	131	137	147	156	163	186
    Sample Input:
    1024
    Sample Output:
    14712706211
    '''
    massTalbe = NumList
    m = len(massTalbe)
    table = [0] * (n+1)
    table[0] = 1

    for i in range(n+1):
        currSum = 0
        for j in range(m):
            if i - massTalbe[j] >= 0:
                currSum += table[i-massTalbe[j]]
        table[i] += currSum
    import pdb;pdb.set_trace()
    return table[n]

def Score(Peptide, Spectrum):
    '''
    Cyclopeptide Scoring Problem: Compute the score of a cyclic peptide against a spectrum.
     Input: An amino acid string Peptide and a collection of integers Spectrum.
     Output: The score of Peptide against Spectrum, Score(Peptide, Spectrum).
    '''
    PeptideSpectrum = Cyclospectrum(Peptide)
    setPepSpectrum = set(PeptideSpectrum)
    score = 0
    for i in setPepSpectrum:
        if PeptideSpectrum.count(i) >= Spectrum.count(i):
            score += Spectrum.count(i)
        else:
            score += PeptideSpectrum.count(i)
    return score

def LinearScore(Peptide, Spectrum):
    PeptideSpectrum = LinearSpectrum(Peptide)
    setPepSpectrum = set(PeptideSpectrum)
    score = 0
    for i in setPepSpectrum:
        if PeptideSpectrum.count(i) >= Spectrum.count(i):
            score += Spectrum.count(i)
        else :
            score += PeptideSpectrum.count(i)

    return score

def LinearScore1(Peptide, Spectrum):
    PeptideSpectrum = LinearSpectrum1(Peptide)
    setPepSpectrum = set(PeptideSpectrum)
    score = 0
    for i in setPepSpectrum:
        if PeptideSpectrum.count(i) >= Spectrum.count(i):
            score += Spectrum.count(i)
        else :
            score += PeptideSpectrum.count(i)
    return score

def Trim(Leaderboard, Spectrum, N):

    ScoreList = []
    for i in range(len(Leaderboard)):
        Peptide = Leaderboard[i]
        ScoreList.append(LinearScore(Peptide,Spectrum))
    #import pdb;pdb.set_trace()
    SortPepList = [i for j,i in sorted(zip(ScoreList, Leaderboard),reverse = True)]
    SortScoreList = [j for j,i in sorted(zip(ScoreList, Leaderboard),reverse = True)]
    #import pdb;pdb.set_trace()
    for i in range(N,len(Leaderboard)):
        if SortScoreList[i] < SortScoreList[N-1]:
            return SortPepList[:i]
    return SortPepList

def Trim1(Leaderboard, Spectrum, N):
    ScoreList = []
    for i in range(len(Leaderboard)):
        Peptide = Leaderboard[i]
        ScoreList.append(LinearScore1(Peptide,Spectrum))
    #import pdb;pdb.set_trace()
    SortPepList = [i for j,i in sorted(zip(ScoreList, Leaderboard),reverse = True)]
    SortScoreList = [j for j,i in sorted(zip(ScoreList, Leaderboard),reverse = True)]
    #import pdb;pdb.set_trace()
    for i in range(N,len(Leaderboard)):
        if SortScoreList[i] < SortScoreList[N-1]:
            return SortPepList[:i]
    return SortPepList

def LeaderboardCyclopeptideSequencing(Spectrum, N, AClist = list(Massdic.keys())):

    def Mass(Peptide):
        return sum([Massdic[i] for i in Peptide])

    LeaderBoard = [[i] for i in AClist if Massdic[i] in Spectrum]
    LeaderRemove = []
    LeaderPeptide = ''
    round = 1
    GatherPep = []
    while LeaderBoard:
        lenth = len(LeaderBoard)
        for i in range(lenth):
            for j in AClist:
                LeaderBoard.append(LeaderBoard[i]+ [j])
        for Pep in LeaderBoard:
            if len(Pep) == round:
                LeaderRemove.append(Pep)
            elif Mass(Pep) == Spectrum[-1]:
                if Score(Pep, Spectrum) > Score(LeaderPeptide, Spectrum):
                    LeaderPeptide = Pep
                #if Score(Pep, Spectrum) == 83:
                    #GatherPep.append(Pep)
                #LeaderRemove.append(Pep)
            elif Mass(Pep) > Spectrum[-1]:
                LeaderRemove.append(Pep)
        for i in LeaderRemove:
            LeaderBoard.remove(i)
        LeaderBoard = Trim(LeaderBoard,Spectrum,N)
        LeaderRemove = []
        round += 1
    return [Massdic[i] for i in LeaderPeptide]

def LeaderboardCyclopeptideSequencing1(Spectrum, N, AClist = list(Massdic.values())):
    def Cyclospectrum1(Peptide):
        PrefixMass = [0]
        for i in range(len(Peptide)):
            PrefixMass.append(PrefixMass[i] + Peptide[i])
        peptideMass = PrefixMass[-1]
        Cyclic = [i for i in Peptide]
        for i in range(len(Peptide)):
            for num in range(i+1,len(Peptide)):
                Cyclic.append(PrefixMass[num+1] - PrefixMass[i])
                if (i > 0) and (num < len(Peptide)):
                    Cyclic.append(peptideMass- PrefixMass[num] + PrefixMass[i])
        Cyclic.insert(0, 0)
        return sorted(Cyclic)

    def Score1(Peptide, Spectrum):
        PeptideSpectrum = Cyclospectrum1(Peptide)
        setPepSpectrum = set(PeptideSpectrum)
        score = 0
        for i in setPepSpectrum:
            if PeptideSpectrum.count(i) >= Spectrum.count(i):
                score += Spectrum.count(i)
            else:
                score += PeptideSpectrum.count(i)
        return score
    def Mass1(Peptide):
        return sum(Peptide)

    LeaderBoard = [[i] for i in AClist if i in Spectrum]
    LeaderRemove = []
    round = 1
    GatherPepDic = {}
    while LeaderBoard:
        lenth = len(LeaderBoard)
        for i in range(lenth):
            for j in AClist:
                LeaderBoard.append(LeaderBoard[i]+ [j])

        for Pep in LeaderBoard:
            if len(Pep) == round:
                LeaderRemove.append(Pep)
            elif Mass1(Pep) == Spectrum[-1]:
                GatherPepDic.update({'-'.join(map(str,Pep)):Score1(Pep, Spectrum)})
                #if Score1(Pep, Spectrum) > Score1(LeaderPeptide, Spectrum):
                    #LeaderPeptide = Pep
                #LeaderRemove.append(Pep)
            elif Mass1(Pep) > Spectrum[-1]:
                LeaderRemove.append(Pep)
        for i in LeaderRemove:
            LeaderBoard.remove(i)
        LeaderBoard = Trim1(LeaderBoard, Spectrum, N)
        LeaderRemove = []
        round += 1
    #import pdb;pdb.set_trace()
    LeaderPepList = [i for i in GatherPepDic if GatherPepDic[i] == max(GatherPepDic.values())]
    return LeaderPepList

def Convolution(spectrum):
    '''
    Spectral Convolution Problem: Compute the convolution of a spectrum.
     Input: A collection of integers Spectrum.
     Output: The list of elements in the convolution of Spectrum. If an element has multiplicity k, it should appear exactly k times;
    you may return the elements in any order.
    '''
    convolut = []
    for i in spectrum:
        for j in spectrum:
            if j - i > 0:
                convolut.append(j-i)
    return convolut

def ConvolutionCyclopeptideSequencing(spectrum, M, N):
    '''
    We now have the outline for a new cyclopeptide sequencing algorithm. Given an experimental spectrum, we first compute the convolution of an experimental spectrum.
    We then select the M most frequent elements between 57 and 200 in the convolution to form an extended alphabet of candidate amino acid masses.
    In order to be fair, we should include the top M elements of the convolution "with ties".
    Finally, we run the algorithm LeaderboardCyclopeptideSequencing, where the amino acid masses are restricted to this alphabet.
    We call this algorithm ConvolutionCyclopeptideSequencing.
    '''
    Spectrum = sorted(spectrum)
    Convo = sorted(Convolution(Spectrum))
    Convolu = []
    for i in Convo:
        if (i>=57) and (i<=200):
            Convolu.append(i)
    #import pdb;pdb.set_trace()
    ConvoluCount = [Convolu.count(i) for i in set(Convolu)]
    SortConvolu = [i for j,i in sorted(zip(ConvoluCount,set(Convolu)),reverse = True)]
    SortConvoluCount = [j for j,i in sorted(zip(ConvoluCount,set(Convolu)),reverse = True)]
    ConvoluPep = SortConvolu
    #import pdb;pdb.set_trace()
    for i in range(M,len(SortConvolu)):
        if SortConvoluCount[i] < SortConvoluCount[i-1]:
            ConvoluPep = SortConvolu[:i]
            break
    import pdb;pdb.set_trace()
    LeaderPepList = LeaderboardCyclopeptideSequencing1(sorted(spectrum),N, ConvoluPep)
    #import pdb;pdb.set_trace()
    return LeaderPepList

starttime = time.time()


fil = open(sys.argv[1],'r')
linlist = []
for lins in fil:
    linlist.append(lins.strip('\r\n'))
#M = linlist[0]
#N = linlist[1]
Spec = linlist[0].split(' ')
Spectrum = []
for i in Spec:
    Spectrum.append(int((float(i) - 1.007) * (1 - 0.0004522)))



'''
Spectrum = [57, 57, 71, 99, 129, 137, 170, 186, 194, 208, 228, 265, 285, 299, 307, 323, 356, 364, 394, 422, 493]
M = 20
N = 60
'''
job = ConvolutionCyclopeptideSequencing(Spectrum, 20, 1000)

print ('job done by {} s.'.format(time.time()-starttime))

cwd = os.getcwd()
#import pdb;pdb.set_trace()
txtFile1 = cwd + '/results.txt'
#txtFile2 = cwd + '/results1.txt'
with open(txtFile1,'w') as f:
    f.write(' '.join(job))

f.close()
'''
with open(txtFile2,'w') as f:
    for i in sorted(job):
        f.write(''.join(i))
        f.write('\n')
f.close()
'''
import pdb;pdb.set_trace()
for i in job:   print (i)