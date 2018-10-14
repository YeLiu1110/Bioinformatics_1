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

def LinearSpectrum(Peptide, Massdic):
    PrefixMass = [0]
    for i in range(len(Peptide)):
        PrefixMass.append(PrefixMass[i] + Massdic[Peptide[i]])
    Linear = [Massdic[Peptide[i]] for i in range(len(Peptide))]
    for i in range(len(Peptide)):
        for num in range(i+1,len(Peptide)):
            Linear.append(PrefixMass[num+1] - PrefixMass[i])
    Linear.insert(0,0)
    #import pdb;pdb.set_trace()
    return sorted(Linear)

def Cyclospectrum(Peptide, Massdic):
    PrefixMass = [0]
    for i in range(len(Peptide)):
        PrefixMass.append(PrefixMass[i] + Massdic[Peptide[i]])
    peptideMass = PrefixMass[-1]
    #import pdb;pdb.set_trace()
    Cyclic = [Massdic[Peptide[i]] for i in range(len(Peptide))]
    for i in range(len(Peptide)):
        for num in range(i + 1, len(Peptide)):
            Cyclic.append(PrefixMass[num+1] - PrefixMass[i])
            if (i > 0) and (num < len(Peptide)):
                Cyclic.append(peptideMass - PrefixMass[num] + PrefixMass[i])

    Cyclic.insert(0, 0)
    #import pdb;pdb.set_trace()
    return sorted(Cyclic)

def CyclopeptideSequencing(Spectrum, Massdic):
    def Consistent(shortlist,longlist):
        shortset = set(shortlist)
        for i in shortset:
            if i not in longlist:
                return False
            #elif shortlist.count(i) > longlist(i):
                #return False
        return True
    def Mass(Peptide,Massdic):
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
            elif Mass(Pep,Massdic) >= Spectrum[-1]:
                if Cyclospectrum(Pep,Massdic) == Spectrum:
                    output.append(Pep)
                PeptideRemove.append(Pep)
            elif Consistent(LinearSpectrum(Pep,Massdic),Spectrum) == False:
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

starttime = time.time()
'''
fil = open(sys.argv[1],'r')
linlist = []
for lins in fil:
    linlist.append(lins.strip('\r\n'))
DNAstring = ''.join(linlist)
print ('Bacteria is read in {} s.'.format(time.time()-starttime))
#results = linlist[2:]

#DNAstring = 'ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA'
#Peptide = 'MA'

fil = open(sys.argv[1],'r')
linlist = []
for lins in fil:
    linlist.append(lins.strip('\r\n'))
#Peptide = linlist[0]
Spec = linlist[0].split(' ')
Spectrum = []
for i in Spec:
    Spectrum.append(int(i))
    
'''
Codondic = CodonDic()
Massdic = MassDic()

global NumList
NumList = sorted(list(set(Massdic.values())))
job = CyclopeptideSequencing([0,71,101,113,131,184,202,214,232,285,303,315,345,416],Massdic)
jobs = []
spe = ['AVQ','QCV','TCE','TVQ','TCQ','ETC']
for i in spe:
    jobs.append(LinearSpectrum(i,Massdic))


#import pdb;pdb.set_trace()
print ('job done by {} s.'.format(time.time()-starttime))
'''
cwd = os.getcwd()

txtFile1 = cwd + '/results.txt'
txtFile2 = cwd + '/results1.txt'
with open(txtFile1,'w') as f:
    for ind,i in enumerate(jobs):
        f.write('{} - > {} '.format(spe[ind],i))
        f.write('\n')

f.close()
with open(txtFile2,'w') as f:
    for i in sorted(job):
        f.write(''.join(i))
        f.write('\n')
f.close()
'''
import pdb;pdb.set_trace()
for i in job:   print (i)