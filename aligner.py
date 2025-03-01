#!/usr/bin/env python
# coding: utf-8

# This was orignally the working notebook created on July 3, 2024. <br>
# I reorganized the code into two other notebooks: aligner_helper (with output evaluation and file generation scripts) and alignment_outout (where I run the script below). <br>
# Reorganization date Jul 18, 2024. 

# In[13]:


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pysam
import pybedtools
import random
import pyBigWig
import enum
import os
import argparse
import sys



def delta(symbol1, symbol2):
    '''
    input: two symbols (nucleotides)
    
        if symbols are the same return 0, else 1
        
    output: score
    '''
    return 0 if symbol1 == symbol2 else 1


# ##### The following functions are only used for NW alignment

# In[16]:


def levenshteinDistance(a,b):
    '''
    input: two sequences a and b
    
        Recursive function to compute the total distance 
        (i.e., num of changes to get from a to b)
        
    output: lowest distance of three possibilities
    '''
    if(len(a) == 0):
        return len(b)
    elif(len(b) == 0):
        return len(a)
    else:
        deletionDistance = levenshteinDistance(a[1:], b) + 1
        insertionDistance= levenshteinDistance(a, b[1:]) + 1
        substituteDistance=levenshteinDistance(a[1:],b[1:]) + delta(a[0], b[0])
        return min(deletionDistance, insertionDistance, substituteDistance)


# In[17]:


def levenshteinMemo(a,b,memo):
    '''
    input: two sequences a and b and memoization matrix
    
        Dynamic programming implementation of Levenshtein distance. 
        Recursive function that builds memoization matrix
        and returns the distance of the last cell of the matrix
        
    output: last item of the memoization dictionary
    '''
    if((a,b) not in memo):
        if(len(a) == 0):
            memo[(a,b)] = len(b)
        elif(len(b) == 0):
            memo[(a,b)] = len(a)
        else:
            delDist = levenshteinMemo(a[1:], b, memo) +1
            insertDist = levenshteinMemo(a, b[1:], memo) +1
            subsDist =   levenshteinMemo(a[1:], b[1:], memo) \
                       + levenshtein(a[0], b[0])
            memo[(a,b)] = min(delDist, insertDist, subsDist)
    return memo[(a,b)]


# In[18]:


def nwDist(a,b, gapPenalty=1):
    '''
    input: sequences a and b and gap penalty score
    
        Needleman-Wunsch algorithm that stores memoization matrix as an array
        rather than a dictionary
        
    output: memoization matrix
    '''
    memoAr = np.zeros((len(a)+1, len(b)+1))
    for i in range(1,len(a)+1):
        memoAr[i,0] = memoAr[i-1,0] + gapPenalty
    for j in range(1,len(b)+1):
        memoAr[0,j] = memoAr[0, j-1] + gapPenalty
    for i in range(1, len(a)+1):
        for j in range(1, len(b)+1):
            memoAr[i,j] = min(memoAr[i-1, j-1] + delta(a[i-1], b[j-1]), 
                              memoAr[i-1, j] + gapPenalty, 
                              memoAr[i, j-1] + gapPenalty)
    return memoAr


# In[19]:


def nwAlign(a,b, distMat, gapPenalty=1):
    '''
    input: sequences a and b, memoization matrix and gap penalty cost
    
        Traceback function that starts from the last cell of the matrix
        and traces all deletion/insertion/substutution decisons made
        
    output: sequence alignment
    '''
    alignA = ""
    alignB = ""
    i = len(a)
    j = len(b)
    while(i>0 or j>0):
        if(i>0 and j>0 and distMat[i,j] == distMat[i-1, j-1] + delta(a[i-1], b[j-1])):
            alignA = a[i-1] + alignA
            alignB = b[j-1] + alignB
            i -= 1
            j -= 1
        elif (i>0 and distMat[i,j] == distMat[i-1,j] + gapPenalty):
            alignA = a[i-1] + alignA
            alignB = "-" + alignB
            i -= 1
        else:
            alignA = "-" + alignA
            alignB = b[j-1] + alignB
            j -= 1
    return (alignA, alignB)


# ##### These functions are used in the profile aligner

# In[20]:


def gotohDist(a, b, u, v, delta):
    '''
    input: sequences a and b, u gap extension cost, v gap insertion cost, and loss function
    
        Affine gap penalty calculation using Gotoh's algorithm
        Recursively fills out three matrices where D keeps track of overall best score
        P keeps track of best scores under the constraint of ending seq b with a gap
        Q keeps track of best scores under the constraint of ending seq a with a gap
        
    output: three memoization matrices
            D - NW matrix, 
            P - matrix when alignment ends in a gap in b
            Q - matrix when alignment ends in a gap in a
    '''
    D = np.zeros((len(a)+1, len(b)+1))
    P = np.zeros((len(a)+1, len(b)+1))
    Q = np.zeros((len(a)+1, len(b)+1))
    for i in range(1,len(a)+1):
        D[i,0] = i*u+v
        Q[i,0] = 10**5
    for j in range(1,len(b)+1):
        D[0,j] = j*u+v
        P[0,j] = 10**5
    for i in range(1, len(a)+1):
        for j in range(1, len(b)+1):
            P[i,j] = min(D[i-1,j] + u+v,
                         P[i-1,j]  + u)
            Q[i,j] = min(D[i,j-1] + u+v,
                         Q[i,j-1]  + u)
            D[i,j] = min(D[i-1,j-1] + delta(a[i-1], b[j-1]),
                         P[i,j],
                         Q[i,j])
    return (D,P,Q) 


# In[12]:


def gotohAlign(a, b, u, v, delta, D, P, Q):
    '''
    input: sequences a and b, 
            u gap extension cost, v gap insertion cost,
            loss function, 
            D, P, and Q matrices

            Traceback of three Gotoh matrices. 
            Generates TD, TQ, and TP matrices by referencing D, Q, and P
            Uses traceback matrices to generate alignment
            
    Output: sequence alignment 
            TQ, TP, and TD traceback matrices
    '''
    moves=enum.IntEnum('moves', ["DIAG_D", "UP_D", "LEFT_D", "UP_P", 
                                 "LEFT_Q", "DOT_P", "DOT_Q"])
    TD = np.zeros((len(a)+1, len(b)+1))
    TP = np.zeros((len(a)+1, len(b)+1))
    TQ = np.zeros((len(a)+1, len(b)+1))
    TD[:,0] = moves.UP_D
    TD[0,:] = moves.LEFT_D
    TP[:,0] = moves.UP_P
    TQ[0,:] = moves.LEFT_Q
    
    for i in range(1, len(a)+1):
        for j in range(1,len(b)+1):
            #set TD
            if(D[i,j] == D[i-1, j-1]+delta(a[i-1], b[j-1])):
                TD[i,j] = moves.DIAG_D
            elif D[i,j] == Q[i,j]:
                TD[i,j] = moves.DOT_Q
            elif D[i,j] == P[i,j]:
                TD[i,j] = moves.DOT_P
            else:
                assert False
            
            #set TP
            if(P[i,j] == P[i-1,j] + u):
                TP[i,j] = moves.UP_P
            elif(P[i,j] == D[i-1,j] + u+v):
                TP[i,j] = moves.UP_D
            else:
                assert False
                
            #set TQ
            if(Q[i,j] == Q[i,j-1] + u):
                TQ[i,j] = moves.LEFT_Q
            elif(Q[i,j] == D[i,j-1] + u+v):
                TQ[i,j] = moves.LEFT_D
            else:

                print(Q[i,j], Q[i,j-1], u)
                print(D[i,j-1]+u+v)
                assert False
    #And now do the traceback. 
    alignA = []
    alignB = []
    i = len(a)
    j = len(b)
    curMat = TD
    while(i>0 or j>0):
        if(curMat[i,j] == moves.DIAG_D.value):
            alignA = [a[i-1]] + alignA
            alignB = [b[j-1]] + alignB
            i -= 1
            j -= 1
            curMat = TD
    
        elif(curMat[i,j] == moves.UP_D.value):
            alignA = [a[i-1]] + alignA
            alignB = ["-"] + alignB
            i -= 1
            curMat = TD
        
        elif(curMat[i,j] == moves.LEFT_D.value):
            alignA = ['-'] + alignA
            alignB = [b[j-1]] + alignB
            j -= 1
            curMat = TD
    
        elif(curMat[i,j] == moves.UP_P.value):
            alignA = [a[i-1]] + alignA
            alignB = ['-'] + alignB
            i -= 1
            curMat = TP
            
        elif(curMat[i,j] == moves.LEFT_Q.value):
            alignA = ['-'] + alignA
            alignB = [b[j-1]] + alignB
            j -= 1
            curMat = TQ
    
        elif(curMat[i,j] == moves.DOT_P.value):
            curMat = TP
        elif(curMat[i,j] == moves.DOT_Q.value):
            curMat=TQ
        
        else:
            print(curMat[i,j])
            print(i)
            print(j)
            
            assert False
            
        
        
    return alignA, alignB, TP, TQ, TD


# In[11]:


def gotoh(a, b, u, v, delta):
    '''
    input: sequences a and b, 
            u gap extension cost, v gap insertion cost,
            loss function

            Puts together previous helper functions to generate alignment and 

    Output: alignments a and b, and total loss
    '''
    D, P, Q = gotohDist(a, b, u, v, delta)
    alnA, alnB, _, _, _ = gotohAlign(a, b, u, v, delta, D, P, Q)
    return (alnA, alnB, D[-1,-1])


# In[12]:


def printAln(alnA, alnB, delta, stringer=lambda x:x): 
    '''
    input: alignments for both sequences, and loss function
    
        Formats the alignment as follows:
        seq a alignment
        loss
        seq b alignment
        
    output: formatted alignment
    '''
    print("".join([stringer(x) for x in alnA]))
    alnStr = ""
    for i in range(len(alnA)):
        if(alnA[i] == alnB[i]):
            alnStr += ' '
        elif(alnA[i] == '-' or alnB[i] == '-'):
            alnStr += '|'
        else:
            alnStr += str(delta(alnA[i], alnB[i]))[0]
    print(alnStr)
    print("".join([stringer(x) for x in alnB]))



# In[13]:


# I updated this function in order to compute relative differences instead of absolute (Jul 29, 2024)

def seqProfile(alpha):
    '''
    input: weight parameter
    
        Computes the sum of sequence and profile losses for each element
        
    output: loss for the given element
    '''
    def comparator(elem1, elem2):
        deltaSequence = 0 if elem1.sequence == elem2.sequence else alpha
        deltaProfile = (1-alpha) * abs(elem1.profile - elem2.profile)
        return deltaSequence + deltaProfile
        
    return comparator


# In[14]:


class Element:
    def __init__(self, sequence, profile):
        self.sequence = sequence
        self.profile = profile
    def __str__(self):
        return "{0:s} {1:f}".format(self.sequence, self.profile)


# In[1]:


#Jul 17: updated the function to take multiple regionIdx
def getSequences(regionIdx, normalize):
    '''
    input: genome regionIdx corresponding to bed file, profile loss normalization function
    
        Retrieves sequences from FASTA files
        and TF binding predictions from BW files
        
    output: nucleotide sequence of interest and corresponding TF binding predictions
    
    '''
    elements = dict()
    region_coords = []
    for genome_index, organism in enumerate(genomes):
        orgBed = pybedtools.BedTool(coordsFnames[organism])
        orgRegion = orgBed[regionIdx[genome_index]]
        orgGenome = pysam.FastaFile(genomeFnames[organism])
        orgCenter = (orgRegion.start + orgRegion.end)//2
        orgStart = orgCenter - 500
        orgEnd = orgStart + 1000
        orgSequence = orgGenome.fetch(orgRegion.chrom, orgStart, orgEnd)
        bw = pyBigWig.open(contributionsFnames[organism])
        orgProfile = bw.values(orgRegion.chrom, orgStart, orgEnd)
        bw.close()
        
        profileAr = normalize(np.array(orgProfile))
        elements[organism] = []
        for i in range(len(orgSequence)):
            elements[organism].append(Element(orgSequence[i].upper(), profileAr[i]))
        region_coords.append(f'{organism}: {orgRegion.chrom} {orgStart}-{orgEnd}')
    
    return elements, region_coords


# In[16]:


def printFolded(strings, foldLen):
    '''
    input: two sequence alignments and profile loss, 
        fold length (i.e., length of each line)

        Formats the output 

    output: formatted alignment 
    '''
    def chunker(string):
        for i in range(0,len(string), foldLen):
            yield string[i:i+foldLen]
    
    chunked = [list(chunker(s)) for s in strings]
    for c in range(len(chunked[0])):
        for s in range(len(chunked)):
            print(chunked[s][c])
        print()
        print()


# In[17]:


def printElemAln(alnA, alnB, delta, stringer=lambda x:x, foldLen=80):
    
    alnStr = ""
    aProfile = ""
    bProfile = ""
    alnStats = []
    aStats = []
    bStats = []
    #Build up arrays of the delta and profiles, so that we can display 
    #that information in a normalized format. 
    for i in range(len(alnA)):
        match (alnA[i], alnB[i]):
            case ('-', bmatch):
                bStats.append(bmatch.profile)
            case (amatch, '-'):
                aStats.append(amatch.profile)
            case (amatch, bmatch):
                aStats.append(amatch.profile)
                bStats.append(bmatch.profile)
                alnStats.append(delta(alnA[i], alnB[i]))

    if len(alnStats) > 0:
        alnQuant = np.quantile(np.array(alnStats), [x/8 for x in range(1,9)])
        alnQuant = np.linspace(np.min(alnStats), np.max(alnStats), num=8)
    else: 
        alnQuant = []
        
    aQuant = np.quantile(np.array(aStats), [x/8 for x in range(1,9)])
    bQuant = np.quantile(np.array(bStats), [x/8 for x in range(1,9)])
    
    aQuant = np.linspace(np.min(aStats), np.max(aStats), num=8)
    bQuant = np.linspace(np.min(bStats), np.max(bStats), num=8)
    
    
    for i in range(len(alnA)):
        match (alnA[i], alnB[i]):
            case ('-', bmatch):
                #There's a gap in A.
                aProfile += ' '
                alnStr += '|'
                pval = bmatch.profile
                #Do an (inefficient, but it's short) linear search through 
                #the quantiles and find where the current value fits in. 
                for i, val in enumerate(bQuant):
                    if(pval >= val):
                        #Create a unicode block character so that I can 
                        #render the alignment as text. 
                        printCode = chr(0x2581+i)
                bProfile = bProfile + printCode
            case (amatch, '-'):
                #There's a gap in B.
                bProfile += ' '
                alnStr += '|'
                pval = amatch.profile
                for i, val in enumerate(aQuant):
                    if(pval >= val):
                        printCode = chr(0x2581+i)
                aProfile = aProfile + printCode
            case (amatch, bmatch):
                #Both sequences are aligned here.
                pval = bmatch.profile
                for i, val in enumerate(bQuant):
                    if(pval >= val):
                        printCode = chr(0x2581+i)
                bProfile = bProfile + printCode
                
                pval = amatch.profile
                for i, val in enumerate(aQuant):
                    if(pval >= val):
                        printCode = chr(0x2581+i)
                aProfile = aProfile + printCode
                
                dist = delta(amatch, bmatch)
                for i, val in enumerate(alnQuant):
                    if(dist >= val):
                        printCode = chr(0x2581+i)
                alnStr += printCode
    #The first line is the profile of A, then the sequence of A, 
    #then the delta scores, then the sequence of B, and finally B's profile. 
    printFolded([aProfile,
                "".join([stringer(x) for x in alnA]),
                alnStr,
                "".join([stringer(x) for x in alnB]),
                bProfile], foldLen)
        
    #print(aProfile)
    #print("".join([stringer(x) for x in alnA]))
    #print(alnStr)
    #print("".join([stringer(x) for x in alnB]))
    #print(bProfile)


# ## Modified Functions

# In[18]:


def ElemAln_return(alnA, alnB, delta, stringer=lambda x:x, foldLen=80):
    
    alnStr = ""
    aProfile = ""
    bProfile = ""
    alnStats = []
    aStats = []
    bStats = []
    #Build up arrays of the stringer and profiles, so that we can display that information in a normalized format. 
    for i in range(len(alnA)):
        match (alnA[i], alnB[i]):
            case ('-', bmatch):
                bStats.append(bmatch.profile)
            case (amatch, '-'):
                aStats.append(amatch.profile)
            case (amatch, bmatch):
                aStats.append(amatch.profile)
                bStats.append(bmatch.profile)
                alnStats.append(delta(alnA[i], alnB[i]))
                
    if len(alnStats) > 0: 
        alnQuant = np.quantile(np.array(alnStats), [x/8 for x in range(1,9)])
        alnQuant = np.linspace(np.min(alnStats), np.max(alnStats), num=8)
    else:
        alnQuant = []
        
    aQuant = np.quantile(np.array(aStats), [x/8 for x in range(1,9)])
    bQuant = np.quantile(np.array(bStats), [x/8 for x in range(1,9)])
    #I was doing quantiles, and I might want to revisit this for the delta,but especially for profiles it's a lot better-looking to use a linear space.
    
    aQuant = np.linspace(np.min(aStats), np.max(aStats), num=8)
    bQuant = np.linspace(np.min(bStats), np.max(bStats), num=8)
    
    
    for i in range(len(alnA)):
        match (alnA[i], alnB[i]):
            case ('-', bmatch):
                #There's a gap in A.
                aProfile += ' '
                alnStr += '|'
                pval = bmatch.profile
                #Do an (inefficient, but it's short) linear search through the quantiles and find where the current value fits in. 
                for i, val in enumerate(bQuant):
                    if(pval >= val):
                        #Create a unicode block character so that I can render the alignment as text. 
                        printCode = chr(0x2581+i)
                bProfile = bProfile + printCode
            case (amatch, '-'):
                #There's a gap in B.
                bProfile += ' '
                alnStr += '|'
                pval = amatch.profile
                for i, val in enumerate(aQuant):
                    if(pval >= val):
                        printCode = chr(0x2581+i)
                aProfile = aProfile + printCode
            case (amatch, bmatch):
                #Both sequences are aligned here.
                pval = bmatch.profile
                for i, val in enumerate(bQuant):
                    if(pval >= val):
                        printCode = chr(0x2581+i)
                bProfile = bProfile + printCode
                
                pval = amatch.profile
                for i, val in enumerate(aQuant):
                    if(pval >= val):
                        printCode = chr(0x2581+i)
                aProfile = aProfile + printCode
                
                dist = delta(amatch, bmatch)
                for i, val in enumerate(alnQuant):
                    if(dist >= val):
                        printCode = chr(0x2581+i)
                alnStr += printCode
    #The first line is the profile of A, then the sequence of A, then the delta scores, then the sequence of B, and finally B's profile. 

    string_lsts = [aProfile,
                "".join([stringer(x) for x in alnA]),
                alnStr,
                "".join([stringer(x) for x in alnB]),
                bProfile]
    #printFolded(string_lsts, foldLen)
    return string_lsts, alnStats

        
    #print(aProfile)print("".join([stringer(x) for x in alnA]))print(alnStr)print("".join([stringer(x) for x in alnB]))print(bProfile)


# In[19]:


#use this as normalization function when using no normalization
def identity(x):
    return x

#use this when using probability transformation
def probability(profile_array):
    total_sum = np.sum(profile_array)
    if total_sum == 0:
        return np.zeros_like(profile_array)
    return profile_array / total_sum


# In[20]:


def seqProfile_modified(alpha, relative = False, threshold = 0):
    '''
    input: alpha value, relative difference (boolean), threshold (boolean)
    
        Computes the sum of profile and sequence loss for aligned bases 
        Loss calculation depends on the relative and threshold settings
        
    output: total loss 
    '''

    #Absolute difference with no threshold for low contribitions
    if relative is False and threshold == 0:
        def comparator(elem1, elem2):
            deltaSequence = 0 if elem1.sequence == elem2.sequence else alpha
            deltaProfile = (1-alpha) * abs(elem1.profile - elem2.profile)
            return deltaSequence + deltaProfile

    #Absolute difference with threshold
    if relative is False and threshold != 0:
        def comparator(elem1, elem2):
            #set threshold to 0.01 for soxn, 0.005 for svb
            if elem1.profile < threshold and elem2.profile < threshold:
                deltaSequence = 0 if elem1.sequence == elem2.sequence else alpha
                deltaProfile = 0
            #if either of the profiles is high contribution, then calculate as usual
            else: 
                deltaSequence = 0 if elem1.sequence == elem2.sequence else alpha
                deltaProfile = (1-alpha) * abs(elem1.profile - elem2.profile)
            return deltaSequence + deltaProfile

    
    #Relative difference with no threshold for high contribitions
    if relative is True and threshold == 0:         
        def comparator(elem1, elem2):
            deltaSequence = 0 if elem1.sequence == elem2.sequence else alpha
            deltaProfile = (1-alpha) * (abs(elem1.profile - elem2.profile)/max([elem1.profile, elem2.profile]))
            return deltaSequence + deltaProfile

    #Relative difference with threshold 
    if relative is True and threshold != 0:         
        def comparator(elem1, elem2):
            if elem1.profile < threshold and elem2.profile < threshold:
                deltaSequence = 0 if elem1.sequence == elem2.sequence else alpha
                deltaProfile = 0
            #if either of the profiles is high contribution, then calculate as usual
            else: 
                deltaSequence = 0 if elem1.sequence == elem2.sequence else alpha
                deltaProfile = (1-alpha) * (abs(elem1.profile - elem2.profile)/max([elem1.profile, elem2.profile]))          
            return deltaSequence + deltaProfile

    
    return comparator
