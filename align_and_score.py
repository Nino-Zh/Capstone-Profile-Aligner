import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pysam
import random
import pyBigWig
import enum
import os
import argparse
import sys

from aligner import gotoh
from aligner import seqProfile_modified
from aligner import printElemAln
from aligner import ElemAln_return


class Element:
    def __init__(self, sequence, profile):
        self.sequence = sequence
        self.profile = profile
    def __str__(self):
        return "{0:s} {1:f}".format(self.sequence, self.profile)
        
def getSequences_synthetic(normalize, genomes, genomeFnames, contributionsFnames):
    '''
    input: profile loss normalization function
    
        Retrieves sequences from FASTA files
        and TF binding predictions from BW files
        
    output: nucleotide sequence of interest and corresponding TF binding predictions
    
    '''
    elements = dict()
    region_coords = []  
    for index, organism in enumerate(genomes):
        # Open the FASTA file
        orgGenome = pysam.FastaFile(genomeFnames[organism])
        
        # Define the chromosome name (assuming a single chromosome)
        chrom = list(orgGenome.references)[0] 
        #print(organism, chrom)
        
        # Define the start and end coordinates for the entire sequence
        orgStart = 0
        orgEnd = orgGenome.get_reference_length(chrom)
        orgSequence = orgGenome.fetch(chrom, orgStart, orgEnd)

        
        # Fetch the sequence for the entire chromosome
        orgSequence = orgGenome.fetch(chrom, orgStart, orgEnd)
        
        # Open the BigWig file
        bw = pyBigWig.open(contributionsFnames[organism])
        
        
        # Fetch the profile values for the entire chromosome
        orgProfile = bw.values(chrom, orgStart, orgEnd)
        #orgProfile = np.nan_to_num(orgProfile, nan=0.0)
        bw.close()
        
        # Normalize the profile
        profileAr = normalize(np.array(orgProfile))
        
        # Store the sequence and profile in the elements dictionary
        elements[organism] = []
        for i in range(len(orgSequence)):
            elements[organism].append(Element(orgSequence[i].upper(), profileAr[i]))
        
        # Store the region coordinates
        region_coords.append(f'{organism}: {chrom} {orgStart}-{orgEnd}')
    
    return elements, region_coords# Format vcf as the alignment


# Format vcf as the alignment
def vcf_to_alignment(vcf_file, reference_seq_path):
    
    raw_seq = ""
    with open(reference_seq_path, 'r') as file:
        for line in file:
            # Skip the header line (starts with >)
            if not line.startswith('>'):
                # Add sequence lines to string, removing whitespace
                raw_seq += line.strip()

    vcf_data = []
    with open(vcf_file, 'r') as f:
        for line in f:
            # Skip headers
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            vcf_data.append({
                'idx': int(fields[1]) - 1,  # Mutation position
                'ref': fields[3],  # Reference sequence
                'alt': fields[4],  # Mutated sequence
                'info': fields[7]  # Indel mutation info
            })

  
    
    alnA = []
    alnB = []

    raw_seq_pos = 0 # index along the raw sequence
    for entry in vcf_data:

        vcf_pos = entry['idx']

        # Copy everything to alnA and alnB until reaching a mutation
        while raw_seq_pos < len(raw_seq) and raw_seq_pos < vcf_pos:
            alnA.append(raw_seq[raw_seq_pos])
            alnB.append(raw_seq[raw_seq_pos])
            raw_seq_pos += 1

        # Insertion mutations
        if 'SVTYPE=INS' in entry['info']:
            # add the first base in both
            alnA.append(raw_seq[raw_seq_pos])
            
            # Add gaps in reference (alnA) and inserted bases in alternative (alnB)
            ins_len = len(entry['alt']) - len(entry['ref'])
            alnA.extend(['-']*ins_len)
            alnB.extend(entry['alt'])

            raw_seq_pos += 1
        
        elif 'SVTYPE=DEL' in entry['info']:
            # add the first base in both
            alnB.append(raw_seq[raw_seq_pos])
            
            del_len = len(entry['ref']) - len(entry['alt'])
            alnA.extend(entry['ref'])
            alnB.extend(['-']*del_len)

            raw_seq_pos += del_len
        # Substitution mutations
        else:
            alnA.append(entry['ref'])
            alnB.append(entry['alt'])
            raw_seq_pos += 1


    return alnA, alnB






# Get mutation data from the alignment 

def get_mutation_data(alnA, alnB):
    '''
    Input: alnA and alnB (list of strings)
    
    Outputs the mutation information for each base in A and B. 
    Iterates over each base and checks if it has been substituted or deleted in the alternative sequence. 
    Info column contains information on base substitution. 

    Output: two tables with three columns (base index, mutation type, info)
    '''
    
    mutation_data_seqA = []
    mutation_data_seqB = []

    # Alignment outputs (contains gaps)
    # Convert all in caps to ensure mismatch is not due to case difference 
    alnA = [str(element).upper() for element in alnA]  
    alnB = [str(element).upper() for element in alnB]
    
    gap = '-'

    # Input sequences
    seqA =  [element for element in alnA if element != gap]
    seqB =  [element for element in alnB if element != gap]

    current_base_idx_A = 0

    # TABLE 1
    # Mutation data of the reference sequence A 
    for element_idx in range(len(alnA)): 
        # Iterate through each aligned pairs

        # If we have a base (rather than an indel)
        if alnA[element_idx] != gap:
    
            # Check for substitution/no mutations
            # When the aligned pair is a base
            if alnB[element_idx] != gap:
                
                # Check if the aligned element matches
                if alnA[element_idx] == alnB[element_idx]:
                    mutation_data_seqA.append([current_base_idx_A, 'No Mutation', 'N/A'])
                    
                # If there is a mismatch between the aligned pairs --> substitution
                else: 
                    # Store the mutation info in the 3rd column
                    mutation_data_seqA.append([current_base_idx_A, 'Substitution', alnB[element_idx]])



            # If aligned pair is a gap (in seqB), it's a deletion mutation:
            else: 
                mutation_data_seqA.append([current_base_idx_A, 'Deletion', 'N/A'])

            current_base_idx_A += 1

        # If we have a gap in seqA, then it's an insertion and we skip


    current_base_idx_B = 0 
    
    # TABLE 2
    # Mutation data of the alternative sequence B
    for element_idx in range(len(alnB)): 
        # Iterate through each aligned pairs

        # If we have a base (rather than an indel)
        if alnB[element_idx] != gap:
            
    
            # Check for substitution/no mutations
            # When the aligned pair is a base
            if alnA[element_idx] != gap:
                
                # Check if the aligned element matches
                if alnB[element_idx] == alnA[element_idx]:
                    mutation_data_seqB.append([current_base_idx_B, 'No Mutation', 'N/A'])
                    
                # If there is a mismatch between the aligned pairs --> substitution
                else: 
                    # Store the mutation info in the 3rd column
                    mutation_data_seqB.append([current_base_idx_B, 'Substitution', alnA[element_idx]])

    

            # If aligned pair is a gap (in seqA), it's a deletion mutation:
            if alnA[element_idx] == gap:
                mutation_data_seqB.append([current_base_idx_B, 'Insertion', 'N/A'])

            current_base_idx_B += 1
    

    return mutation_data_seqA, mutation_data_seqB

            
# Compare table 1s and table 2s of true and algorithm mutation data to score the alignment

def score_mutation_tables(algo_table, true_table, is_table2=False):
    """
    Compare mutation tables (Table 1 or Table 2) and calculate the score.
    
    Args:
        algo_table: Algorithm-generated mutation data
        true_table: Ground-truth mutation data
        is_table2: If True, only compare insertions (ignores substitutions/no mutations)
    
    Returns:
        Score (int)
    """
    # Convert tables to dictionaries: {base_idx: (mutation_type, info)}
    algo_dict = {entry[0]: (entry[1], entry[2]) for entry in algo_table}
    true_dict = {entry[0]: (entry[1], entry[2]) for entry in true_table}
    
    # All indices (union of algo and truth indices)
    all_indices = set(algo_dict.keys()).union(true_dict.keys())
    score = 0
    
    for idx in all_indices:
        algo_entry = algo_dict.get(idx, None)
        true_entry = true_dict.get(idx, None)
        
        # Skip if comparing Table 2 and mutation is not an insertion
        if is_table2:
            if algo_entry and algo_entry[0] != "Insertion":
                algo_entry = None  # Treat as non-existent for scoring
            if true_entry and true_entry[0] != "Insertion":
                true_entry = None
        
        # Case 1: Both entries exist
        if algo_entry and true_entry:
            if algo_entry == true_entry:
                score += 1  # Correct match
            else:
                score -= 1  # Mismatch
        # Case 2: Missing entry in one of the tables
        elif algo_entry != true_entry:
            score -= 1  # Penalize missing/mismatched entries
    
    return score

def calculate_total_score(algo_table1, algo_table2, true_table1, true_table2):
    """Compute total score by comparing Table 1 and Table 2."""
    # Score Table 1 (substitutions/deletions)
    table1_score = score_mutation_tables(algo_table1, true_table1)
    #print(table1_score)
    # Score Table 2 (insertions only)
    table2_score = score_mutation_tables(algo_table2, true_table2, is_table2=True)
    #print(table2_score)
    return table1_score + table2_score

