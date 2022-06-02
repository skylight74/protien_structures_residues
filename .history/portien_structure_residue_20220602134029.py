# In this assignment, your goal is to implement a method to find amino acid pairs that are in close
# proximity in the structure (i.e., in three dimensions) but far away in the sequence. In particular,
# your program will get input a protein structure (which may be composed of multiple amino acid
# chains) as a single PDB file and two numbers, one integer, S, to be used as the distance threshold
# in the sequence and a real number, D, to be used as a distance threshold in the stucture. You willl
# then report amino-acid pairs (for each chain of the structure) whose distance in 3-dimensions is
# less than or equal to D, and which are at least S amino acids apart in the sequence (i.e., if one
# amino acid is the ith amino acid, the other has to be i+Sth or more, or vice versa).
# Here is a step by step description of the algorithm you are expected to implement:
# 1) Get a protein structure as input from the user as a PDB (Protein Data Bank) file.
# 2) Read the PDB file and get the sequence of each chain from the SEQRES records and the
# coordinates of each aminno acid from the ATOM records. Use the alpha carbon
# coordinate of amino acids as their coordinates.
# 3) For each chain, compute pairwise distances, as the the Euclidean distance, between all
# amino acid pairs in that chain. Output amino acid pairs whose distance in 3-dimensions is
#less than or equal to D, and which are at least S amino acids apart in the sequence.
# 4) Repeat the process for all chains in the structure.
# 5) Output the results to the user.

from cmath import pi
import sys
import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec

# However, you will only need to read the ATOM and SEQRES records of PDB files. The ATOM
# record contains the coordinates of the atoms that make up the structure. For each amino acid, you
# are only going to use the CA atom (alpha-Carbon) coordinates. The atom records look like below:
# ATOM      1  N   ASP A  30      31.904  -0.904  -0.904  1.00  0.00           N
# ATOM      2  CA  ASP A  30      31.904  -0.904  -0.904  1.00  0.00           C
# ATOM      3  C   ASP A  30      32.904  -0.904  -0.904  1.00  0.00           C
# ATOM      4  O   ASP A  30      33.904  -0.904  -0.904  1.00  0.00           O

# The SEQRES record contains the sequence of the amino acids in the chain. For each chain, you
# are only going to use the first 20 amino acids. The SEQRES record looks like below:
# SEQRES   1 A   20

# the following fuction reads the PDB file get the sequence of each chain from the SEQRES records and the coordinates of each amino acid as the ATOM records and return the amino acid sequences and the coordinates as a dictionary
def read_pdb(pdb_file):
    # initialize sequences as a dictionary with chain id as the key and amino acid sequence list as the value
    seq = {}
    # initialize coordinates as a dictionary with amino acid id as the key and coordinates as the value
    coords = {}
    # open the pdb file
    pdb_file = open(pdb_file, 'r')

    # read the pdb file line by line
    for line in pdb_file:
        # if the line is the SEQRES record, get the amino acid sequence
        if line[0:6] == 'SEQRES':
            # get the chain id
            chain_id = line[11]
            chain_seq =line[19:].rstrip()
            #split the amino acid sequence into a list
            line_seq = chain_seq.split()
            # check if the chain id is in the dictionary of sequences
            if chain_id not in seq:
                # if not, initialize the sequence as an empty string
                seq[chain_id] = ''
            # get the amino acid sequence
                seq[chain_id] += line_seq
            else :
                # if the chain id is in the dictionary of sequences, append the amino acid sequence
                seq[chain_id] += line_seq

        if line[0:4] == 'ATOM':
            # check if the atom is an alpha carbon
            if line[12:16] == 'CA  ':
                # get the amino acid id
                aa_id = line[22:27]
                # get the coordinates
                coords[aa_id] = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
    # close the pdb file
    pdb_file.close()

    print('The amino acid sequences are:')
    for chain_id in seq:
        print(chain_id, seq[chain_id])
    print('The amino acid coordinates are:')
    for aa_id in coords:
        print(aa_id, coords[aa_id])
    # return the amino acid sequences and the coordinates as a dictionary
    return seq, coords
# the following fuction computes the distance between two amino acids
def distance(a, b):
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)


# Compute pairwise distances between all amino acids in a chain. The distance is the Euclidean and return the amino acid pairs whose distance in 3d is less than or equal to D and which are at least S amino acids apart in the sequence.
def compute_distances(seq, coords, D, S):
    # initialize the distance matrix
    dist = np.zeros((len(seq), len(seq)))

    # compute the distance matrix
    for i in range(len(seq)):
        for j in range(i+1, len(seq)):
            dist[i, j] = distance(coords[i], coords[j])
            dist[j, i] = dist[i, j]

    # find the amino acid pairs whose distance in 3d is less than or equal to D and which are at least S amino acids apart in the sequence
    pairs = []
    for i in range(len(seq)):
        for j in range(i+S, len(seq)):
            if dist[i, j] <= D:
                pairs.append([seq[i], seq[j]])

    return pairs

# the following fuction prints the amino acid pairs into the stdout and the output file
def print_pairs(pairs, output_file):
    # print the amino acid pairs into the stdout
    for pair in pairs:
        print(pair[0] + ' ' + pair[1])

    # print the amino acid pairs into the output file
    output_file = open(output_file, 'w')
    for pair in pairs:
        output_file.write(pair[0] + ' ' + pair[1] + '\n')
    output_file.close()

# run the program
def main():
    # check if we are in the main function
    # get the input arguments
    pdb_file = sys.argv[1]
    D = float(sys.argv[2])
    S = int(sys.argv[3])
    output_file = "output.txt"

    # read the pdb file
    seq, coords = read_pdb(pdb_file)

    # compute the amino acid pairs
    pairs = compute_distances(seq, coords, D, S)

    # print the amino acid pairs into the stdout and the output file
    print_pairs(pairs, output_file)

# run the main function
if __name__ == '__main__':
    main()
    print('Done')
