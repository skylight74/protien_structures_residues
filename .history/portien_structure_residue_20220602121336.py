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

# The PDB file is a text file that contains the ATOM and SEQRES records. The ATOM record
# contains the coordinates of the atoms that make up the structure. The SEQRES record contains
# the sequence of the amino acids in the chain.


# The following function reads the PDB file and returns a dictionary of the coordinates of the
# atoms in the PDB file. The keys of the dictionary are the atom names. The values of the
# dictionary are the coordinates of the atoms.
def read_pdb(pdb_file):
    atom_coordinates = {}
    with open(pdb_file, "r") as f:
        for line in f:
            if line[0:4] == "ATOM":
                atom_name = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atom_coordinates[atom_name] = [x, y, z]
    print("The coordinates of the atoms in the PDB file are:")
    print(atom_coordinates)
    return atom_coordinates

# The following function reads the SEQRES record and returns the sequence of the amino acids in the
# chain.
def read_seqres(seqres_record):
    chain_sequence = seqres_record[11:].strip()
    return chain_sequence

# the following function computes the Euclidean distance between two points in 3-dimensions.
def compute_distance(point1, point2):
    distance = math.sqrt((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2 + (point1[2] - point2[2])**2)
    return distance

# the following function computes the pairwise distances between all amino acids in a chain for all chains in the structure.
def compute_pairwise_distances(atom_coordinates, chain_sequence):
    pairwise_distances = {}
    for chain_id in atom_coordinates:
        for i in range(len(chain_sequence)):
            for j in range(i+1, len(chain_sequence)):
                atom_name1 = chain_sequence[i] + chain_id
                atom_name2 = chain_sequence[j] + chain_id
                distance = compute_distance(atom_coordinates[chain_id][i], atom_coordinates[chain_id][j])
                pairwise_distances[(atom_name1, atom_name2)] = distance
    return pairwise_distances

# the following function returns the amino acid pairs whose distance in 3-dimensions is less than or equal to D, and which are at least S amino acids apart in the sequence.
def find_close_pairs(pairwise_distances, D, S):
    close_pairs = []
    for pair in pairwise_distances:
        if pairwise_distances[pair] <= D:
            close_pairs.append(pair)
    close_pairs_filtered = []
    for pair in close_pairs:
        if pair[0][1] - pair[1][1] >= S:
            close_pairs_filtered.append(pair)
    return close_pairs_filtered

# main function
def main():
    # get the name of the PDB file from the user
    pdb_file = input("Enter the name of the PDB file: ")
    # read the PDB file
    atom_coordinates = read_pdb(pdb_file)
    # read the SEQRES record
    with open(pdb_file, "r") as f:
        for line in f:
            if line[0:6] == "SEQRES":
                seqres_record = line
    chain_sequence = read_seqres(seqres_record)
    # compute the pairwise distances between all amino acids in a chain for all chains in the structure
    pairwise_distances = compute_pairwise_distances(atom_coordinates, chain_sequence)
    # find the amino acid pairs whose distance in 3-dimensions is less than or equal to D, and which are at least S amino acids apart in the sequence
    close_pairs = find_close_pairs(pairwise_distances, D, S)
    # output the results to the user
    print("The amino acid pairs whose distance in 3-dimensions is less than or equal to D, and which are at least S amino acids apart in the sequence are:")
    print(close_pairs)
