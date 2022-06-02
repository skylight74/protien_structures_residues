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

# This function will read the PDB file and get the sequence of each chain from the SEQRES records and the coordinates of each aminno acid from the ATOM records.
def read_pdb_file(pdb_file):
    # Initialize the lists of amino acid sequences and coordinates.
    amino_acid_sequences = []
    amino_acid_coordinates = []
    # Initialize the list of chains.
    chains = []
    # Initialize the list of amino acid names.
    amino_acid_names = ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"]
    # Initialize the list of amino acid coordinates.
    amino_acid_coordinates = [[0, 0, 0]]
    # Initialize the chain identifier.
    chain_identifier = 1
    # Open the PDB file.
    with open(pdb_file) as f:
        # Read the lines in the PDB file.
        lines = f.readlines()
    # For each line in the PDB file.
    for line in lines:
        # If the line begins with "ATOM", get the amino acid name.
        if line[:4] == "ATOM":
            # Get the amino acid name.
            amino_acid_name = line[17:20]
            # If the amino acid name is not in the list of amino acid names, exit.
            if amino_acid_name not in amino_acid_names:
                sys.exit("The amino acid name is not in the list of amino acid names.")
            # If the amino acid name is in the list of amino acid names, get the chain identifier.
            chain_identifier = line[21]
        # If the line begins with "SEQRES", get the amino acid sequence.
        elif line[:6] == "SEQRES":
            # Get the amino acid sequence.
            amino_acid_sequence = line[11:].rstrip()
            # If the amino acid sequence contains an X, exit.
            if "X" in amino_acid_sequence:
                sys.exit("The amino acid sequence contains an X.")
            # Add the amino
            amino_acid_sequences.append(amino_acid_sequence)
        # If the line begins with "ATOM" and the amino acid name is in the list of amino acid names, get the coordinates.
        elif line[:4] == "ATOM" and amino_acid_name in amino_acid_names:
            # Get the coordinates.
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            # Add the coordinates to the list of coordinates.
            amino_acid_coordinates.append([x, y, z])
    # For each amino acid sequence.
    for amino_acid_sequence in amino_acid_sequences:
        # If the amino acid sequence is not the same length as the list of coordinates, exit.
        if len(amino_acid_sequence) != len(amino_acid_coordinates):
            sys.exit("The amino acid sequence is not the same length as the list of coordinates.")
        # If the amino acid sequence is the same length as the list of coordinates, add the amino acid sequence to the list of chains.
        else:
            chains.append(amino_acid_sequence)
    # Return the list of chains.
    return chains




# This function computes the Euclidean distance between two points.
def euclidean_distance(point1, point2):
    # Compute the distance.
    distance = math.sqrt((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2 + (point1[2] - point2[2])**2)
    # Return the distance.
    return distance

# This function computes the pairwise distances between all amino acids in a chain.
def compute_pairwise_distances(sequence, coordinates):
    # Initialize the pairwise distances.
    pairwise_distances = []
    # For each amino acid in the chain, compute the pairwise distances.
    for i in range(len(sequence)):
        # Initialize the list of distances.
        distances = []
        # For each amino acid in the chain, compute the distance.
        for j in range(len(sequence)):
            # Compute the distance.
            distance = euclidean_distance(coordinates[i], coordinates[j])
            # Add the distance to the list.
            distances.append(distance)
        # Add the list of distances to the list of pairwise distances.
        pairwise_distances.append(distances)
    # Return the pairwise distances.
    return pairwise_distances


# This function will read the user input from the command line.
def read_user_input():
    # Get the PDB file.
    pdb_file = sys.argv[1]
    # Get the S value.
    s = float(sys.argv[2])
    # Get the D value.
    d = float(sys.argv[3])
    # Return the PDB file, S value, and D value.
    return pdb_file, s, d

# This function will . Output amino acid pairs whose distance in 3-dimensions is less than or equal to D, and which are at least S amino acids apart in the sequence.
def output_results(sequence, pairwise_distances, s, d):
    # For each chain, compute the amino acid pairs whose distance in 3-dimensions is less than or equal to D, and which are at least S amino acids apart in the sequence.
    for i in range(len(sequence)):
        # Initialize the list of amino acid pairs.
        amino_acid_pairs = []
        # For each amino acid in the chain, compute the amino acid pairs whose distance in 3-dimensions is less than or equal to D, and which are at least S amino acids apart in the sequence.
        for j in range(len(sequence)):
            # If the distance is less than or equal to D, and the distance is at least S amino acids apart in the sequence, add the amino acid pair to the list.
            if pairwise_distances[i][j] <= d and j - i >= s:
                amino_acid_pairs.append([sequence[i], sequence[j]])
        # Print the list of amino acid pairs.
        print(amino_acid_pairs)
# main function
def main():
    # Read the user input.
    pdb_file, s, d = read_user_input()
    # Read the PDB file and get the sequence and the coordinates of each amino acid.
    sequence, coordinates = read_pdb(pdb_file)
    # Compute the pairwise distances between all amino acids in a chain.
    pairwise_distances = compute_pairwise_distances(sequence, coordinates)
    # Output amino acid pairs whose distance in 3-dimensions is less than or equal to D, and which are at least S amino acids apart in the sequence.
    output_results(sequence, pairwise_distances, s, d)
    # Return the sequence and the coordinates.
    return sequence, coordinates

# Call the main function.
main()
