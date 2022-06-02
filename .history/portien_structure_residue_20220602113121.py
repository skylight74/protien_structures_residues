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
    # Initialize the sequence of each chain.
    sequence = []
    # Initialize the coordinates of each amino acid.
    coordinates = []
    # Open the PDB file.
    with open(pdb_file, 'r') as pdb_file:
        # For each line in the PDB file,
        for line in pdb_file:
            # If the line starts with the ATOM record,
            if line[0:4] == 'ATOM':
                # If the line contains the CA atom,
                if line[12:16] == 'CA  ':
                    # Get the chain ID.
                    chain_id = line[21:22]
                    # Get the amino acid.
                    amino_acid = line[17:20]
                    # Get the x coordinate.
                    x = float(line[30:38])
                    # Get the y coordinate.
                    y = float(line[38:46])
                    # Get the z coordinate.
                    z = float(line[46:54])
                    # Add the amino acid to the sequence of the chain.
                    sequence[chain_id].append(amino_acid)
                    # Add the x, y, and z coordinates to the coordinates of the amino acid.
                    coordinates[chain_id].append([x, y, z])
            # If the line starts with the SEQRES record,
            if line[0:6] == 'SEQRES':
                # Get the chain ID.
                chain_id = line[11:12]
                # Get the amino acid sequence.
                amino_acid_sequence = line[19:70]
                # Add the amino acid sequence to the sequence of the chain.
                sequence[chain_id] = amino_acid_sequence.split()
                # Initialize the coordinates of the amino acid.
                coordinates[chain_id] = []
    # Return the sequence and the coordinates.
    return sequence, coordinates
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
    sequence, coordinates = read_pdb_file(pdb_file)
    # Compute the pairwise distances between all amino acids in a chain.
    pairwise_distances = compute_pairwise_distances(sequence, coordinates)
    # Output amino acid pairs whose distance in 3-dimensions is less than or equal to D, and which are at least S amino acids apart in the sequence.
    output_results(sequence, pairwise_distances, s, d)
    # Return the sequence and the coordinates.
    return sequence, coordinates

# Call the main function.
main()
