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


#Read the PDB file and get the sequence of each chain from the SEQRES records and the coordinates of each aminno acid from the ATOM records. Use the alpha carbon coordinate of amino acids as their coordinates.
def read_pdb_file(pdb_file):
    # Initialize the sequence and coordinates.
    sequence = []
    coordinates = []
    # Open the PDB file.
    with open(pdb_file) as f:
        # For each line in the PDB file,
        for line in f:
            # If the line starts with ATOM,
            if line.startswith('ATOM'):
                # Get the atom name.
                atom_name = line[12:16].strip()
                # If the atom name is CA,
                if atom_name == 'CA':
                    # Get the residue name.
                    residue_name = line[17:20].strip()
                    # Get the residue number.
                    residue_number = int(line[22:26].strip())
                    # Get the x coordinate.
                    x = float(line[30:38].strip())
                    # Get the y coordinate.
                    y = float(line[38:46].strip())
                    # Get the z coordinate.
                    z = float(line[46:54].strip())
                    # Append the residue name and residue number to the sequence.
                    sequence.append([residue_name, residue_number])
                    # Append the x, y, and z coordinates to the coordinates.
                    coordinates.append([x, y, z])
            # If the line starts with SEQRES,
            elif line.startswith('SEQRES'):
                # Get the chain ID.
                chain_id = line[11]
                # Get the residues.
                residues = line[19:].strip().split()
                # For each residue,
                for residue in residues:
                    # Append the chain ID and residue to the sequence.
                    sequence.append([chain_id, residue])
    # Return the sequence and coordinates.
    return sequence, coordinates

# This function will compute the pairwise distances between all amino acids in a chain.
def compute_pairwise_distances_from_sequence(sequence):
    # Initialize the pairwise distances.
    pairwise_distances = []
    # For each amino acid in the chain, compute the pairwise distances.
    for i in range(len(sequence)):
        # Initialize the list of distances.
        distances = []
        # For each amino acid in the chain, compute the distance.
        for j in range(len(sequence)):
            # Compute the distance.
            distance = euclidean_distance(sequence[i], sequence[j])
            # Add the distance to the list.
            distances.append(distance)
        # Add the list of distances to the list of pairwise distances.
        pairwise_distances.append(distances)
    # Return the pairwise distances.
    return pairwise_distances


# This function will compute the Euclidean distance between two amino acids.
def euclidean_distance(amino_acid_1, amino_acid_2):
    # Get the x, y, and z coordinates of the first amino acid.
    x1, y1, z1 = amino_acid_1[2:]
    # Get the x, y, and z coordinates of the second amino acid.
    x2, y2, z2 = amino_acid_2[2:]
    # Compute the distance.
    distance = math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
    # Return the distance.
    return distance

# This function will compute the pairwise distances between all amino acids in a chain.
def compute_pairwise_distances_from_coordinates(coordinates):
    # Initialize the pairwise distances.
    pairwise_distances = []
    # For each amino acid in the chain, compute the pairwise distances.
    for i in range(len(coordinates)):
        # Initialize the list of distances.
        distances = []
        # For each amino acid in the chain, compute the distance.
        for j in range(len(coordinates)):
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
def output_amino_acid_pairs(pairwise_distances, s, d):
    # Initialize the amino acid pairs.
    amino_acid_pairs = []
    # For each amino acid in the chain,
    for i in range(len(pairwise_distances)):
        # For each amino acid in the chain,
        for j in range(len(pairwise_distances)):
            # If the distance is less than or equal to D,
            if pairwise_distances[i][j] <= d:
                # If the amino acids are at least S amino acids apart in the sequence,
                if i - j >= s:
                    # Append the amino acid pair to the amino acid pairs.
                    amino_acid_pairs.append([i, j])
    # Return the amino acid pairs.
    return amino_acid_pairs
# main function
def main():
    # Read the user input.
    pdb_file, s, d = read_user_input()
    # Read the PDB file and get the sequence and coordinates.
    sequence, coordinates = read_pdb_file(pdb_file)
    # Compute the pairwise distances between all amino acids in the sequence.
    pairwise_distances_from_sequence = compute_pairwise_distances_from_sequence(sequence)
    # Compute the pairwise distances between all amino acids in the coordinates.
    pairwise_distances_from_coordinates = compute_pairwise_distances_from_coordinates(coordinates)
    # Output the amino acid pairs whose distance in 3-dimensions is less than or equal to D, and which are at least S amino acids apart in the sequence.
    amino_acid_pairs = output_amino_acid_pairs(pairwise_distances_from_sequence, s, d)
    # Print the amino acid pairs.
    print(amino_acid_pairs)

# Call the main function.
main()
