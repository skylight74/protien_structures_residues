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


# The following function reads the PDB file and returns a dictionary of the coordinates of the atoms in the PDB file for each chain. The dictionary is keyed by the chain ID.
def read_pdb(pdb_file):
    # Open the PDB file and read the lines
    f = open(pdb_file, 'r')
    lines = f.readlines()
    f.close()
    # Initialize the dictionary to store the coordinates in
    coordinates = {}
    # Loop over the lines and read the coordinates
    for line in lines:
        # Skip the line if it is not an ATOM line
        if line[0:4] != 'ATOM':
            continue
        # Get the chain ID
        chain_ID = line[21:22]
        # Get the residue number
        residue_number = line[22:26]
        # Get the atom type
        atom_type = line[13:14]
        # Get the coordinates
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        # Add the coordinates to the dictionary
        if chain_ID not in coordinates:
            coordinates[chain_ID] = {}
        coordinates[chain_ID][residue_number] = {
            'type': atom_type,
            'coordinates': np.array([x, y, z])
        }
    # Return the coordinates
    return coordinates


# The following function reads the SEQRES record and returns a dictionary of the sequence of the chain in the SEQRES record. The dictionary is keyed by the chain ID.
def read_seqres(seqres_file):
    # Open the SEQRES file and read the lines
    f = open(seqres_file, 'r')
    lines = f.readlines()
    f.close()
    # Initialize the dictionary to store the sequence in
    sequence = {}
    # Loop over the lines and read the sequence
    for line in lines:
        # Skip the line if it is not a SEQRES line
        if line[0:6] != 'SEQRES':
            continue
        # Get the chain ID
        chain_ID = line[11:12]
        # Get the sequence
        sequence[chain_ID] = line[19:].strip()
    # Return the sequence
    return sequence


# the following function computes the Euclidean distance between two points in 3-dimensions
def euclidean_distance(point1, point2):
    # Compute the Euclidean distance between the two points
    distance = math.sqrt((point1[0] - point2[0])**2 +
                         (point1[1] - point2[1])**2 +
                         (point1[2] - point2[2])**2)
    # Return the distance
    return distance


# The following function computes the pairwise distances between all amino acids in a chain. The function takes as input the coordinates of the amino acids in the chain and returns a dictionary of the pairwise distances. The dictionary is keyed by the chain ID.
def compute_pairwise_distances(coordinates):
    # Initialize the dictionary to store the pairwise distances in
    pairwise_distances = {}
    # Loop over the chain IDs
    for chain_ID in coordinates:
        # Initialize the dictionary to store the pairwise distances for this chain in
        pairwise_distances[chain_ID] = {}
        # Loop over the residue numbers
        for residue_number in coordinates[chain_ID]:
            # Get the coordinates of the current amino acid
            current_coordinates = coordinates[chain_ID][residue_number][
                'coordinates']
            # Initialize the dictionary to store the pairwise distances for this amino acid in
            pairwise_distances[chain_ID][residue_number] = {}
            # Loop over the residue numbers
            for other_residue_number in coordinates[chain_ID]:
                # Get the coordinates of the other amino acid
                other_coordinates = coordinates[chain_ID][
                    other_residue_number]['coordinates']
                # Compute the Euclidean distance between the two amino acids
                distance = euclidean_distance(current_coordinates,
                                              other_coordinates)
                # Add the distance to the dictionary
                pairwise_distances[chain_ID][residue_number][
                    other_residue_number] = distance
    # Return the pairwise distances
    return pairwise_distances


# the following function will report the amino acid pairs whose distance in 3-dimensions is
#less than or equal to D, and which are at least S amino acids apart in the sequence (i.e., if one
#amino acid is the ith amino acid, the other has to be i+Sth or more, or vice versa).
def report_pairs(pairwise_distances, D, S):
    # Initialize the list to store the pairs in
    pairs = []
    # Loop over the chain IDs
    for chain_ID in pairwise_distances:
        # Loop over the residue numbers
        for residue_number in pairwise_distances[chain_ID]:
            # Loop over the residue numbers
            for other_residue_number in pairwise_distances[chain_ID]:
                # Get the pairwise distance
                distance = pairwise_distances[chain_ID][residue_number][
                    other_residue_number]
                # Check if the distance is less than or equal to D
                if distance <= D:
                    # Check if the residues are at least S amino acids apart
                    if abs(int(residue_number) - int(other_residue_number)) >= S:
                        # Add the pair to the list
                        pairs.append((chain_ID, residue_number,
                                      other_residue_number))
    # Return the pairs
    return pairs

# the following function will output the protien pairs using FASTA format
def output_pairs(pairwise_distances, coordinates, sequence):
    # Initialize the list to store the pairs in
    pairs = []
    # Loop over the chain IDs
    for chain_ID in pairwise_distances:
        # Loop over the residue numbers
        for residue_number in pairwise_distances[chain_ID]:
            # Loop over the residue numbers
            for other_residue_number in pairwise_distances[chain_ID]:
                # Get the pairwise distance
                distance = pairwise_distances[chain_ID][residue_number][
                    other_residue_number]
                # Check if the distance is less than or equal to D
                if distance <= D:
                    # Check if the residues are at least S amino acids apart
                    if abs(int(residue_number) - int(other_residue_number)) >= S:
                        # Add the pair to the list
                        pairs.append((chain_ID, residue_number,
                                      other_residue_number))
    # Return the pairs
    return pairs

# the following function will compute the RMSD between two chains.
def compute_rmsd(coordinates1, coordinates2):
    # Initialize the list to store the pairwise distances in
    rmsd = []
    # Loop over the chain IDs
    for chain_ID in coordinates1:
        # Loop over the residue numbers
        for residue_number in coordinates1[chain_ID]:
            # Get the coordinates of the current amino acid
            current_coordinates = coordinates1[chain_ID][residue_number][
                'coordinates']
            # Get the coordinates of the other amino acid
            other_coordinates = coordinates2[chain_ID][residue_number][
                'coordinates']
            # Compute the Euclidean distance between the two amino acids
            distance = euclidean_distance(current_coordinates,
                                          other_coordinates)
            # Add the distance to the list
            rmsd.append(distance)
    # Return the RMSD
    return math.sqrt(sum(rmsd) / len(rmsd))

# the following fuction will use computer RMSD results to mesure the similarity of two chains.
def compare_chains(coordinates1, coordinates2):
    # Compute the pairwise distances between the two chains
    pairwise_distances = compute_pairwise_distances(coordinates1)
    # Report the pairs whose distance is less than or equal to D
    pairs = report_pairs(pairwise_distances, D, S)
    # Compute the RMSD between the two chains
    rmsd = compute_rmsd(coordinates1, coordinates2)
    # Return the RMSD and the pairs
    return rmsd, pairs


def main():
    # Get the command line arguments
    pdb_file = sys.argv[1]
    seqres_file = sys.argv[1]
    D = float(sys.argv[2])
    S = float(sys.argv[3])
    # Read the PDB file
    coordinates = read_pdb(pdb_file)
    # Read the SEQRES file
    sequence = read_seqres(seqres_file)
    # Compute the pairwise distances between the amino acids in the PDB file
    pairwise_distances = compute_pairwise_distances(coordinates)
    # Report the pairs whose distance is less than or equal to D
    pairs = report_pairs(pairwise_distances, D, S)
    # Output the pairs in the format required by the Rosetta script
    output_pairs(pairs, 'pairs.txt')

main()