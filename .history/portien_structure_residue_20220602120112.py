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

# this function computer the euclidean distance between two amino acids
def compute_euclidean_distance(coordinates1, coordinates2):
    # Get the x, y, and z coordinates of the first amino acid.
    x1, y1, z1 = coordinates1
    # Get the x, y, and z coordinates of the second amino acid.
    x2, y2, z2 = coordinates2
    # Compute the euclidean distance.
    distance = math.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
    # Return the euclidean distance.
    return distance

# main function
def main():
    # Get the PDB file.
    pdb_file = sys.argv[1]
    # Get the D.
    D = float(sys.argv[2])
    # Get the S.
    S = int(sys.argv[3])
    # Read the PDB file and get the sequence and coordinates.
    sequence, coordinates = read_pdb_file(pdb_file)
    # Compute the pairwise distances.
    pairwise_distances = compute_pairwise_distances(sequence, coordinates, D, S)
    # Print the pairwise distances.
    print(pairwise_distances)

# Call the main function.
if __name__ == '__main__':
    main()

