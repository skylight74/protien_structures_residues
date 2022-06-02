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


# This function reads the PDB file and returns the sequence of each chain and the coordinates of each
# amino acid.
def read_pdb(pdb_file):
    # Open the PDB file.
    pdb_file = open(pdb_file, 'r')
    # Initialize the sequence and the coordinates of each amino acid.
    sequence = []
    coordinates = []
    # Read the PDB file line by line.
    for line in pdb_file:
        # If the line starts with ATOM, get the amino acid sequence and the coordinates.
        if line[0:4] == 'ATOM':
            # Get the amino acid sequence.
            sequence.append(line[13:15])
            # Get the coordinates.
            coordinates.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    # Close the PDB file.
    pdb_file.close()
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
    s = int(sys.argv[2])
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

# This function will plot the pairwise distances between all amino acids in a chain.
def plot_pairwise_distances(sequence, pairwise_distances):
    # Initialize the figure.
    fig = plt.figure()
    # Initialize the subplot.
    ax = fig.add_subplot(111)
    # Initialize the list of colors.
    colors = ['red', 'green', 'blue', 'yellow', 'orange', 'purple', 'brown', 'pink', 'gray', 'black']
    # For each amino acid in the chain, plot the pairwise distances.
    for i in range(len(sequence)):
        # Initialize the list of distances.
        distances = []
        # For each amino acid in the chain, compute the distance.
        for j in range(len(sequence)):
            # Compute the distance.
            distance = pairwise_distances[i][j]
            # Add the distance to the list.
            distances.append(distance)
        # Plot the distances.
        ax.plot(range(len(sequence)), distances, color=colors[i])
    # Set the x-axis label.
    ax.set_xlabel('Amino acid')
    # Set the y-axis label.
    ax.set_ylabel('Distance')
    # Set the title.
    ax.set_title('Pairwise distances')
    # Show the plot.
    plt.show()
    # Save the plot.
    fig.savefig('pairwise_distances.png')
    # Close the figure.
    plt.close(fig)
    # Return the figure.
    return fig

# This function will plot the pairwise distances between all amino acids in a chain.
def plot_pairwise_distances_2(sequence, pairwise_distances):
    # Initialize the figure.
    fig = plt.figure()
    # Initialize the subplot.
    ax = fig.add_subplot(111)
    # Initialize the list of colors.
    colors = ['red', 'green', 'blue', 'yellow', 'orange', 'purple', 'brown', 'pink', 'gray', 'black']
    # For each amino acid in the chain, plot the pairwise distances.
    for i in range(len(sequence)):
        # Initialize the list of distances.
        distances = []
        # For each amino acid in the chain, compute the distance.
        for j in range(len(sequence)):
            # Compute the distance.
            distance = pairwise_distances[i][j]
            # Add the distance to the list.
            distances.append(distance)
        # Plot the distances.
        ax.plot(range(len(sequence)), distances, color=colors[i])
    # Set the x-axis label.
    ax.set_xlabel('Amino acid')
    # Set the y-axis label.
    ax.set_ylabel('Distance')
    # Set the title.
    ax.set_title('Pairwise distances')
    # Show the plot.
    plt.show()
    # Save the plot.
    fig.savefig('pairwise_distances_2.png')
    # Close the figure.
    plt.close(fig)
    # Return the figure.
    return fig

# main
if __name__ == '__main__':
    # the main function will run the program.
    # Read the user input from the command line.
    pdb_file, s, d = read_user_input()
    # Read the PDB file.
    sequence, coordinates = read_pdb(pdb_file)
    # Compute the pairwise distances.
    pairwise_distances = compute_pairwise_distances(sequence, coordinates)
    # Output the amino acid pairs whose distance in 3-dimensions is less than or equal to D, and which are at least S amino acids apart in the sequence.
    output_results(sequence, pairwise_distances, s, d)
    # Plot the pairwise distances between all amino acids in a chain.
    plot_pairwise_distances(sequence, pairwise_distances)
    # Plot the pairwise distances between all amino acids in a chain.
    plot_pairwise_distances_2(sequence, pairwise_distances)
    # Exit the program.
    sys.exit()

