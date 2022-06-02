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
import math
from tabulate import tabulate
# the following fuction takes a protein structure as input and returns the secondary structure of each chain
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
    # open the pdb file
    pdb_file = open(pdb_file, 'r')

    index = 0
    # read the pdb file line by line
    for line in pdb_file:
        # if the line is the SEQRES record, get the amino acid sequence
        if line[0:4] == 'ATOM' and line[13:16] == 'CA ':
            # for each chain id, get the amino acid id and then add a coordinate to the dictionary
            chain_id = line[21]
            amino_acid_residue_name = line[17:20]
            amino_acid_coordinate = [
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54])
            ]
            if chain_id in seq:
                seq[chain_id].append(
                    [amino_acid_residue_name, amino_acid_coordinate])
            else:
                seq[chain_id] = [[
                    amino_acid_residue_name, amino_acid_coordinate
                ]]
    # close the pdb file
    pdb_file.close()
    # print the amino acid sequences and the coordinates
    # print(seq['A'])
    # print(coords)
    # return the amino acid sequences and the coordinates
    return seq


# get protien coordinates in all chains using protien name and sequnce as input
def get_protein_coordinates_in_seq(protein_sequence, chain_id):
    protien_coords = []
    for i in protein_sequence[chain_id]:
        protien_coords.append(i[1])
    return protien_coords


# the following fuction calculates the euclidian distance between two points in 3-dimensions
def distance(point1, point2):
    return math.sqrt((point1[0] - point2[0]) ** 2 + (point1[1] - point2[1]) ** 2 + (point1[2] - point2[2]) ** 2)




# the follwoing fuction takes a chain id and a list of amino acid coordinates as input and returns the distance between all amino acids in the chain
def get_distance_between_all_amino_acids(protein_coordinates, chain_id, D,S):

    distance_list = []
    for i in range(len(protein_coordinates[chain_id])):
        for j in range(i+1, len(protein_coordinates[chain_id])):
            eucladian_distance=distance(protein_coordinates[chain_id][i][1], protein_coordinates[chain_id][j][1])
            abs_distance=abs(i-j)
            if eucladian_distance <= D and abs_distance >= S:
                distance_list.append([protein_coordinates[chain_id][i][0], protein_coordinates[chain_id][j][0], eucladian_distance, abs_distance])
    return distance_list


# the following fuction computes the pairwise distance between all amino acids in a chain whose distance is less than or equal to D and which are at least S amino acids apart in the sequence
def get_pairwise_distance(protein_sequence, S, D):
    # initialize the pairwise distance as a dictionary with chain id as the key and amino acid pairwise distance as the value
    pairwise_distance = {}
    # initialize the pairwise distance as a dictionary with chain id as the key and amino acid pairwise distance as the value
    # return the pairwise distance
    return pairwise_distance
# run the program
def main():
    # check if we are in the main function
    # get the input arguments
    pdb_file = sys.argv[1]
    D = float(sys.argv[2])
    S = int(sys.argv[3])
    output_file = "output.txt"
    # print("reading the pdb file")
    # read the pdb file
    seq = read_pdb(pdb_file)
    match_list = {}
    for chain_id in seq:
        match_list[chain_id]=get_distance_between_all_amino_acids(seq,chain_id,D,S)
    # write the output to the output file
    # tabluate the pairwise distance
    table = tabulate(match_list['A'], headers=['Amino Acid 1', 'Amino Acid 2', 'Distance', 'Absolute Distance'])
    print(table)
    with open("output_file.out", 'w') as f:
        # f.write(str(match_list))
        f.write("Chain_id\tAmino_acid_1\tAmino_acid_2\tDistance\tAbsolute_distance\n")
        for chain_id in match_list:
            for i in match_list[chain_id]:
                f.write(chain_id+"\t"+str(i[0])+"\t"+str(i[1])+"\t"+str(i[2])+"\t"+str(i[3])+"\n")
    print("done")

# run the main function
if __name__ == '__main__':
    main()
    print('Done')
