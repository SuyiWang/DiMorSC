# 2D-3D_1-stable_manifold
This is a software that extract 1-stable manifold from arbitrary simplicial complex. The input is a simplicial complex with function value defined on vertices. The output is a graph that represents the 1-stable manifold.

# Folder Content Highlights

## Main Folder

1. densityRips.cpp - Main program. Computes 1-stable manifold from function defined on simplicial complex. input: BIN file; output: vert.txt edge.txt. There is only 1 parameter (persistence threshold) to adjust.

2. Makecommand - sample command for building densityRips.

## Matlab Helper (list only the most frequently used ones)
1. LoadAllen.m/LoadPartha.m - Load particular data file and generate BIN file for densityRips. There is one adjustable parameter - cut off density threshold. In LoadAllen.m it is at line 93. In Load Partha, it is at Line 92. The threshold removes a data point when its density is lower than the threshold, therefore reduce the size of simplicial complex.

2. MorsePost.m - This file converts the output of densityRips, which is usually a graph, to a tree using maximum spanning tree. It also smoothes the branch a little and removes low persistence branches. Then write tree data to .swc file. There is one adjustable parameter - branch persistence threshold.

3. Draw1stable.m - This file visualizes the output of densityRips, without any post processing.

# Prerequisites
1. main code DensityRips depends on PHAT to compute persistence pairs. PHAT is available on github here (https://github.com/blazs/phat) - put PHAT file in the project root folder.
2. Matlab code depends on matlab_bgl (https://www.mathworks.com/matlabcentral/fileexchange/10922-matlabbgl) - put matlab_bgl in matlab_helper folder
3. vaa3d module depends on 'vaa3d_matlab_io', which is a folder in vaa3d source file - this is provided in case it might not be available online.

# Creating input file from raw data
Helper files "LoadPartha" and "LoadAllen" in <matlabHelper> folder are provided to generate input from raw data. 

# Compilation
## densityRips
g++ densityRips.cpp -I./phat/include -std=c++11 -o density3D
## Triangulation
g++ Triangulation.cpp -std=c++11 -o triangulate

# Running densityRips

format:

./density3D [inputfile] [output_file_vertex] [output_file_edge] [persistence_threshold]

example:

./density3D Allen.bin output_vert.txt output_edge.txt 20

# Test files
Test files are provided for checking if the software is successfully setup. 
https://drive.google.com/open?id=0B_mktdY-pu10MHRlaUZneUhJMEE

If the output on your machine is the same as the ones given in the output folder, congratulations, the setup is complete and you are good to go! (Test data good in Ubuntu only)


# Post processing
Helper file "Morse_post" post process the "vert.txt" "edge.txt" by 1) keeping the maximum spanning tree of the largest connected component of the input; 2) assigning function value based on the density on each remaining vertex and their distance to the root; 3) simplify branch using persistence.


# Contact
wang.3066@buckeyemail.osu.edu

