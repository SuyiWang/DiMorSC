# 2D-3D_1-stable_manifold - DiMoSC
This is software (Discrete Morse on Simplicial Complex) that extract 1-stable manifold from a function defined on arbitrary simplicial complex. The input is a binary file describing the function on simplicial complex. The output is a graph representing the 1-stable manifold.

# Folder Content Highlights

## Main Folder

1. DiMoSC.cpp - Main program.
* input format (binary):
[\<int32\> * 1] [\<double\> * 4 * n] [\<int32\> * 1] [\<int\> * 2 * m] [\<int32\>*1] [\<int32\> * 3 * k]

  * How to interpret it: There are altogether 6 blocks in the binary file. The first block contains a single integer called n represented with 32-bit. The second block contains double (64-bit) with total number of 4n - every 4 double describe the location and function value of a vertex (for 3D case). For 2D case, the this block should contain 3n of doubles. The third block contains a single 32-bit integer called m. The fourth block contains 2m of 32-bit integers specifying the two vertex of an edge, assuming the indices of vertices start from zero. The fifth block contains a single 32-bit integer k. The sixth block contains 3k of 32-bit integers specifying the three vertices of a triangle.

  * output format: (3D)
vert.txt:
Each line specifies a vertex: x y z f critical. xyz are the coordinates of the vertex. f denotes the function value on the vertex. critical decalres the criticality of the vertex - (-1) means ordinary.
edge.txt:
Each line specifies an edge: v1 v2 critial density. v1 v2 are two indices of the edge - index starts from 1. (This is DIFFERENT! It is designed this way because the post processing scripts are all in MATLAB). critical decalres the criticality of the edge - (-1) means ordinary. density specify the maximum intensity of its supporting saddle - later, the edges with low density can be removed using a threshold.

2. Makecommand - sample command for building DiMoSC.

## Matlab Helper (list only the most frequently used ones)
1. LoadAllen.m/LoadPartha.m - Load particular data file and generate BIN file for densityRips.

2. MorsePost.m - This file converts the output of densityRips, which is usually a graph, to a tree using maximum spanning tree. It also smoothes the branch a little and removes low persistence branches. Then write tree data to .swc file.

3. Draw1stable.m - This file visualizes the output of densityRips, without any post processing.

# Prerequisites
1. main code DensityRips depends on PHAT to compute persistence pairs. PHAT is available on github here (https://github.com/blazs/phat)  or bitbucket (https://bitbucket.org/phat-code/phat) (tends to have newer version)- put PHAT file in the project root folder.
2. Matlab code depends on matlab_bgl (https://www.mathworks.com/matlabcentral/fileexchange/10922-matlabbgl) - put matlab_bgl in matlab_helper folder
3. vaa3d module depends on 'vaa3d_matlab_io', which is a folder in vaa3d source file - this is provided in case it might not be available online.

# Creating input file from raw data
Helper files "LoadPartha" and "LoadAllen" in \<matlabHelper\> folder are provided to generate input from raw data. 

# Compilation
## densityRips
g++ DiMoSC.cpp -I./phat/include -std=c++11 -o DiMoSC
## Triangulation
g++ Triangulation.cpp -std=c++11 -o triangulate

# Running densityRips

format:

./DiMoSC \<input_file> \<output_file> \<persistence_threshold> \<dimension> [use_previous]

example:

./DiMoSC Allen.bin output 20 3

# Test files
Test files are provided for checking if the software is successfully setup. 
https://drive.google.com/open?id=0B_mktdY-pu10MHRlaUZneUhJMEE

If the output on your machine is the same as the ones given in the output folder, congratulations, the setup is complete and you are good to go! (Test data good in Ubuntu only)


# Post processing
Helper file "Morse_post" post process the "vert.txt" "edge.txt" by 1) keeping the maximum spanning tree of the largest connected component of the input; 2) assigning function value based on the density on each remaining vertex and their distance to the root; 3) simplify branch using persistence.


# Contact
wang.3066@buckeyemail.osu.edu

