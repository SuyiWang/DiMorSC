# 2D-3D_1-stable_manifold
This is a software that extract 1-stable manifold from arbitrary simplicial complex. The input is a simplicial complex with function value defined on vertices. The output is a graph that represents the 1-stable manifold.

# Folder Content

## Main Folder

1. densityRips.cpp           (main: Compute 1-stable manifold. input: BIN file; output: vert.txt edge.txt)
2. Simplicial2Complex.h      (contains classes of simplicial complex and discrete vector field)
3. persistence.h             (persistence class)
4. Makecommand               (sample command for building densityRips)
5. gridComplex.txt           (sample input for 2D case)

## Matlab Helper (list only the most frequently used ones)
1. Draw1stable.m (Helper that draws the vert.txt and edge.txt file)
2. MorsePost.m (Postprocessing steps)
3. LoadAllen.m (Load particular data file and generate BIN file for densityRips)

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


# Post processing
Helper file "Morse_post" post process the "vert.txt" "edge.txt" by 1) keeping the maximum spanning tree of the largest connected component of the input; 2) assigning function value based on the density on each remaining vertex and their distance to the root; 3) simplify branch using persistence.


# Contact
wang.3066@buckeyemail.osu.edu

