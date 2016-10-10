# 2D-3D_1-stable_manifold
This is a software that extract 1-stable manifold from arbitrary simplicial complex. The input is a simplicial complex with function value defined on vertices. The output is a graph that represents the 1-stable manifold.

# Folder Content

## Main Folder

1. densityRips.cpp           (main: Compute 1-stable manifold. input: BIN file; output: vert.txt edge.txt)
2. Simplicial2Complex.h      (contains classes of simplicial complex and discrete vector field)
3. persistence.h             (persistence class)
4. Makecommand               (sample command for building densityRips)
5. gridComplex.txt (sample input for 2D case)

## Matlab Helper (list only the most frequently used ones)
1. Draw1stable.m (Helper that draws the vert.txt and edge.txt file)
2. MorsePost.m (Postprocessing steps)
3. LoadAllen.m (Load particular data file and generate BIN file for densityRips)



# Compile
g++ densityRips.cpp -I./phat/include -std=c++11 -o density3D


# Running command
./density3D Allen.txt output_vert.txt output_edge.txt 20
