# 2D-3D_1-stable_manifold
Extract 1-stable manifold from arbitrary simplicial complex

1. Folder Content

<Main Folder>
|
|-- densityRips.cpp (main: Compute 1-stable manifold. input: BIN file; output: vert.txt edge.txt)
|-- Simplicial2Complex.h (contains classes of simplicial complex and discrete vector field)
|-- persistence.h (persistence class)
|-- Makecommand (sample command for building densityRips)
|-- gridComplex.txt (sample input for 2D case)
|
|-- <Matlab Helper> (list only the most frequently used ones)
|   ..
|   ..
|   |-- Draw1stable.m (Helper that draws the vert.txt and edge.txt file)
|   |-- MorsePost.m (Postprocessing steps)
|   |-- LoadAllen.m (Load particular data file and generate BIN file for densityRips)
<end>


2. Compile
g++ densityRips.cpp -I./phat/include -std=c++11 -o density3D


3. Running command
./density3D Allen.txt output_vert.txt output_edge.txt 20
