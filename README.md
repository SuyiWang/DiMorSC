## DiMoSC - Discrete Morse on Simplicial Complex
A Toolbox for skeletonization

The software DiMorSC extracts near-linear structures from simplicial complex. The near-linear structure is often referred to as skeleton or 1D non-linear structure (e.g. trees and graphs that can be decomposed into lines). The simplicial complex is a general representation of domains potentially with extra data. Most common data structures such as image stacks and point cloud can be conveniently converted to simplicial complex and it can further represent more irregular and complicated domains.


## Overview
The tool box is a commandline based software. It has one core module - DiMorSC and two extension modules - Triangulate and graph2tree. The core module DiMorSC extracts skeleton from simplicial complex. The module Triangulate converts point cloud on regular grid (i.e. all points have integer coordinates) to simplicial complex, which can be read by DiMorSC. The module graph2tree converts the output of DiMorSC to trees.

Helper functions to visualize objects such as volume, graph and simplicial complex and to process image data are provided in python.


## Prerequisites
1. Core code DiMorSC depends on PHAT to compute persistence pairs. PHAT is available on github here (https://github.com/blazs/phat)  or bitbucket (https://bitbucket.org/phat-code/phat) (tends to have newer version)- download and put PHAT file in "extern" folder.
2. tree simplification and triangulation requires boost (https://www.boost.org/) - download and put boost in 'extern' folder or build and install it.
3. Python helper is for Python3.5+ and it depends on Numpy, sci-kit and Vispy. A pip installer will be provided for simpler installation.

## Compile DiMorSC
In general simply execute "make" should compile the code. 

## Running DiMorSC
./bin/DiMorSC \<input_file> \<output_prefix> \<persistence_threshold> \<dimension> [use_previous]
./bin/Triangulate \<density_file\> \<fill\> \<2 (2D)/3 (3D)\>
./bin/graph2tree \<graphfile.ini\>

## Test data
example for running input in data folder

./bin/DiMorSC data/OP_7_trunc.bin output/OP_7 5 3

## More Test files
Test files are provided for checking if the software is successfully setup. 
https://drive.google.com/open?id=0B_mktdY-pu10MHRlaUZneUhJMEE

## Input and output description
## DiMorSC input 

input format (binary):
[\<int32\> * 1] [\<double\> * 4 * n] [\<int32\> * 1] [\<int\> * 2 * m] [\<int32\>*1] [\<int32\> * 3 * k]

  * How to interpret it: There are altogether 3 blocks in the binary file, specifying vertex, edge and triangle information, respectively. The first block contains vertex information. It starts with a 32-bit integer n and is followed by 4*n double (64-bit float) - every 4 double describe the coordinate and function value of a vertex (for 3D case). For 2D case, the this block should contain 3*n doubles. The second block contains starts with 32-bit integer m and is followed by 2*m 32-bit integers. The integers are the vertex indices, which start from zero. Similar to the second block, the third block contains triangle information that starts with a 32-bit integer k and is followed by 3*k 32-bit integers specifying the indices of triangle vertices.
  
## DiMorSC output

output format: (3D)
_vert.txt:
Each line specifies a vertex: x y z f c. xyz are the coordinates of the vertex. f denotes the function value on the vertex. c decalres the criticality of the vertex - (-1) means ordinary.
_edge.txt:
Each line specifies an edge: v1 v2 critial density. v1 v2 are two indices of the edge - index starts from 1. (This is DIFFERENT! It is designed this way because the post processing scripts are all in MATLAB). critical decalres the criticality of the edge - (-1) means ordinary. density specify the maximum intensity of its supporting saddle - later, the edges with low density can be removed using a threshold.

## Triangulate input
input format (binary):
[\<int32\> * 1] [\<double\> * 4 * n]
Similar to DiMorSC input, it starts with a 32-bit integer n and is followed by 4*n double (64-bit float) - every 4 double describe the coordinate and function value of a point (for 3D case).


## Terminology
The skeleton is mathematically modelled as 1-stable manifold from Discrete Morse Theory. See reference for more mathematical details. 

## Reference
https://doi.org/10.1101/321489

This work is credited to Computer Science and Engineering Department, The Ohio State University and 
 Cold Spring Harbor Laboratory


## Contact
suyi.wang@ini.usc.edu


