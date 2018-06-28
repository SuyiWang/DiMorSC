## DiMoSC - Discrete Morse on Simplicial Complex
A Toolbox for skeletonization

The software DiMorSC extracts near-linear structures from simplicial complex. The near-linear structure is often referred to as skeleton or 1D non-linear structure (e.g. trees and graphs that can be decomposed into lines). The simplicial complex is a general representation of domains. Most common data structure such as image stacks and point cloud can be conveniently converted to simplicial complex and it can further represent more irregular and complicated domains.


## Overview
The tool box is a commandline based software. It has one core module - DiMorSC and two extension modules - Triangulate and graph2tree. The core module DiMorSC extracts skeleton from simplicial complex. The module Triangulate converts point cloud on regular grid (i.e. all points have integer coordinates) to simplicial complex. The module graph2tree converts the output of DiMorSC to trees.


## Prerequisites
1. Core code DiMorSC depends on PHAT to compute persistence pairs. PHAT is available on github here (https://github.com/blazs/phat)  or bitbucket (https://bitbucket.org/phat-code/phat) (tends to have newer version)- download and put PHAT file in "core" folder.
2. Matlab code depends on matlab_bgl (https://www.mathworks.com/matlabcentral/fileexchange/10922-matlabbgl) - put matlab_bgl in tree_simplification folder

## Compile DiMorSC
In general simply execute "make" should compile the code (g++ 11 is required). If not, the code can be compiled using the following command.
## densityRips
g++ DiMoSC.cpp -I./phat/include -std=c++11 -o DiMoSC
## Triangulation
g++ Triangulation.cpp -std=c++11 -o triangulate

## Running DiMorSC
format:

./DiMorSC \<input_file> \<output_prefix> \<persistence_threshold> \<dimension> [use_previous]

## On test data
example for running input in data folder

./bin/DiMorSC data/OP_7_trunc.bin output/OP_7 5 3

## More Test files
Test files are provided for checking if the software is successfully setup. 
https://drive.google.com/open?id=0B_mktdY-pu10MHRlaUZneUhJMEE

If the output on your machine is the same as the ones given in the output folder, congratulations, the setup is complete and you are good to go! (Test data good in Ubuntu only)

## Input and output
## input 

input format (binary):
[\<int32\> * 1] [\<double\> * 4 * n] [\<int32\> * 1] [\<int\> * 2 * m] [\<int32\>*1] [\<int32\> * 3 * k]

  * How to interpret it: There are altogether 6 blocks in the binary file. The first block contains a single integer called n represented with 32-bit. The second block contains double (64-bit) with total number of 4n - every 4 double describe the location and function value of a vertex (for 3D case). For 2D case, the this block should contain 3n of doubles. The third block contains a single 32-bit integer called m. The fourth block contains 2m of 32-bit integers specifying the two vertex of an edge, assuming the indices of vertices start from zero. The fifth block contains a single 32-bit integer k. The sixth block contains 3k of 32-bit integers specifying the three vertices of a triangle.
  
## output

output format: (3D)
vert.txt:
Each line specifies a vertex: x y z f critical. xyz are the coordinates of the vertex. f denotes the function value on the vertex. critical decalres the criticality of the vertex - (-1) means ordinary.
edge.txt:
Each line specifies an edge: v1 v2 critial density. v1 v2 are two indices of the edge - index starts from 1. (This is DIFFERENT! It is designed this way because the post processing scripts are all in MATLAB). critical decalres the criticality of the edge - (-1) means ordinary. density specify the maximum intensity of its supporting saddle - later, the edges with low density can be removed using a threshold.

## Terminology
The skeleton is mathematically modelled as 1-stable manifold from Discrete Morse Theory. See reference for more mathematical details. 

## Reference
https://doi.org/10.1101/321489

This work is credited to Computer Science and Engineering Department, The Ohio State University and 
 Cold Spring Harbor Laboratory


## Contact
suyi.wang@ini.usc.edu


