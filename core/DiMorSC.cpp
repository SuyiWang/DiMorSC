/*
The software DiMoSC computes 1-stable manifold in arbitrary dimension from function defined on 2-dimensional simplicial complex.


Dependency:
	PHAT, C++11, openMP

	
Compile:
	// g++ DiMoSC.cpp -O3 -I./phat/include -openMP -std=c++11 -o DiMoSC -w
	// Debug switch - slower but output more information

	
Execute command:
	./DiMorSC <input_file> <output_prefix> <persistence_threshold> <DIM> [saved_persis_pair]

	
Parameters:
	// argv[1] - BIN file, format specified below.
	// argv[2] - Prefix of output files
	// argv[3] - Persistence threshold for simplification
	// argv[4] - Dimension of points (2 or 3)
	// argv[5] - If specified, the program will load previously 
				 computed persistence pairing.
	// argv[6] - reserved.

	
Input specification:
	The input file is composed of six blocks with the follow meaning:
	<1. num_vertices>
	<2. vertex list>
	<3. num_edges>
	<4. edge list>
	<5. num_triangles>
	<6. triangle list>
	
	1: contains a single int32 specifying the total number of vertices
	2: contains geographic information of vertices in blocks.
	   Each of the <num_vertices> blocks has the following information:
	   [v.x1 v.x2 ... v.xn v.f], xn bounded by DIM.
	   Currently, supports DIM = 2 and 3. For higher dimensions, please adjust the constant MAX_DIM in simplicial2complex.h
	   v.x_i and v.f are all of type double
	3: contains a single int32 specifying the total number of edges
	4: specifies how edges are connected using index pairs.
	   A pair [v1 v2] specifies an edge connecting vertices of 
	   index v1 and v2.
	5: contains a single int32 specifying the total number of triangles
	6: specifies how triangles are connected using index tuples.
	   A tuple [v1 v2 v3] specifies an triangles formed by connecting vertices of index v1, v2 and v3.

The output are two ascii files <output_prefix>_vert.txt, [output_prefix]_edge.txt.
See example for more input and output details.

Additional Assumption: 
	Number of simplices would not be more than MAX_INT / 3.
	If it exceed the maximum allowed numnber, please use split and merge algorithm.

*/

#define DEBUG 0

int DIM;			// Dimension of data - 3 for 3D, 2 for 2D

// For higher dimension, change MAX_DIM
#define MAX_DIM 3			// - Will be used in Simplex.h
#define EPS_compare 1e-8	// - used in comparison functions

#include <ctime>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <unordered_map>
using namespace std;

#include "Simplex.h"
#include "persistence.h"
#include "DiscreteVField.h"
#include "Simplicial2Complex.h"


int main(int argc, char* argv[]){
	// Output filenames
	string output_file[2];
	// pre-save filename
	string pre_save;
	bool use_pre_save = false;
	
	
	// Threshold for simplification
	double et_delta = 0;
	double ve_delta = 0;
	
	
	//  Resolving parameters
    if (argc < 5){
    	// must provide 4 parameters.
    	// argv[1] - inut_file
    	// argv[2] - output_file
    	// argv[3] - persistence_threshold
    	// argv[4] - dimension
    	// argv[5] - use_previous - optional
    	// argv[6] - triangle threshold - under experiment
		cout << "Usage: ./DiMorSC <input_file> <output_file> <persistence_threshold> <dimension> [use_previous]"
		<<endl;
		return 0;
    }else{
		output_file[0] = string(argv[2]) + "_vert.txt";
		output_file[1] = string(argv[2]) + "_edge.txt";
		ve_delta = atof(argv[3]);
		DIM = atoi(argv[4]);
	}
    if (argc >= 6){
		pre_save = string(argv[5]);
		use_pre_save = true;
    }
    if (argc >= 7){
    	et_delta = atof(argv[6]);
    }else{
    	et_delta = ve_delta;
    }
	
	cout << argc-1 << " parameters detected"<< endl;
	
	
	//  Output memory cost if in debug mode.
	if (DEBUG){
		cout << "debug mode\n";
		// memory usage info.
		cout << sizeof(Simplex) << " " << sizeof(Vertex) << " " << sizeof(Edge) << " "
			 << sizeof(Triangle) << " " << sizeof(Triangle*)<< endl;
	}
	
	
	//  Main pipeline
	Simplicial2Complex K;
	clock_t startTime;
    clock_t testTime;
    double time_passed = 0;
	if (!use_pre_save){
		
		//  Loading input file
		cout << "Reading in simplicial complex...\n";
		startTime = clock();
		K.buildComplexFromFile2_BIN(argv[1]);
		testTime = clock(); time_passed = (startTime - testTime) / (double) CLOCKS_PER_SEC;
		cout << "Done in " << time_passed << " \n";
		cout.flush();
		// cin.get(); // has test
	

		// Build psudo morse function
		cout << "Building pseudo-Morse function...\n";
		startTime = clock();
		K.buildPsuedoMorseFunction();
		testTime = clock(); time_passed = (startTime - testTime) / (double) CLOCKS_PER_SEC;
		cout << "Done in " << time_passed << " \n";
		cout.flush();
		// cin.get(); 

		
		//  Build filtration
		cout << "Building filtration...\n";
		startTime = clock();
		K.buildFiltrationWithLowerStar();
		testTime = clock(); time_passed = (startTime - testTime) / (double) CLOCKS_PER_SEC;
                cout << "Done in " << time_passed << " \n";
		cout.flush();
		// cin.get(); // has test

		
		//  Computing persistence pairs using PHAT
		cout << "Computing persistence pairs...\n";
		startTime = clock();
		K.PhatPersistence();
		testTime = clock(); time_passed = (startTime - testTime) / (double) CLOCKS_PER_SEC;
                cout << "Done in " << time_passed << " \n";
		cout.flush();
		
		
		//  Writing persistence info.
		cout << "Writing pre_saved_data...\n";
		startTime = clock();
		K.write_presave(argv[2]);
		testTime = clock(); time_passed = (startTime - testTime) / (double) CLOCKS_PER_SEC;
                cout << "Done in " << time_passed << " \n";
		cout.flush();
	}else{
		cout << "Reading in pre_saved_data...\n";
		startTime = clock();
		K.Load_Presaved(argv[1], pre_save);
		testTime = clock(); time_passed = (startTime - testTime) / (double) CLOCKS_PER_SEC;
                cout << "Done in " << time_passed << " \n";
		cout.flush();
		// cin.get();
		
		// Build psudo morse function
		cout << "Building pseudo-Morse function...\n";
                startTime = clock();
		K.buildPsuedoMorseFunction();
		testTime = clock(); time_passed = (startTime - testTime) / (double) CLOCKS_PER_SEC;
                cout << "Done in " << time_passed << " \n";
		cout.flush();
	}

	
	//  Cancelling persistence pairs
	cout << "Cancelling persistence pairs with delta " << ve_delta << "\n";
	startTime = clock();
	// Cancellation does not use function values on simplicies
	K.cancelPersistencePairs(ve_delta);
	testTime = clock(); time_passed = (startTime - testTime) / (double) CLOCKS_PER_SEC;
        cout << "Done in " << time_passed << " \n";
	cout.flush();
	
	
	//  Writing output
	startTime = clock();
	K.outputArcs(output_file[0], output_file[1], et_delta);
	testTime = clock(); time_passed = (startTime - testTime) / (double) CLOCKS_PER_SEC;
        cout << "Results written in " << time_passed << " \n";
	return 0;
}
