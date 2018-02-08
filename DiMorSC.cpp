/*
DiMoSC computes 1-stable manifold (max-sadlle path) from function 
defined on simplicial complex.

Input: Simplicial complex - binary version.
Output: [file_name]_vert.txt, [file_name]_edge.txt - ascii version.

See example for more input and output details.

Assumption: Number of simplices would not be more than MAX_INT / 3.
If it exceed the maximum allowed numnber, please use split and merge algorithm.

Usage: ./densityRips <input_file> <output_file> <persistence_threshold> <dimension> [use_previous]
*/

// g++ DiMoSC.cpp -I./phat/include -std=c++11 -o DiMoSC -w
// Debug switch - slower but output more information
#define DEBUG 0


int DIM;			// Dimension of data - 3 for 3D, 2 for 2D

// For higher dimension, change MAX_DIM
#define MAX_DIM 3			// - Will be used in Simplex.h
#define EPS_compare 1e-8	// - used in comparison functions

#include <ctime>
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
	string output_file[2];
	string pre_save;
	bool use_pre_save = false;
	
	double et_delta = 0;
	double ve_delta = 0;
	
    if (argc < 5){
    	// must provide 4 parameters.
    	// argv[1] - inut_file
    	// argv[2] - output_file
    	// argv[3] - persistence_threshold
    	// argv[4] - dimension
    	// argv[5] - use_previous - optional
    	// argv[6] - triangle threshold - under experiment
		cout << "Usage: ./DiMOSC <input_file> <output_file> <persistence_threshold> <dimension> [use_previous]"
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
	if (DEBUG){
		cout << "debug mode\n";
		// memory usage info.
		cout << sizeof(Simplex) << " " << sizeof(Vertex) << " " << sizeof(Edge) << " "
			 << sizeof(Triangle) << " " << sizeof(Triangle*)<< endl;
	}
	
	Simplicial2Complex K;
	if (!use_pre_save){
		cout << "Reading in simplicial complex...\n";
		K.buildComplexFromFile2_BIN(argv[1]);
		cout << "Done\n";
		cout.flush();
		// cin.get(); // has test
	

		// Build psudo morse function
		cout << "Building pseudo-Morse function...\n";
		K.buildPsuedoMorseFunction();
		cout << "Done\n";
		cout.flush();
		// cin.get(); 


		cout << "Building filtration...\n";
		K.buildFiltrationWithLowerStar();
		cout << "Done\n";
		cout.flush();
		// cin.get(); // has test


		cout << "Computing persistence pairs...\n";
		K.PhatPersistence();
		cout << "Done!\n";
		cout.flush();
		
		
		cout << "Writing pre_saved_data...\n";
		K.write_presave(argv[2]);
		cout << "Done!\n";
		cout.flush();
	}else{
		cout << "Reading in pre_saved_data...\n";
		K.Load_Presaved(argv[1], pre_save);
		cout << "Done\n";
		cout.flush();
		// cin.get();
		
		// Build psudo morse function
		cout << "Building pseudo-Morse function...\n";
		K.buildPsuedoMorseFunction();
		cout << "Done\n";
		cout.flush();
	}

	cout << "Cancelling persistence pairs with delta " << ve_delta << "\n";
	// Cancellation does not use function values on simplicies
	K.cancelPersistencePairs(ve_delta);
	cout << "Done\n";
	cout.flush();

	K.outputArcs(output_file[0], output_file[1], et_delta);
	return 0;
}
