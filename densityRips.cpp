#define DEBUG 1

#include "persistence.h"
//#include "bitmap_image.hpp"
#include <ctime>
#include <fstream>
#include <iostream>
//#define PERSISTENCE 0.75
//#define THRESHOLD
#define REBUILD_COMPLEX false
#define REBUILD_COMPLEX_DISCRETE false
#define RECONSTRUCT_IMAGES false
#define BUILD_SUBSET_MAP false
#define MAKE_SINGLE_SUBSETS false

using namespace std;

int main(int argc, char* argv[]){
//
// Computing 1-stable manifold from arbitrary simplicial complex
// Input: Simplicial complex - vertex list / edge list / triangle list
// Output: A graph - vertex list / edge list
//
    /*
    if (argc == 1){
        argv[1] = "testComplex_ms.txt";
        argv[2] = "outvert2d.txt";
        argv[3] = "outedge2d.txt";
        argv[4] = "24";
    }
    */
	Simplicial2Complex K;

	cout << "Reading in simplicial complex...\n";
	// 2D defined by edge. 3D defined by vertex.
	K.buildComplexFromFile2(argv[1]); //<<<<<<<<<<< Change this if necessary [File: Triangle defined by edge][File2: by vertex]
	//K.outputComplex("testcomplex0.txt");
	cout << "Done\n";

	K.flipAndTranslateVertexFunction();


	cout << "Building pseudo-Morse function...\n";
	K.buildPsuedoMorseFunction();
	cout << "Done\n";

	PersistencePairs P(&K);
	cout << "Building filtration...\n";
	P.buildFiltration();
	cout << "Done\n";
	
	if (DEBUG){
		ofstream filtration("filtration.txt", ios_base::out | ios_base::trunc);
		filtration << setprecision(16);
		for (vector<Simplex*>::iterator it = P.filtration.begin(); it != P.filtration.end(); it++){
			filtration << (*it)->funcValue << " ";
			if ((*it)->dim == 0){
				filtration << "Vertex " << ((Vertex*)(*it))->getVPosition() << endl;
			}
			else if ((*it)->dim == 1){
				filtration << "+ " << ((Edge*)(*it))->getSymPerturb() << "e " <<"Edge " << ((Edge*)(*it))->getEPosition() << endl;
			}
			else if ((*it)->dim == 2){
				filtration <<"+ " << get<0>(((Triangle*)(*it))->getSymPerturb()) << "e + " << get<1>(((Triangle*)(*it))->getSymPerturb()) << "e2 " << "Triangle " << ((Triangle*)(*it))->getTPosition() << endl;
			}
		}
		filtration.close();
	}



	cout << "Computing persistence pairs...\n";
	P.computePersistencePairsWithClear();
	cout << "Done!\n";
	
	if (DEBUG){
		cout << "Outputing persistence pairs...\n";
		P.outputPersistencePairs("persistencePairs.txt");
		cout << "Done\n";
	}

	cout << "Cancelling persistence pairs...\n";
	P.cancelPersistencePairs(atof(argv[4]));
	cout << "Done\n";

	K.outputArcs(argv[2], argv[3]);
}
