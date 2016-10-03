#define DEBUG 1

int DIM;
double delta;
#include "persistence.h"
#include <ctime>
#include <fstream>
#include <iostream>


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
    if (argc == 1){
        argv[1] = "testComplex_ms.txt";
        argv[2] = "outvert2d.txt";
        argv[3] = "outedge2d.txt";
        argv[4] = "24";
    }
    if (argc <=5){
    	DIM = 3;
    }else{
    	DIM = atoi(argv[5]);
    }
	cout << sizeof(Simplex) << " " << sizeof(Vertex) << " " << sizeof(Edge) << " "
		 << sizeof(Triangle) << endl;
	cout << sizeof(Simplex*) << " " << sizeof(Vertex*) << " " << sizeof(Edge*) << " "
		 << sizeof(Triangle*) << endl;

	// 2D defined by edge. 3D defined by vertex.
	// Change this if necessary [FromFile: Triangle defined by edge][FromFile2: by vertex]
	// 3D support binary input, 2D does not.
	cout << "Reading in simplicial complex...\n";
	Simplicial2Complex K;
	if (DIM == 3)
		K.buildComplexFromFile2_BIN(argv[1]);
	else
		K.buildComplexFromFile(argv[1]);
	// K.outputComplex("testcomplex0.txt");
	cout << "Done\n";
	cout.flush();


	// Build psudo morse function
	cout << "Building pseudo-Morse function...\n";
	K.buildPsuedoMorseFunction();
	cout << "Done\n";
	cout.flush();


	cout << "Building filtration...\n";
	PersistencePairs P(&K);
	// P.buildFiltration();
	P.buildFiltrationWithLowerStar();
	if (DEBUG){
		ofstream filtration("filtration.txt", ios_base::out | ios_base::trunc);
		filtration << setprecision(16);
		for (vector<Simplex*>::iterator it = P.filtration.begin(); it != P.filtration.end(); it++){
			filtration << (*it)->funcValue << " ";
			if ((*it)->dim == 0){
				filtration << "Vertex " << ((Vertex*)(*it))->getoriPosition() << endl;
			}
			else if ((*it)->dim == 1){
				filtration <<"Edge " << ((Edge*)(*it))->getEPosition() << endl;
			}
			else if ((*it)->dim == 2){
				filtration << "Triangle " << ((Triangle*)(*it))->getTPosition() << endl;
			}
		}
		filtration.close();
	}
	cout << "Done\n";
	cout.flush();



	cout << "Computing persistence pairs...\n";
	P.PhatPersistence();
	if (DEBUG){
		cout << "Outputing persistence pairs...\n";
		P.outputPersistencePairs("persistencePairs.txt", DEBUG);
		cout << "Done\n";
	}
	cout << "Done!\n";
	cout.flush();


	delta = atof(argv[4]);
	cout << "Cancelling persistence pairs with delta " << delta << "\n";
	P.cancelPersistencePairs(delta);
	cout << "Done\n";
	cout.flush();

	K.outputArcs(argv[2], argv[3]);
}
