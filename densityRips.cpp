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

	//if (MAKE_SINGLE_SUBSETS) {
	//	cout << "Making subsets...\n";
	//	ifstream lengths("lengths.txt");
	//	ifstream points("pca_10.txt");
	//	ifstream funcValues("kde_10.txt");
	//	//for each digit 0-9
	//	for (int i = 0; i < 10; i++) {
	//		ofstream subset("subset" + std::to_string(i) + ".txt");
	//		ofstream subset_kde("subset" + std::to_string(i) + "_kde.txt");
	//		ofstream subset_labels("subset" + std::to_string(i) + "_labels.txt");
	//		subset << setprecision(20);
	//		subset_kde << setprecision(20);
	//		int length;
	//		lengths >> length;
	//		for (int j = 0; j < length; j++) {
	//			double coords[10];
	//			double funcValue;
	//			for (int k = 0; k < 10; k++) {
	//				points >> coords[k];
	//			}
	//			funcValues >> funcValue;
	//			for (int k = 0; k < 10; k++) {
	//				subset << coords[k] << " ";
	//			}
	//			subset << "\n";
	//			subset_kde << funcValue << "\n";
	//			subset_labels << i << "\n";
	//		}
	//		subset.close();
	//		subset_kde.close();
	//	}
	//	lengths.close();
	//	points.close();
	//	funcValues.close();
	//	cout << "Done\n";
	//}

	if (BUILD_SUBSET_MAP) {
		cout << "Building subset map...\n";
		ifstream subset("subset_pca.txt");
		ifstream full("pca_10.txt");
		ofstream subsetMap("subset_map.txt");
		subsetMap << "indices\n";
		double *points[70000];
		cout << "\tReading points...\n";
		for (int i = 0; i < 70000; i++) {
			double *coords = new double[10];
			for (int j = 0; j < 10; j++) {
				full >> coords[j];
			}
			points[i] = coords;
		}
		cout << "\tDone\n";
		int count = 0;
		while (!subset.eof()) {
			cout << "\tHandling subset point " << count << "\n";
			double coords[10];
			for (int j = 0; j < 10; j++) {
				subset >> coords[j];
				if (subset.eof()) break;
			}
			if (subset.eof()) break;
			int index = 0;
			bool found = false;
			while (!found) {
				for (int j = 0; j < 10; j++) {
					if (coords[j] != points[index][j]) {
						index++;
						break;
					}
					else if (j == 9) {
						found = true;
					}
				}
			}
			subsetMap << index << "\n";
			cout << "\tDone\n";
			count++;
		}
		cout << "Done\n";
	}

	/*if (RECONSTRUCT_IMAGES) {
		ifstream images("digits.txt");
		int count = 0;
		while (!images.eof()) {
			char byteArray[784];
			for (int i = 0; i < 784; i++) {
				int intensity;
				images >> intensity;
				if (images.eof()) break;
				byteArray[i] = (char)intensity;
			}
			if (images.eof()) break;
			bitmap_image img(28, 28);
			for (int i = 0; i < 784; i++) {
				img.set_pixel(i % 28, i / 28, byteArray[i], byteArray[i], byteArray[i]);
			}
			img.save_image("IMAGES\\IMG" + std::to_string(count) + ".bmp");
			count++;
		}
	}*/



	double RADIUS = .0;
	Simplicial2Complex K;
	if (REBUILD_COMPLEX)
	{
		ifstream coordinates("pca_10.txt");
		ifstream funcValues("kde_10.txt");
		int vPosition = 0;
		cout << "Reading density data...\n";

		while (!coordinates.eof()) {
			double coords[10];
			coordinates >> coords[0] >> coords[1] >> coords[2] >> coords[3] >> coords[4] >> coords[5] >> coords[6] >> coords[7] >> coords[8] >> coords[9];
			double funcValue;
			funcValues >> funcValue;
			if (coordinates.eof() || funcValues.eof()) break;
#ifdef THRESHOLD
			if (funcValue < THRESHOLD) continue;
#endif
			Vertex *v = new Vertex(coords, funcValue);
			vPosition = K.addVertex(v);
			v->setVposition(vPosition);
		}

		coordinates.close();
		funcValues.close();
		cout << "Done\n";
		
		/*cout << "Removing duplicate vertices...\n";
		vector<Vertex*> *unique = new vector<Vertex*>();
		for(vector<Vertex*>::iterator it = K.vBegin(); it != K.vEnd(); it++){
			if(unique->empty()){
				unique->push_back(*it);
			}else{
				bool insert = true;
				for(vector<Vertex*>::iterator itU = unique->begin(); itU != unique->end(); itU++){
					if((*it)->getFuncValue() == (*itU)->getFuncValue()){
						insert = false;
						break;
					}
				}
				if(insert == true){
					unique->push_back(*it);
				}
			}
		}
		K.setVertices(unique);
		cout << "Done\n";*/
		
		cout << "Cutting down data set...\n";
		vector<Vertex*>* newVertList = new vector<Vertex*>();
		for(int i = 0; i < K.order();i+=3){
			Vertex *v = K.getVertex(i);
			newVertList->push_back(v);
			v->setVposition(newVertList->size() - 1);
			//cout << newVertList->size() -1 << endl;
		}
		K.setVertices(newVertList);
		cout << "Done\n";

		cout << "Building Rips complex...\n";
		K.buildRipsComplex(RADIUS, 0.0);
		cout << "Done\n";

		

		cout << "Outputting Rips complex...\n";
		K.outputComplex("complex.txt");
		cout << "Done\n";
	}
	else {
		cout << "Reading in simplicial complex...\n";
		K.buildComplexFromFile2(argv[1]); //<<<<<<<<<<< Change this if necessary [File: Triangle defined by edge][File2: by vertex]
		K.outputComplex("testcomplex0.txt");
		cout << "Done\n";
	}
	
	K.flipAndTranslateVertexFunction();

	/*if (REBUILD_COMPLEX_DISCRETE) {
		cout << "Building complex from discrete metric space...\n";
		K.buildComplexDiscrete("liam-20newsgroups\\cosine_similarity.mat", "liam-20newsgroups\\kde.mat", RADIUS);
		cout << "Done\n";

		cout << "Outputting Rips complex...\n";
		K.outputComplex("complex.txt");
		cout << "Done\n";
	}
	else {
		cout << "Reading in simplicial complex...\n";
		K.buildComplexFromFile("complex.txt");
		cout << "Done\n";
	}*/

	cout << "Building pseudo-Morse function...\n";
	K.buildPsuedoMorseFunction();
	cout << "Done\n";

/*
	cout << "Outputing simplices...\n";
	K.outputVertices("vertices.txt");
	K.outputEdges("edges.txt");
	K.outputTriangles("triangle.txt");
	cout << "Done\n";*/

	/*output vertices and triangles for MATLAB*/
	//cout << "Outputing vertices and triangles for MATLAB...\n";
	//ofstream MLvertoutput("MATLABvertices.txt", ios_base::trunc | ios_base::out);
	//MLvertoutput << "xCoord yCoord funcValue\n";
	//for (vector<Vertex*>::iterator it = K.vBegin(); it != K.vEnd(); it++){
	//	Vertex* v = *it;
	//	double xCoord = get<0>(v->getCoords());
	//	double yCoord = get<1>(v->getCoords());
	//	double funcValue = v->getFuncValue();
	//	MLvertoutput << xCoord << " " << yCoord << " " << funcValue << "\n";
	//}
	//MLvertoutput.close();
	//ofstream MLtrioutput("MATLABtriangles.txt", ios_base::trunc | ios_base::out);
	//MLtrioutput << "vert1 vert2 vert3 funcValue\n";
	//for (vector<Triangle*>::iterator it = K.tBegin(); it != K.tEnd(); it++){
	//	Triangle *t = *it;
	//	Vertex *v1 = get<0>(t->getVertices());
	//	Vertex *v2 = get<1>(t->getVertices());
	//	Vertex *v3 = get<2>(t->getVertices());

	//	/*Add 1 to each position because MATLAB uses 1-indexing*/
	//	MLtrioutput << v1->getVPosition() + 1 << " " << v2->getVPosition() + 1 << " " << v3->getVPosition() + 1 << " " << t->getFuncValue() << "\n";
	//}
	//MLtrioutput.close();
	//cout << "Done\n";



	
	PersistencePairs P(&K);
	cout << "Building filtration...\n";
	P.buildFiltration();
	cout << "Done\n";
	
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


	cout << "Computing persistence pairs...\n";
	P.computePersistencePairsWithClear();
	cout << "Done!\n";

	/*cout << "Outputing persistence pairs...\n";
	P.outputPersistencePairs("persistencePairs.txt");
	cout << "Done\n";*/

	cout << "Cancelling persistence pairs...\n";
	P.cancelPersistencePairs(atof(argv[4]));
	cout << "Done\n";
	
	K.outputArcs(argv[2], argv[3]);
	
/*
	ofstream triangleIndices("MATLABtriangles.txt", ios_base::trunc | ios_base::out);
	triangleIndices << "v1 v2 v3\n";
	for (vector<Triangle*>::iterator it = K.tBegin(); it != K.tEnd(); it++) {
		Triangle *t = *it;
		Vertex *v1 = get<0>(t->getVertices());
		Vertex *v2 = get<1>(t->getVertices());
		Vertex *v3 = get<2>(t->getVertices());

		triangleIndices << v1->getVPosition() + 1 << " " << v2->getVPosition() + 1 << " " << v3->getVPosition() + 1 << "\n";
	}

	ofstream vertexIndices("vertices.txt", ios_base::trunc | ios_base::out);
	ofstream edgeIndices("MATLABedges.txt", ios_base::trunc | ios_base::out);
	ofstream particularEdges("MATLABparticularedges.txt", ios_base::trunc | ios_base::out);
	vertexIndices << "index\n";
	edgeIndices << "index1 index2\n";
	set<Simplex*> *arcs = new set<Simplex*>();
	for (unordered_set<Simplex*>::iterator it = K.cBegin(); it != K.cEnd(); it++) {
		Simplex *s = *it;
		if (s->dim == 1) {
			Edge *e = (Edge*)s;
			set<Simplex*> *descMan = K.descendingManifold(s);
			arcs->insert(descMan->begin(), descMan->end());
			delete descMan;
		}
	}

	if (arcs->empty()) {
		cout << "No critical edges\n";
	}
	ofstream manifoldtest("manifoldtest.txt");
	for (set<Simplex*>::iterator it = arcs->begin(); it != arcs->end(); it++) {
		Simplex *s = *it;
		if (s->dim == 0) {
			Vertex *v = (Vertex*)s;
			int vPosition = v->getVPosition();
			vertexIndices << vPosition + 1 << "\n";
		}
		else {
			Edge *e = (Edge*)s;
			Vertex *v1 = get<0>(e->getVertices());
			Vertex *v2 = get<1>(e->getVertices());
			int vPosition1 = v1->getVPosition();
			int vPosition2 = v2->getVPosition();
			edgeIndices << vPosition1 + 1 << " " << vPosition2 + 1 << "\n";
			double *coords1 = v1->getCoords();
			double *coords2 = v2->getCoords();
			for (int i = 0; i < 10; i++) {
				manifoldtest << coords1[i] << " ";
			}
			for (int i = 0; i < 10; i++) {
				manifoldtest << coords2[i] << " ";
			}
			manifoldtest << "\n";
		}
	}
	vertexIndices.close();
	edgeIndices.close();

	ofstream complextest("complextest.txt");
	complextest << K.order() << " " << 10 << "\n";
	for (vector<Vertex*>::iterator it = K.vBegin(); it != K.vEnd(); it++) {
		Vertex *v = *it;
		double *coords = v->getCoords();
		for (int i = 0; i < 10; i++) {
			complextest << coords[i] << " ";
		}
		complextest << v->getFuncValue() << "\n";
	}
	int numOfTriangles = 0;
	for (vector<Triangle*>::iterator it = K.tBegin(); it != K.tEnd(); it++) {
		numOfTriangles++;
	}
	complextest << numOfTriangles << "\n";
	for (vector<Triangle*>::iterator it = K.tBegin(); it != K.tEnd(); it++) {
		Triangle *t = *it;
		Vertex *v1 = get<0>(t->getVertices());
		Vertex *v2 = get<1>(t->getVertices());
		Vertex *v3 = get<2>(t->getVertices());
		complextest << v1->getVPosition() << " " << v2->getVPosition() << " " << v3->getVPosition() << "\n";
	}

	char temp;



	cin >> temp;
	*/

	//ofstream criticalV("criticalVertices.txt");
	//ofstream criticalE("criticalEdges.txt");
	//ofstream criticalT("criticalTris.txt");
	//criticalV << "x y funcValue\n";
	//criticalE << "x1 y1 x2 y2 funcValue\n";
	////criticalT << "x1 y1 x2 y2 x3 y3 funcValue\n";
	//for (unordered_set<Simplex*>::iterator it = K.cBegin(); it != K.cEnd(); it++){
	//	Simplex* s = *it;
	//	if (s->dim == 0){
	//		Vertex *v = (Vertex*)s;
	//		criticalV << get<0>(v->getCoords()) << " " << get<1>(v->getCoords()) << " " << v->getFuncValue() << "\n";
	//	}
	//	else if (s->dim == 1){
	//		Edge *e = (Edge*)s;
	//		Vertex *v1 = get<0>(e->getVertices());
	//		Vertex *v2 = get<1>(e->getVertices());
	//		criticalE << get<0>(v1->getCoords()) << " " << get<1>(v1->getCoords()) << " " << get<0>(v2->getCoords()) << " " << get<1>(v2->getCoords()) << " " << e->getFuncValue() << "\n";
	//	}
	//	/*else{
	//		Triangle *t = (Triangle*)s;
	//		Vertex *v1 = get<0>(t->getVertices());
	//		Vertex *v2 = get<1>(t->getVertices());
	//		Vertex *v3 = get<2>(t->getVertices());
	//		criticalT << get<0>(v1->getCoords()) << " " << get<1>(v1->getCoords()) << " " << get<0>(v2->getCoords()) << " " << get<1>(v2->getCoords()) << " " << get<0>(v3->getCoords()) << " " << get<1>(v3->getCoords()) << " " << t->getFuncValue() << "\n";
	//	}*/
	//}
	//criticalV.close();
	//criticalE.close();
	//criticalT.close();

	/*ofstream MLedges("MATLABedges.txt");
	MLedges << "x1 y1 x2 y2 funcValue\n";
	set<Edge*> *edges = new set<Edge*>();
	for (unordered_set<Simplex*>::iterator it = K.cBegin(); it != K.cEnd(); it++){
		Simplex *s = *it;
		if (s->dim == 1){
			set<Simplex*>* unstable1Manifold = K.descendingManifold(s);
			for (set<Simplex*>::iterator it2 = unstable1Manifold->begin(); it2 != unstable1Manifold->end(); it2++){
				if ((*it2)->dim == 1){
					edges->insert((Edge*)*it2);
				}
			}
		}
	}
	for (set<Edge*>::iterator it = edges->begin(); it != edges->end(); it++){
		Edge *e = *it;
		Vertex *v1 = get<0>(e->getVertices());
		Vertex *v2 = get<1>(e->getVertices());
		MLedges << get<0>(v1->getCoords()) << " " << get<1>(v1->getCoords()) << " " << get<0>(v2->getCoords()) << " " << get<1>(v2->getCoords()) << " " << e->getFuncValue() << "\n";
	}
	MLedges.close();*/
}
