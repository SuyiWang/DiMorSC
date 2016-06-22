#include <vector>
#include <tuple>
#include <unordered_set>
#include <set>
#include <stack>
//#include <ANN/ANN.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <unordered_map>
#include <map>
//#include "C:\Program Files\MATLAB\R2015b\extern\include\mat.h"
#define DIM 3
//#define min(a,b) (a <= b)? a : b
using namespace std;

/*namespace std{
	namespace
	{

		// Code from boost
		// Reciprocal of the golden ratio helps spread entropy
		//     and handles duplicates.
		// See Mike Seymour in magic-numbers-in-boosthash-combine:
		//     http://stackoverflow.com/questions/4948780

		template <class T>
		inline void hash_combine(std::size_t& seed, T const& v)
		{
			seed ^= std::hash<T>()(v)+0x9e3779b9 + (seed << 6) + (seed >> 2);
		}

		// Recursive template code derived from Matthieu M.
		template <class Tuple, size_t Index = std::tuple_size<Tuple>::value - 1>
		struct HashValueImpl
		{
			static void apply(size_t& seed, Tuple const& tuple)
			{
				HashValueImpl<Tuple, Index - 1>::apply(seed, tuple);
				hash_combine(seed, std::get<Index>(tuple));
			}
		};

		template <class Tuple>
		struct HashValueImpl<Tuple, 0>
		{
			static void apply(size_t& seed, Tuple const& tuple)
			{
				hash_combine(seed, std::get<0>(tuple));
			}
		};
	}

	template <typename ... TT>
	struct hash<std::tuple<TT...>>
	{
		size_t
			operator()(std::tuple<TT...> const& tt) const
		{
			size_t seed = 0;
			HashValueImpl<std::tuple<TT...> >::apply(seed, tt);
			return seed;
		}

	};
}*/

class Simplex;
class Simplicial2Complex;
class Triangle;
class Edge;
class Vertex;
class DiscreteVField;

class Simplex{
public:
	/*Simplex may be cast to Vertex, Edge, or Triangle based on whether
	dim = 0, 1, or 2, respectively*/
	int dim;
	double funcValue;
	static bool simplexPointerCompare(const Simplex* s, const Simplex* t);
	unsigned int filtrationPosition;
};

class Simplicial2Complex{
	vector<Vertex*> vertexList;
	vector<Edge*> edgeList;
	vector<Triangle*> triList;
	unordered_set<Simplex*> criticalSet;
	DiscreteVField *V;

public:
	Simplicial2Complex();
	void reserveVertexVector(int size){
		this->vertexList.reserve(size);
	}
	void reserveEdgeVector(int size){
		this->edgeList.reserve(size);
	}
	void reserveTriangleVector(int size){
		this->triList.reserve(size);
	}

	int addVertex(Vertex *v);
	int setVertices(vector<Vertex*>* list){
		this->vertexList = *list;
	}
	int addEdge(Edge *e);
	int addTriangle(Triangle* t);
	bool multiEdgeCheck();
	bool multiTriangleCheck();
	Vertex* getVertex(int position);
	Edge* getEdge(int position);
	Triangle* getTriangle(int position);
	void addCriticalPoint(Simplex *s);
	void removeCriticalPoint(Simplex *s);
	void setDiscreteVField(DiscreteVField *V);
	int order();
	void buildRipsComplex(double radius, double eps);
	void outputArcs(string, string);
	void buildComplexFromFile(string pathname);
	void buildComplexFromFile2(string pathname);
	void buildComplexDiscrete(string, string, double);
	void outputComplex(string pathname);
	void buildPsuedoMorseFunction();
	void cancel0PersistencePairs();
	void cancel0PersistencePairs2();
	/*Requires: s is a critical simplex*/
	set<Simplex*>* ascendingManifold(Simplex *s);
	set<Simplex*>* descendingManifold(Simplex *s);
	void flipAndTranslateVertexFunction();
	DiscreteVField* getDiscreteVField();
	vector<Vertex*>::iterator vBegin(){
		return this->vertexList.begin();
	}
	vector<Vertex*>::iterator vEnd(){
		return this->vertexList.end();
	}
	vector<Edge*>::iterator eBegin(){
		return this->edgeList.begin();
	}
	vector<Edge*>::iterator eEnd(){
		return this->edgeList.end();
	}
	vector<Triangle*>::iterator tBegin(){
		return this->triList.begin();
	}
	vector<Triangle*>::iterator tEnd(){
		return this->triList.end();
	}
	unordered_set<Simplex*>::iterator cBegin(){
		return this->criticalSet.begin();
	}
	unordered_set<Simplex*>::iterator cEnd(){
		return this->criticalSet.end();
	}
	void outputVertices(string pathname);
	void outputEdges(string pathname);
	void outputTriangles(string pathname);
};

class Vertex:public Simplex{
	double coords[DIM];
	vector<Edge*> incidenceList;
	int vPosition;
public:
	Vertex(double *coords, double funcValue);
	Vertex(double funcValue) {
		this->funcValue = funcValue;
		this->dim = 0;
	}
	void addEdge(Edge* e);
	void removeEdge(Edge* e);
	Vertex* getAdjacentVertex(Edge* e);
	bool isAdjacent(Vertex *v);
	bool isDeleted();
	void deleteVertex();
	double* getCoords();
	double getFuncValue();
	void setFuncValue(double funcValue){
		this->funcValue = funcValue;
	}
	int degree();
	int getVPosition();
	void setVposition(int p);
	Edge* findEdge(Vertex *v);
	bool hasEdge(Edge *e){
		for(vector<Edge*>::iterator it = this->incidenceList.begin(); it != this->incidenceList.end(); it++){
			if ((*it) == e){
				return true;
			}
		}
		return false;
	}
	vector<Edge*>::iterator begin(){
		return this->incidenceList.begin();
	}
	vector<Edge*>::iterator end(){
		return this->incidenceList.end();
	}
	static bool vertexPointerCompare(Vertex *v1, Vertex *v2){
		if (v1->funcValue < v2->funcValue){
			return true;
		}
		else if (v1->funcValue == v2->funcValue){
			return v1->vPosition < v2->vPosition;
		}
		else{
			return false;
		}
	}
};

class Edge:public Simplex{
	tuple<Vertex*, Vertex*> vertices;
	vector<Triangle*> incidenceList;
	int ePosition;
	/*This symbolic perturbation represents a multiple of some infinitesimal eps>0
	This is used to distinguish between edges with the same function value. It is equal to the
	sum of the function values of the two vertices*/
	double symbolicPerturbation1;

public:
	Edge(Vertex* v1, Vertex* v2);
	Edge(Vertex* v1, Vertex* v2, double funcValue);
	tuple<Vertex*, Vertex*> getVertices();
	void addTriangle(Triangle *t);
	void removeTriangle(Triangle *t);
	/*requires Triangle *t is in this->incidenceList*/
	Vertex* oppsiteVertex(Triangle *t);
	bool hasTriangle(Triangle *t){
		for (vector<Triangle*>::iterator it = this->incidenceList.begin(); it != this->incidenceList.end(); it++){
			if ((*it) == t){
				return true;
			}
		}
		return false;
	}
	int degree();
	bool isDeleted();
	void deleteEdge();
	double getFuncValue();
	void setFuncValue(double funcValue);
	int getEPosition();
	void setEposition(int p);
	/*requires that this and e share a Vertex*/
	Vertex* findVertex(Edge* e);
	double getSymPerturb(){ return this->symbolicPerturbation1; }
	void setSymPerturb(double f){ this->symbolicPerturbation1 = f; }
	vector<Triangle*>::iterator begin(){
		return this->incidenceList.begin();
	}
	vector<Triangle*>::iterator end(){
		return this->incidenceList.end();
	}
	static bool edgePointerCompare(Edge *e1, Edge *e2){
		if (e1->funcValue < e2->funcValue){
			return true;
		}
		else if (e1->funcValue == e2->funcValue){
			if (e1->symbolicPerturbation1 < e2->symbolicPerturbation1){
				return true;
			}
			else if (e1->symbolicPerturbation1 == e2->symbolicPerturbation1){
				return e1->ePosition < e2->ePosition;
			}
			else{
				return false;
			}
		}
		else{
			return false;
		}
	}
};

class Triangle:public Simplex{
	tuple<Edge*, Edge*, Edge*> edges;
	tuple<Vertex*, Vertex*, Vertex*> vertices;
	int tPosition;
	/*These symbolic perturbations represent multiples of some infinitesimal eps > 0, where
	f(this) = this->funcValue + symbolicPerturbation1 * eps + symbolicPerturbation2 * eps^2
	This is used to raise the triangle above its edges and distinguish it from other triangles
	with the same function value.
	symbolicPerturbation1 = the symbolic perturbation of the edge with the highest function value
	symbolicPerturbation2 = the sum of the vertices*/
	double symbolicPerturbation1;
	double symbolicPerturbation2;

public:
	Triangle(Edge *e1, Edge *e2, Edge *e3);
	Triangle(Vertex *v1, Vertex *v2, Vertex *v3);
	Triangle(Edge *e1, Edge *e2, Edge *e3, double funcValue);
	Triangle(Vertex *v1, Vertex *v2, Vertex *v3, double funcValue);
	tuple<Edge*, Edge*, Edge*> getEdges();
	tuple<Vertex*, Vertex*, Vertex*> getVertices();
	bool isDeleted();
	void deleteTriangle();
	double getFuncValue();
	void setFuncValue(double funcValue);
	int getTPosition();
	void setTposition(int p);
	Vertex* oppositeVertex(Edge* e);
	tuple<double, double> getSymPerturb(){ return make_tuple(this->symbolicPerturbation1, this->symbolicPerturbation2); }
	void setSymPerturb(double f1, double f2){ this->symbolicPerturbation1 = f1; this->symbolicPerturbation2 = f2; }
	static bool triPointerCompare(Triangle *t1, Triangle *t2){
		if (t1->funcValue < t2->funcValue){
			return true;
		}
		else if (t1->funcValue == t2->funcValue){
			if (t1->symbolicPerturbation1 < t2->symbolicPerturbation1){
				return true;
			}
			else if (t1->symbolicPerturbation1 == t2->symbolicPerturbation1){
				if (t1->symbolicPerturbation2 < t2->symbolicPerturbation2){
					return true;
				}
				else if (t1->symbolicPerturbation2 == t2->symbolicPerturbation2){
					return t1->tPosition < t2->tPosition;
				}
				else{
					return false;
				}
			}
			else{
				return false;
			}
		}
		else{
			return false;
		}
	}
};

class DiscreteVField{
private:
	unordered_map<Vertex*, Edge*> *VEmap;
	unordered_map<Edge*, Triangle*> *ETmap;

public:
	DiscreteVField();
	bool containsPair(Vertex *v, Edge *e);
	bool containsPair(Edge *e, Triangle *t);
	void addPair(Vertex *v, Edge *e);
	void addPair(Edge *e, Triangle *t);
	void removePair(Vertex *v, Edge *e);
	void removePair(Edge *e, Triangle *t);
	void outputVEmap();
	void outputETmap();
};

Simplicial2Complex::Simplicial2Complex(){
	this->V = new DiscreteVField();
}

int Simplicial2Complex::addVertex(Vertex *v){
	int position = this->vertexList.size();
	this->vertexList.push_back(v);
	return position;
}

int Simplicial2Complex::addEdge(Edge *e){
	int position = this->edgeList.size();
	this->edgeList.push_back(e);
	return position;
}

int Simplicial2Complex::addTriangle(Triangle *t){
	int position = this->triList.size();
	this->triList.push_back(t);
	return position;
}

bool Simplicial2Complex::multiEdgeCheck() {
	for (int i = 0; i < this->edgeList.size(); i++) {
		if (i % 1000 == 0) cout << i << "\n";
		for (int j = i + 1; j < this->edgeList.size(); j++) {
			Edge *e1 = this->edgeList[i];
			Edge *e2 = this->edgeList[j];
			Vertex *v11 = get<0>(e1->getVertices());
			Vertex *v12 = get<1>(e1->getVertices());
			Vertex *v21 = get<0>(e2->getVertices());
			Vertex *v22 = get<1>(e2->getVertices());
			if (v11 == v21 && v12 == v22 || v11 == v22 && v12 == v21) {
				return true;
			}
		}
	}
	return false;
}

Vertex* Simplicial2Complex::getVertex(int position){
	return this->vertexList.at(position);
}

Edge* Simplicial2Complex::getEdge(int position){
	return this->edgeList.at(position);
}

Triangle* Simplicial2Complex::getTriangle(int position){
	return this->triList.at(position);
}

void Simplicial2Complex::addCriticalPoint(Simplex *s){
	this->criticalSet.insert(s);
}

void Simplicial2Complex::removeCriticalPoint(Simplex *s){
	this->criticalSet.erase(s);
}

void Simplicial2Complex::setDiscreteVField(DiscreteVField *V){
	this->V = V;
}

DiscreteVField* Simplicial2Complex::getDiscreteVField(){
	return this->V;
}

int Simplicial2Complex::order(){
	return this->vertexList.size();
}

void Simplicial2Complex::outputComplex(string pathname) {
	ofstream file(pathname);
	file << std::setprecision(20);
	file << this->vertexList.size() << "\n";
	for (vector<Vertex*>::iterator it = this->vertexList.begin(); it != this->vertexList.end(); it++) {
		Vertex *v = *it;
		double *coords = v->getCoords();
		double funcValue = v->getFuncValue();
		for (int j = 0; j < DIM; j++) {
			file << coords[j] << " ";
		}
		file << funcValue << "\n";
	}

	file << this->edgeList.size() << "\n";
	for (vector<Edge*>::iterator it = this->edgeList.begin(); it != this->edgeList.end(); it++) {
		Edge *e = *it;
		Vertex *v1 = get<0>(e->getVertices());
		Vertex *v2 = get<1>(e->getVertices());
		file << v1->getVPosition() << " " << v2->getVPosition() << "\n";
	}

	file << this->triList.size() << "\n";
	for (vector<Triangle*>::iterator it = this->triList.begin(); it != this->triList.end(); it++) {
		Triangle *t = *it;
		Edge *e1 = get<0>(t->getEdges());
		Edge *e2 = get<1>(t->getEdges());
		Edge *e3 = get<2>(t->getEdges());
		file << e1->getEPosition() << " " << e2->getEPosition() << " " << e3->getEPosition() << "\n";
	}

	file.close();
}

void Simplicial2Complex::outputArcs(string vertexFile, string edgeFile){
	ofstream vFile(vertexFile);
	ofstream eFile(edgeFile);
	set<Simplex*> manifolds;
	for(unordered_set<Simplex*>::iterator it = this->cBegin(); it != this->cEnd(); it++){
		Simplex *s = *it;
		if(s->dim == 1){
			Edge *e = (Edge*)s;
			set<Simplex*> *manifold = this->descendingManifold((Simplex*)e);
			for (set<Simplex*>::iterator it2 = manifold->begin(); it2 != manifold->end(); it2++){
				manifolds.insert(*it2);
			}
			delete manifold;
		}
	}
	vector<Vertex*> vertices;
	vector<Edge*>edges;
	for(set<Simplex*>::iterator it = manifolds.begin(); it != manifolds.end(); it++){
		Simplex* s = *it;
		if(s->dim == 0){
			vertices.push_back((Vertex*)s);
		}
		else{
			edges.push_back((Edge*)s);
		}
	}
	std::map<Vertex*, int> map;
	for(int i = 0; i < vertices.size(); i++){
		Vertex *v = vertices[i];
		map.insert( std::pair<Vertex*,int>(v, i + 1));
		for(int j = 0; j < DIM; j++){
			vFile << v->getCoords()[j] << " ";
		}
		vFile << v->getFuncValue() << " ";
		if (this->criticalSet.count((Simplex*)v) > 0){
			vFile << "0";
		}else{
			vFile << "-1";
		}
		vFile << endl;
	}
	for(int i = 0; i < edges.size(); i++){
		Edge *e = edges[i];
		Vertex* v1 = get<0>(e->getVertices());
		Vertex* v2 = get<1>(e->getVertices());
		eFile << map.find(v1)->second << " " << map.find(v2)->second << " ";
		if(this->criticalSet.count((Simplex*)e) > 0){
			eFile << "1";
		}else{
			eFile << "-1";
		}
		eFile << endl;
	}
}

void Simplicial2Complex::buildComplexFromFile(string pathname) {
	ifstream file(pathname);
	int numOfVertices;
	file >> numOfVertices;
	for (int i = 0; i < numOfVertices; i++) {
		double coords[DIM];
		double funcValue;
		for (int j = 0; j < DIM; j++) {
			file >> coords[j];
		}
		file >> funcValue;
		Vertex *v = new Vertex(coords, funcValue);
		int vPosition = this->addVertex(v);
		v->setVposition(vPosition);
	}

	int numOfEdges;
	file >> numOfEdges;
	for (int i = 0; i < numOfEdges; i++) {
		int vIndex1, vIndex2;
		file >> vIndex1 >> vIndex2;
		Vertex *v1 = this->getVertex(vIndex1);
		Vertex *v2 = this->getVertex(vIndex2);
		Edge *e = new Edge(v1, v2);
		v1->addEdge(e);
		v2->addEdge(e);
		int ePosition = this->addEdge(e);
		e->setEposition(ePosition);
	}

	int numOfTris;
	file >> numOfTris;
	for (int i = 0; i < numOfTris; i++) {
		int eIndex1, eIndex2, eIndex3;
		file >> eIndex1 >> eIndex2 >> eIndex3;
		Edge *e1 = this->getEdge(eIndex1);
		Edge *e2 = this->getEdge(eIndex2);
		Edge *e3 = this->getEdge(eIndex3);
		Triangle *t = new Triangle(e1, e2, e3);
		e1->addTriangle(t);
		e2->addTriangle(t);
		e3->addTriangle(t);
		int tPosition = this->addTriangle(t);
		t->setTposition(tPosition);
	}

	file.close();
}

void Simplicial2Complex::buildComplexFromFile2(string pathname) {
	ifstream file(pathname);
	int numOfVertices;
	file >> numOfVertices;
	for (int i = 0; i < numOfVertices; i++) {
		double coords[DIM];
		double funcValue;
		for (int j = 0; j < DIM; j++) {
			file >> coords[j];
		}
		file >> funcValue;
		Vertex *v = new Vertex(coords, funcValue);
		int vPosition = this->addVertex(v);
		v->setVposition(vPosition);
	}

	int numOfEdges;
	file >> numOfEdges;
	for (int i = 0; i < numOfEdges; i++) {
		int vIndex1, vIndex2;
		file >> vIndex1 >> vIndex2;
		Vertex *v1 = this->getVertex(vIndex1);
		Vertex *v2 = this->getVertex(vIndex2);
		Edge *e = new Edge(v1, v2);
		v1->addEdge(e);
		v2->addEdge(e);
		int ePosition = this->addEdge(e);
		e->setEposition(ePosition);
	}

	int numOfTris;
	file >> numOfTris;
	for (int i = 0; i < numOfTris; i++) {
		int vIndex1, vIndex2, vIndex3;
		file >> vIndex1 >> vIndex2 >> vIndex3;
		Vertex *v1 = this->getVertex(vIndex1);
		Vertex *v2 = this->getVertex(vIndex2);
		Vertex *v3 = this->getVertex(vIndex3);
		Triangle *t = new Triangle(v1, v2, v3);
		Edge *e1 = get<0>(t->getEdges());
		Edge *e2 = get<1>(t->getEdges());
		Edge *e3 = get<2>(t->getEdges());
		e1->addTriangle(t);
		e2->addTriangle(t);
		e3->addTriangle(t);
		int tPosition = this->addTriangle(t);
		t->setTposition(tPosition);
	}

	file.close();
}

/*
void Simplicial2Complex::buildRipsComplex(double radius, double eps){
	this->edgeList.clear();
	this->triList.clear();
	
	//Build kd-tree
	int numOfVertices = this->vertexList.size();
	ANNpointArray data = annAllocPts(numOfVertices, DIM);
	for (int i = 0; i < numOfVertices; i++){
		ANNpoint currentPoint = data[i];
		Vertex *v = this->vertexList.at(i);
		double *coords = v->getCoords();
		for (int j = 0; j < DIM; j++) {
			currentPoint[j] = coords[j];
		}
	}

	ANNkd_tree *kdTree = new ANNkd_tree(data, numOfVertices, DIM);

	ANNdist sqRad = (double)radius * radius;
	for (int i = 0; i < numOfVertices; i++){
		if(i % 1000 == 0) cout << "Handling vertex " << i << "\n";
		ANNpoint queryPt = data[i];
		Vertex *queryVert = this->vertexList.at(i);
		ANNidxArray nnIdx = new ANNidx[100];
		ANNdistArray dists = new ANNdist[100];

		int numOfNN = kdTree->annkFRSearch(queryPt, sqRad, 100, nnIdx, dists, eps);

		for (int j = 0; j < ((numOfNN < 100)? numOfNN : 100); j++){

			Vertex *NNVert = this->vertexList.at(nnIdx[j]);
			//If NNVert is not the query Vertex and is not already adjacent
			if (!(NNVert == queryVert || queryVert->isAdjacent(NNVert))){
				//Create Edge between them
				Edge *e = new Edge(queryVert, NNVert);
				//But don't add it yet, because we need to test paths of length 2 from the query Vertex
				//If there is already a path of length 2 from the query vertex to the NNVert, then add a triangle
				for (vector<Edge*>::iterator it = queryVert->begin(); it != queryVert->end(); it++){
					Vertex *midVert = queryVert->getAdjacentVertex(*it);
					for (vector<Edge*>::iterator it2 = midVert->begin(); it2 != midVert->end(); it2++){
						Vertex *endVert = midVert->getAdjacentVertex(*it2);
						//If the end Vertex is equal to the NN Vertex, then we need to add a triangle corresponding to the path we took
						if (endVert == NNVert){
							Triangle *t = new Triangle(e, *it, *it2);
							e->addTriangle(t);
							(*it)->addTriangle(t);
							(*it2)->addTriangle(t);
							int tPosition = this->addTriangle(t);
							t->setTposition(tPosition);
						}
					}
				}

				//Now add the new Edge
				queryVert->addEdge(e);
				NNVert->addEdge(e);
				int ePosition = this->addEdge(e);
				e->setEposition(ePosition);
			}
		}
		delete[] nnIdx;
		delete[] dists;
	}
	annDeallocPts(data);
}
*/

//void Simplicial2Complex::buildComplexDiscrete(string distPathname, string funcPathname, double radius) {
//	const char *distPname = distPathname.c_str();
//	const char *funcPname = funcPathname.c_str();
//	
//	MATFile *kdeFile = matOpen(funcPname, "r");
//	mxArray *kdeValues = matGetVariable(kdeFile, "data");
//	int numOfVertices = (int)mxGetN(kdeValues);
//	double *arrayPtr = mxGetPr(kdeValues);
//
//	for (int i = 0; i < numOfVertices; i++) {
//		double funcValue = arrayPtr[i];
//		Vertex *v = new Vertex(funcValue);
//		int vPosition = this->addVertex(v);
//		v->setVposition(vPosition);
//	}
//
//	mxDestroyArray(kdeValues);
//
//	MATFile *distFile = matOpen(distPname, "r");
//	mxArray *dists = matGetVariable(distFile, "data");
//	arrayPtr = mxGetPr(dists);
//	int numOfRows = (int)mxGetM(dists);
//
//	for (int i = 0; i < numOfVertices; i++) {
//		for (int j = i + 1; j < numOfVertices; j++) {
//			double distance = 1 - arrayPtr[i + j * numOfRows];
//			if (distance < radius) {
//				Vertex *v1 = this->vertexList[i];
//				Vertex *v2 = this->vertexList[j];
//				Edge *e = new Edge(v1, v2);
//				for (vector<Edge*>::iterator it = v1->begin(); it != v1->end(); it++) {
//					Vertex *midVert = v1->getAdjacentVertex(*it);
//					for (vector<Edge*>::iterator it2 = midVert->begin(); it2 != midVert->end(); it2++) {
//						Vertex *endVert = midVert->getAdjacentVertex(*it2);
//						if (endVert == v2) {
//							Triangle *t = new Triangle(e, *it, *it2);
//							int tPosition = this->addTriangle(t);
//							t->setTposition(tPosition);
//						}
//					}
//				}
//				int ePosition = this->addEdge(e);
//				e->setEposition(ePosition);
//			}
//		}
//	}
//
//	mxDestroyArray(dists);
//	
//}

void Simplicial2Complex::buildPsuedoMorseFunction(){
	for (unsigned int i = 0; i < this->edgeList.size(); i++){
		Edge *e = this->edgeList.at(i);
		tuple<Vertex*, Vertex*> vertices = e->getVertices();
		Vertex *max = get<0>(vertices);
		Vertex *v2 = get<1>(vertices);
		double total = max->getFuncValue() + v2->getFuncValue();
		if (v2->getFuncValue() > max->getFuncValue()){
			max = v2;
		}
		e->setFuncValue(max->getFuncValue());
		e->setSymPerturb(total);
	}
	for (unsigned int i = 0; i < this->triList.size(); i++){
		Triangle *t = this->triList.at(i);
		tuple<Edge*, Edge*, Edge*> edges = t->getEdges();
		Edge *max = get<0>(edges);
		Edge *e2 = get<1>(edges);
		Edge *e3 = get<2>(edges);
		tuple<Vertex*, Vertex*, Vertex*> vertices = t->getVertices();
		double total = get<0>(vertices)->getFuncValue() + get<1>(vertices)->getFuncValue() + get<2>(vertices)->getFuncValue();
		if (e2->getFuncValue() > max->getFuncValue() || e2->getFuncValue() == max->getFuncValue() && e2->getSymPerturb() > max->getSymPerturb()){
			max = e2;
		}
		if (e3->getFuncValue() > max->getFuncValue() || e3->getFuncValue() == max->getFuncValue() && e3->getSymPerturb() > max->getSymPerturb()){
			max = e3;
		}
		t->setFuncValue(max->getFuncValue());
		t->setSymPerturb(max->getSymPerturb(), total);
	}
}

void Simplicial2Complex::cancel0PersistencePairs(){
	//Mark every simplex as a critical point first
	for (int i = 0;(unsigned) i < this->vertexList.size(); i++){
		this->criticalSet.insert((Simplex*)this->vertexList.at(i));
	}
	for (int i = 0;(unsigned) i < this->edgeList.size(); i++){
		this->criticalSet.insert((Simplex*)this->edgeList.at(i));
	}
	for (int i = 0;(unsigned) i < this->triList.size(); i++){
		this->criticalSet.insert((Simplex*)this->triList.at(i));
	}

	//Start with the vertices and edges
	for (int i = 0;(unsigned) i < this->vertexList.size(); i++){
		Vertex *v = this->vertexList.at(i);
		Vertex *min = v;
		for (vector<Edge*>::iterator j = v->begin(); j != v->end(); j++){
			Vertex *w = v->getAdjacentVertex(*j);
			if (w->getFuncValue() < v->getFuncValue()){
				min = w;
			}
		}
		if (min != v){
			Edge *e = v->findEdge(min);
			this->V->addPair(v, e);
			this->criticalSet.erase((Simplex*)v);
			this->criticalSet.erase((Simplex*)e);
		}
	}
	
	//Now do edges and triangles
	for (int i = 0; (unsigned)i < this->triList.size(); i++){
		Triangle *t = this->triList.at(i);
		tuple<Edge*, Edge*, Edge*> edges = t->getEdges();
		Edge *e1 = get<0>(edges);
		Edge *e2 = get<1>(edges);
		Edge *e3 = get<2>(edges);
		/*Only consider the critical edges*/
		vector<Edge*> criticalEdges;
		if (this->criticalSet.count((Simplex*)e1) == 1){ criticalEdges.push_back(e1); }
		if (this->criticalSet.count((Simplex*)e2) == 1){ criticalEdges.push_back(e2); }
		if (this->criticalSet.count((Simplex*)e3) == 1){ criticalEdges.push_back(e3); }
		/*Now figure out which critical edge has the greatest difference between it and the opposite vertex*/
		if (criticalEdges.size() > 0){
			double maxDif = 0;
			Edge* maxDifEdge = criticalEdges.at(0);
			for (int i = 0; (unsigned)i < criticalEdges.size(); i++){
				Edge *critEdge = criticalEdges.at(i);
				Vertex *opVertex = critEdge->oppsiteVertex(t);
				if (critEdge->getFuncValue() == t->getFuncValue() && critEdge->getFuncValue() - opVertex->getFuncValue() > maxDif){
					maxDifEdge = critEdge;
				}
			}
			/*Check again that maxDifEdge and t are a valid 0-persistence pair*/
			if (maxDifEdge->getFuncValue() == t->getFuncValue()){
				this->V->addPair(maxDifEdge, t);
				this->criticalSet.erase((Simplex*)maxDifEdge);
				this->criticalSet.erase((Simplex*)t);
			}
		}
	}
}

void Simplicial2Complex::cancel0PersistencePairs2(){
	/*Start with the vertices in ascending order*/
	vector<Vertex*> sortedV(this->vertexList);
	sort(sortedV.begin(), sortedV.end(), Vertex::vertexPointerCompare);
	for (vector<Vertex*>::iterator it = sortedV.begin(); it != sortedV.end(); it++){
		Vertex *v = *it;
		/*For each vertex v with nonzero degree, find the minimum edge differing from v by an infinitesimal*/
		if (v->degree() > 0){
			Edge *min = *(v->begin());
			for (vector<Edge*>::iterator it2 = v->begin(); it2 != v->end(); it2++){
				Edge* e = *it2;
				if (Edge::edgePointerCompare(e, min)){
					min = e;
				}
			}
			/*If the minimum edge differs infinitesimally, add the pair to the discrete gradient vector field*/
			if (min->getFuncValue() == v->getFuncValue()){
				this->V->addPair(v, min);
			}
			/*Otherwise, v is critical*/
			else{
				this->criticalSet.insert((Simplex*)v);
			}
		}
		/*If v has degree zero, then it is critical*/
		else{
			this->criticalSet.insert((Simplex*)v);
		}
	}

	/*Now do the same with the edges*/
	vector<Edge*> sortedE(this->edgeList);
	sort(sortedE.begin(), sortedE.end(), Edge::edgePointerCompare);
	for (vector<Edge*>::iterator it = sortedE.begin(); it != sortedE.end(); it++){
		Edge *e = *it;
		/*For each edge e with at least 1 cofacet, find the minimum triangle differing from e by an second order infinitesimal*/
		if (e->degree() > 0){
			Triangle *min = *(e->begin());
			for (vector<Triangle*>::iterator it2 = e->begin(); it2 != e->end(); it2++){
				Triangle *t = *it2;
				if (Triangle::triPointerCompare(t, min)){
					min = t;
				}
			}
			/*Now you have the minimum triangle in min, check if it differs only by a second order infitesimal*/
			if (min->getFuncValue() == e->getFuncValue() && get<0>(min->getSymPerturb()) == e->getSymPerturb()){
				this->V->addPair(e, min);
			}
			/*Otherwise, e is critical*/
			else{
				this->criticalSet.insert((Simplex*)e);
			}
		}
		/*If e has no cofacets, then it is critical*/
		else{
			this->criticalSet.insert((Simplex*)e);
		}
	}

	/*Finally the triangles. We don't need to sort them since any non-critical triangle is already in the discrete vector field*/
	for (vector<Triangle*>::iterator it = this->triList.begin(); it != this->triList.end(); it++){
		Triangle *t = *it;
		tuple<Edge*, Edge*, Edge*> edges = t->getEdges();
		if (!(this->V->containsPair(get<0>(edges), t) || this->V->containsPair(get<1>(edges), t) || this->V->containsPair(get<2>(edges), t))){
			/*Then t is not in any pair in V, so it is critical*/
			this->criticalSet.insert((Simplex*)t);
		}
	}
}

set<Simplex*>* Simplicial2Complex::descendingManifold(Simplex* s){
	set<Simplex*> *manifold;
	/*If s is a minimum, the only simplex in the decending manifold is s itself*/
	if (s->dim == 0){
		manifold = new set<Simplex*>({ s });
	}
	else if (s->dim == 1){
		manifold = new set<Simplex*>({ s });
		stack<Simplex*> *st = new stack<Simplex*>();
		Edge *e = (Edge*)s;
		Vertex *v1 = get<0>(e->getVertices());
		Vertex *v2 = get<1>(e->getVertices());

		st->push((Simplex*)v1);
		st->push((Simplex*)v2);
		manifold->insert((Simplex*)v1);
		manifold->insert((Simplex*)v2);
		while (!st->empty()){
			Simplex* simplex = st->top();
			st->pop();
			/*simplex could either be an edge or a vertex*/
			if (simplex->dim == 0){
				/*If it's a vertex, look at all incident edges and see where you can go. There'll be at most one direction*/
				Vertex *vertex = (Vertex*)simplex;
				for (vector<Edge*>::iterator it = vertex->begin(); it != vertex->end(); it++){
					if (this->V->containsPair(vertex, *it)){
						manifold->insert(*it);
						st->push(*it);
						break;
					}
				}
			}
			else{
				/*An edge can only be entered from a vertex, so that uses up its arrow, and it can only leave via the other vertex*/
				Edge *edge = (Edge*)simplex;
				Vertex *vert1 = get<0>(edge->getVertices());
				Vertex *vert2 = get<1>(edge->getVertices());
				Vertex *exitVert = vert1;
				if (!this->V->containsPair(vert2, edge)){
					exitVert = vert2;
				}
				manifold->insert((Simplex*)exitVert);
				st->push((Simplex*)exitVert);
			}

		}
		delete st;
	}
	else{
		manifold = new set<Simplex*>({ s });
		stack<Simplex*> *st = new stack<Simplex*>();
		Triangle *t = (Triangle*)s;
		Edge *e1 = get<0>(t->getEdges());
		Edge *e2 = get<1>(t->getEdges());
		Edge *e3 = get<2>(t->getEdges());

		st->push((Simplex*)e1);
		st->push((Simplex*)e2);
		st->push((Simplex*)e3);
		manifold->insert((Simplex*)e1);
		manifold->insert((Simplex*)e2);
		manifold->insert((Simplex*)e3);

		while (!st->empty()){
			Simplex* simplex = st->top();
			st->pop();
			if (simplex->dim == 0){
				/*If it's a vertex, look at all incident edges and see where you can go. There'll be at most one direction*/
				Vertex *vertex = (Vertex*)simplex;
				for (vector<Edge*>::iterator it = vertex->begin(); it != vertex->end(); it++){
					if (this->V->containsPair(vertex, *it)){
						manifold->insert(*it);
						st->push(*it);
						break;
					}
				}
			}
			else if (simplex->dim == 1){
				/*An edge can exit via a vertex or a triangle*/
				Edge *edge = (Edge*)simplex;
				Vertex *vert1 = get<0>(edge->getVertices());
				Vertex *vert2 = get<1>(edge->getVertices());
				if (!this->V->containsPair(vert1, edge)){
					manifold->insert((Simplex*)vert1);
					st->push((Simplex*)vert1);
				}
				if (!this->V->containsPair(vert2, edge)){
					manifold->insert((Simplex*)vert2);
					st->push((Simplex*)vert2);
				}

				for (vector<Triangle*>::iterator it = edge->begin(); it != edge->end(); it++){
					Triangle *tri = *it;
					if (this->V->containsPair(edge, tri)){
						manifold->insert((Simplex*)tri);
						st->push((Simplex*)tri);
						break;
					}
				}
			}
			else{
				/*A triangle can exit out of two of its three edges*/
				Triangle *tri = (Triangle*)simplex;
				Edge *edge1 = get<0>(tri->getEdges());
				Edge *edge2 = get<1>(tri->getEdges());
				Edge *edge3 = get<2>(tri->getEdges());

				Edge *exitEdge1, *exitEdge2;
				if (this->V->containsPair(edge1, tri)){
					exitEdge1 = edge2;
					exitEdge2 = edge3;
				}
				else if (this->V->containsPair(edge2, tri)){
					exitEdge1 = edge1;
					exitEdge2 = edge3;
				}
				else{
					exitEdge1 = edge1;
					exitEdge2 = edge2;
				}

				manifold->insert((Simplex*)exitEdge1);
				manifold->insert((Simplex*)exitEdge2);
				st->push((Simplex*)exitEdge1);
				st->push((Simplex*)exitEdge2);
			}
		}
		delete st;
	}

	return manifold;
}

void Simplicial2Complex::flipAndTranslateVertexFunction(){
	/*Flip the function and find the maximum function value*/
	double max = 0;
	for (int i = 0; i < this->vertexList.size(); i++){
		if (this->vertexList[i]->getFuncValue() > max){
			max = this->vertexList[i]->getFuncValue();
		}
		this->vertexList[i]->setFuncValue(0 - this->vertexList[i]->getFuncValue());
	}
	/*Translate by max*/
	for (int i = 0; i < this->vertexList.size(); i++){
		this->vertexList[i]->setFuncValue(this->vertexList[i]->getFuncValue() + max);
	}
}

//void Simplicial2Complex::outputVertices(string pathname){
//	ofstream output(pathname);
//	output << std::setprecision(10);
//	for (vector<Vertex*>::iterator it = this->vertexList.begin(); it != this->vertexList.end(); it++){
//		Vertex *v = *it;
//		double xCoord = get<0>(v->getCoords());
//		double yCoord = get<1>(v->getCoords());
//		double funcValue = v->getFuncValue();
//		output << xCoord << " " << yCoord << " " << funcValue << endl;
//	}
//	output.close();
//}
//
//void Simplicial2Complex::outputEdges(string pathname){
//	ofstream output(pathname);
//	for (vector<Edge*>::iterator it = this->edgeList.begin(); it != this->edgeList.end(); it++){
//		Edge *e = *it;
//		Vertex *v1 = get<0>(e->getVertices());
//		Vertex *v2 = get<1>(e->getVertices());
//		double xCoord1 = get<0>(v1->getCoords());
//		double yCoord1 = get<1>(v1->getCoords());
//		double xCoord2 = get<0>(v2->getCoords());
//		double yCoord2 = get<1>(v2->getCoords());
//		double funcValue = e->getFuncValue();
//
//		output << xCoord1 << " " << yCoord1 << " " << xCoord2 << " " << yCoord2 << " " << funcValue << endl;
//	}
//	output.close();
//}
//
//void Simplicial2Complex::outputTriangles(string pathname){
//	ofstream output(pathname);
//	for (vector<Triangle*>::iterator it = this->triList.begin(); it != this->triList.end(); it++){
//		Triangle *t = *it;
//		Vertex *v1 = get<0>(t->getVertices());
//		Vertex *v2 = get<1>(t->getVertices());
//		Vertex *v3 = get<2>(t->getVertices());
//		double xCoord1 = get<0>(v1->getCoords());
//		double yCoord1 = get<1>(v1->getCoords());
//		double xCoord2 = get<0>(v2->getCoords());
//		double yCoord2 = get<1>(v2->getCoords());
//		double xCoord3 = get<0>(v3->getCoords());
//		double yCoord3 = get<1>(v3->getCoords());
//		double funcValue = t->getFuncValue();
//
//		output << xCoord1 << " " << yCoord1 << " " << xCoord2 << " " << yCoord2 << " " << xCoord3 << " " << yCoord3 << " " << funcValue << endl;
//	}
//	output.close();
//}
//

Vertex::Vertex(double *coords, double funcValue){
	for (int i = 0; i < DIM; i++) {
		this->coords[i] = coords[i];
	}
	this->funcValue = funcValue;
	this->dim = 0;
}

void Vertex::addEdge(Edge *e){
	this->incidenceList.push_back(e);
}

void Vertex::removeEdge(Edge *e){
	for (vector<Edge*>::iterator i = this->incidenceList.begin(); i != this->incidenceList.end(); i++){
		if (*i == e){
			this->incidenceList.erase(i);
			break;
		}
	}
}

Vertex* Vertex::getAdjacentVertex(Edge *e){
	tuple<Vertex*, Vertex*> vertices = e->getVertices();
	if (get<0>(vertices) == this){
		return get<1>(vertices);
	}
	else{
		return get<0>(vertices);
	}
}

bool Vertex::isAdjacent(Vertex *v){
	for (vector<Edge*>::iterator i = this->incidenceList.begin(); i != this->incidenceList.end(); i++){
		if (this->getAdjacentVertex(*i) == v){
			return true;
		}
	}
	return false;
}

//bool Vertex::isDeleted(){
//	return this->deleted;
//}

//void Vertex::deleteVertex(){
//	this->deleted = true;
//}

double* Vertex::getCoords(){
	return this->coords;
}

double Vertex::getFuncValue(){
	return this->funcValue;
}

int Vertex::degree(){
	return this->incidenceList.size();
}

int Vertex::getVPosition(){
	return this->vPosition;
}

void Vertex::setVposition(int p){
	this->vPosition = p;
}
/* Requires that v is adjacent to this*/
Edge* Vertex::findEdge(Vertex *v){
	for (vector<Edge*>::iterator i = this->incidenceList.begin(); i != this->incidenceList.end(); i++){
		if (this->getAdjacentVertex(*i) == v){
			return *i;
		}
	}
	/*This line should never execute as long as requirements are met*/
	return nullptr;
}

Edge::Edge(Vertex *v1, Vertex *v2){
	this->vertices = make_tuple(v1, v2);
	this->dim = 1;
}

Edge::Edge(Vertex *v1, Vertex *v2, double funcValue){
	this->vertices = make_tuple(v1, v2);
	this->funcValue = funcValue;
	this->dim = 1;
}

tuple<Vertex*, Vertex*> Edge::getVertices(){
	return this->vertices;
}

void Edge::addTriangle(Triangle *t){
	this->incidenceList.push_back(t);
}

void Edge::removeTriangle(Triangle *t){
	for (vector<Triangle*>::iterator i = this->incidenceList.begin(); i != this->incidenceList.end(); i++){
		if (*i == t){
			this->incidenceList.erase(i);
		}
	}
}

Vertex* Edge::oppsiteVertex(Triangle *t){
	tuple<Vertex*, Vertex*, Vertex*> tVertices = t->getVertices();
	tuple<Vertex*, Vertex*> eVertices = this->vertices;
	if (get<0>(tVertices) != get<0>(eVertices) && get<0>(tVertices) != get<1>(eVertices)){
		return get<0>(tVertices);
	}
	else if (get<1>(tVertices) != get<0>(eVertices) && get<1>(tVertices) != get<1>(eVertices)){
		return get<1>(tVertices);
	}
	else{
		return get<2>(tVertices);
	}
}

int Edge::degree(){
	return this->incidenceList.size();
}

//bool Edge::isDeleted(){
//	return this->deleted;
//}

//void Edge::deleteEdge(){
//	this->deleted = true;
//}

double Edge::getFuncValue(){
	return this->funcValue;
}

void Edge::setFuncValue(double funcValue){
	this->funcValue = funcValue;
}

int Edge::getEPosition(){
	return this->ePosition;
}

void Edge::setEposition(int p){
	this->ePosition = p;
}

Vertex* Edge::findVertex(Edge *e){
	if (get<0>(this->vertices) == get<0>(e->vertices) || get<0>(this->vertices) == get<1>(e->vertices)){
		return get<0>(this->vertices);
	}
	else{
		return get<1>(this->vertices);
	}
}

Triangle::Triangle(Edge *e1, Edge* e2, Edge *e3){
	this->edges = make_tuple(e1, e2, e3);
	Vertex *v1 = e1->findVertex(e2);
	Vertex *v2 = e2->findVertex(e3);
	Vertex *v3 = e3->findVertex(e1);
	this->vertices = make_tuple(v1, v2, v3);
	this->dim = 2;
}

Triangle::Triangle(Vertex *v1, Vertex *v2, Vertex *v3){
	this->vertices = make_tuple(v1, v2, v3);
	Edge *e1 = v1->findEdge(v2);
	Edge *e2 = v2->findEdge(v3);
	Edge *e3 = v3->findEdge(v1);
	this->edges = make_tuple(e1, e2, e3);
	this->dim = 2;
}

Triangle::Triangle(Edge *e1, Edge* e2, Edge *e3, double funcValue){
	this->edges = make_tuple(e1, e2, e3);
	Vertex *v1 = e1->findVertex(e2);
	Vertex *v2 = e2->findVertex(e3);
	Vertex *v3 = e3->findVertex(e1);
	this->vertices = make_tuple(v1, v2, v3);
	this->funcValue = funcValue;
	this->dim = 2;
}

Triangle::Triangle(Vertex *v1, Vertex *v2, Vertex *v3, double funcValue){
	this->vertices = make_tuple(v1, v2, v3);
	Edge *e1 = v1->findEdge(v2);
	Edge *e2 = v2->findEdge(v3);
	Edge *e3 = v3->findEdge(v1);
	this->edges = make_tuple(e1, e2, e3);
	this->funcValue = funcValue;
	this->dim = 2;
}

tuple<Edge*, Edge*, Edge*> Triangle::getEdges(){
	return this->edges;
}

tuple<Vertex*, Vertex*, Vertex*> Triangle::getVertices(){
	return this->vertices;
}

//bool Triangle::isDeleted(){
//	return this->deleted;
//}

//void Triangle::deleteTriangle(){
//	this->deleted = true;
//}

double Triangle::getFuncValue(){
	return this->funcValue;
}
void Triangle::setFuncValue(double funcValue){
	this->funcValue = funcValue;
}
int Triangle::getTPosition(){
	return this->tPosition;
}
void Triangle::setTposition(int p){
	this->tPosition = p;
}


DiscreteVField::DiscreteVField(){
	this->VEmap = new unordered_map<Vertex*, Edge*>();
	this->ETmap = new unordered_map<Edge*, Triangle*>();
}

bool DiscreteVField::containsPair(Vertex *v, Edge *e){
	if (this->VEmap->count(v) == 1){
		if (this->VEmap->at(v) == e){
			return true;
		}
	}
	return false;
}

bool DiscreteVField::containsPair(Edge *e, Triangle *t){
	if (this->ETmap->count(e) == 1){
		if (this->ETmap->at(e) == t){
			return true;
		}
	}
	return false;
}

void DiscreteVField::addPair(Vertex *v, Edge *e){
	std::pair<Vertex*, Edge*> pair = std::make_pair(v, e);
	this->VEmap->insert(pair);
}

void DiscreteVField::addPair(Edge *e, Triangle *t){
	std::pair<Edge*, Triangle*> pair = std::make_pair(e, t);
	this->ETmap->insert(pair);
}

void DiscreteVField::removePair(Vertex *v, Edge *e){
	this->VEmap->erase(v);
}

void DiscreteVField::removePair(Edge *e, Triangle *t){
	this->ETmap->erase(e);
}

void DiscreteVField::outputVEmap(){
	ofstream output("VEmap.txt");
	for (unordered_map<Vertex*, Edge*>::iterator it = this->VEmap->begin(); it != this->VEmap->end(); it++){
		std::pair<Vertex*, Edge*> pair = *it;
		output << "Vertex " << pair.first->getVPosition() << " Edge " << pair.second->getEPosition() << " ("<< get<0>(pair.second->getVertices())->getVPosition()<<", "<<get<1>(pair.second->getVertices())->getVPosition()<<")\n";
	}
	output.close();
}

void DiscreteVField::outputETmap(){
	ofstream output("ETmap.txt");
	for (unordered_map<Edge*, Triangle*>::iterator it = this->ETmap->begin(); it != this->ETmap->end(); it++){
		std::pair<Edge*, Triangle*> pair = *it;
		output << "Edge " << pair.first->getFuncValue() << " Triangle " << pair.second->getFuncValue() << "\n";
	}
	output.close();
}

bool Simplex::simplexPointerCompare(const Simplex *s, const Simplex *t){
	/*We sort simplices lexicographically according to the following scheme
	1st		by function value
	2nd		by dimension
	3rd		3 cases
			Vertex
				1st by vPosition
			Edge
				1st by symbolic perturbation
				2nd by ePosition
			Triangle
				1st by first-order symbolic perturbation
				2nd by second-order symbolic perturbation
				3rd by tPosition

	These are each broken down into if-elseif-else statements
	if (<less than>){ return true;}
	else if (<equal to>){ next on list }
	else { return false;}
	*/



	//By function value
	if (s->funcValue < t->funcValue){
		return true;
	}
	else if (s->funcValue == t->funcValue){
		//By dimension
		if (s->dim < t->dim){
			return true;
		}
		else if (s->dim == t->dim){
			//3 cases
			if (s->dim == 0){
				//by vPosition
				return ((Vertex*)s)->getVPosition() < ((Vertex*)t)->getVPosition();
			}
			else if (s->dim == 1){
				//by symbolic perturbation
				if (((Edge*)s)->getSymPerturb() < ((Edge*)t)->getSymPerturb()){
					return true;
				}
				else if (((Edge*)s)->getSymPerturb() == ((Edge*)t)->getSymPerturb()){
					//by ePosition
					return ((Edge*)s)->getEPosition() < ((Edge*)t)->getEPosition();
				}
				else{
					return false;
				}
			}
			else{
				//by first-order symbolic perturbation
				if (get<0>(((Triangle*)s)->getSymPerturb()) < get<0>(((Triangle*)t)->getSymPerturb())){
					return true;
				}
				else if (get<0>(((Triangle*)s)->getSymPerturb()) == get<0>(((Triangle*)t)->getSymPerturb())){
					//by second-order symbolic perturbation
					if (get<1>(((Triangle*)s)->getSymPerturb()) < get<1>(((Triangle*)t)->getSymPerturb())){
						return true;
					}
					else if (get<1>(((Triangle*)s)->getSymPerturb()) == get<1>(((Triangle*)t)->getSymPerturb())){
						//by tPosition
						return ((Triangle*)s)->getTPosition() < ((Triangle*)t)->getTPosition();
					}
					else{
						return false;
					}
				}
				else{
					return false;
				}
			}
		}
		else{
			return false;
		}
	}
	else{
		return false;
	}
}
