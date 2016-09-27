#include <vector>
#include <tuple>
#include <unordered_set>
#include <set>
#include <stack>
#include <unordered_map>
#include <map>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>

#define MAX_DIM 3
#define EPS_compare 1e-8

using namespace std;

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
	
	// This is original - deprecated
	static bool simplexPointerCompare(const Simplex* s, const Simplex* t);
	static bool simplexPointerCompare2(const Simplex *s, const Simplex *t);
	unsigned int filtrationPosition;
	Simplex* prev = NULL;
};

class Simplicial2Complex{
	vector<Vertex*> vertexList;
	vector<Edge*> edgeList;
	vector<Triangle*> triList;
	unordered_set<Simplex*> criticalSet;
	DiscreteVField *V;

public:
	Simplicial2Complex();

	int addVertex(Vertex *v);
	int addEdge(Edge *e);
	int addTriangle(Triangle* t);
	// Require, no duplicate edge
	// no duplicate triangles
	Vertex* getVertex(int position);
	Edge* getEdge(int position);
	Triangle* getTriangle(int position);
	void addCriticalPoint(Simplex *s);
	void removeCriticalPoint(Simplex *s);
	bool isCritical(Simplex *s);
	void setDiscreteVField(DiscreteVField *V);
	int order();
	void buildRipsComplex(double radius, double eps);
	void outputArcs(string, string);
	void buildComplexFromFile(string pathname);
	void buildComplexFromFile2(string pathname);
	void buildComplexFromFile2_BIN(string pathname);
	void outputComplex(string pathname);
	void buildPsuedoMorseFunction();
	/*Requires: s is a critical simplex*/
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
	void sortVertices(){
		sort(vertexList.begin(), vertexList.end(), Simplex::simplexPointerCompare2);
	}
};

class Vertex:public Simplex{
	double coords[MAX_DIM];
	vector<Edge*> incidenceList;
	int vPosition;
	int oriPosition;
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
	double* getCoords();
	double getFuncValue();
	void setFuncValue(double funcValue){
		this->funcValue = funcValue;
	}
	int degree();
	int getVPosition();
	void setVposition(int p);
	int getoriPosition(){
		return oriPosition;
	}
	void setoriposition(int p){
		oriPosition = p;
	}
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

	void output(){
	    cout<< this->funcValue << "\t";
	}
};

class Edge:public Simplex{
	tuple<Vertex*, Vertex*> vertices;
	vector<Triangle*> incidenceList;
	int ePosition;


public:
	double persistence;
	int critical_type;
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
	double getFuncValue();
	void setFuncValue(double funcValue);
	int getEPosition();
	void setEposition(int p);
	/*requires that this and e share a Vertex*/
	Vertex* findVertex(Edge* e);

	vector<Triangle*>::iterator begin(){
		return this->incidenceList.begin();
	}
	vector<Triangle*>::iterator end(){
		return this->incidenceList.end();
	}
	
	void output(){
	    cout<< get<0>(vertices)->funcValue << " " << get<1>(vertices)->funcValue << "\t";
	}
};

class Triangle:public Simplex{
	tuple<Edge*, Edge*, Edge*> edges;
	tuple<Vertex*, Vertex*, Vertex*> vertices;
	int tPosition;

public:
	Triangle(Edge *e1, Edge *e2, Edge *e3);
	Triangle(Vertex *v1, Vertex *v2, Vertex *v3);
	Triangle(Edge *e1, Edge *e2, Edge *e3, double funcValue);
	Triangle(Vertex *v1, Vertex *v2, Vertex *v3, double funcValue);
	tuple<Edge*, Edge*, Edge*> getEdges();
	tuple<Vertex*, Vertex*, Vertex*> getVertices();
	double getFuncValue();
	void setFuncValue(double funcValue);
	int getTPosition();
	void setTposition(int p);
	
	void output(){
	    cout<< get<0>(vertices)->funcValue << " " << get<1>(vertices)->funcValue << " " << get<2>(vertices)->funcValue << "\t";
	}
	
	void SortVertex(){
		Vertex *v1, *v2, *v3;
		v1 = get<0>(vertices); v2 = get<1>(vertices); v3 = get<2>(vertices); 
		if (simplexPointerCompare2(v1, v2)) swap(v1, v2);
		if (simplexPointerCompare2(v1, v3)) swap(v1, v3);
		if (simplexPointerCompare2(v2, v3)) swap(v2, v3);
		vertices = make_tuple(v1, v2, v3);
		// the first is the one with highest density
	}
};

class DiscreteVField{
private:
	unordered_map<Vertex*, Edge*> *VEmap;
	unordered_map<Edge*, Triangle*> *ETmap;
	unordered_set<Edge*> *Eflag;
	unordered_map<Simplex*, vector<Simplex*>* > *jumpmap;

public:
	DiscreteVField();
	Edge* containsPair(Vertex *v);
	Triangle* containsPair(Edge *e);
	bool containsEdge(Edge *e);
	void addPair(Vertex *v, Edge *e);
	void addPair(Edge *e, Triangle *t);
	void removePair(Vertex *v, Edge *e);
	void removePair(Edge *e, Triangle *t);
	void outputVEmap();
	void outputETmap();
	bool CanJump(Simplex*);
	void UnWarp(Simplex* s, bool flip, vector<Simplex*>* path);
	void addJump(Simplex*, vector<Simplex*>*);
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

bool Simplicial2Complex::isCritical(Simplex *s){
	return (this->criticalSet.count(s) == 1);
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
	cout<< "Writing 1-stable manifold\n";
	ofstream output_info("output_info.txt", ios_base::out | ios_base::trunc);
	int counter = 0;
	for(unordered_set<Simplex*>::iterator it = this->cBegin(); it != this->cEnd(); it++){
		Simplex *s = *it;
		// && s->critical_type == 1 ---> try this?
		if(s->dim == 1){
			Edge *e = (Edge*)s;
			// Some edge not critical are inserted?!
			if (e->critical_type == 2 && e->persistence < delta - EPS_compare) continue;
			output_info << e->persistence << " " << e->getEPosition() << " "
						<< get<0>(e->getVertices())->funcValue << "\n";
			set<Simplex*> *manifold = this->descendingManifold((Simplex*)e);
			for (set<Simplex*>::iterator it2 = manifold->begin(); it2 != manifold->end(); it2++){
				manifolds.insert(*it2);
			}
			delete manifold;
			counter++;
		}
	}
	output_info.close();
	cout << "Written " << counter << "arcs\n";
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
	cout << "\tReading " << numOfVertices << "vertices" << endl;
	for (int i = 0; i < numOfVertices; i++) {
		double coords[MAX_DIM];
		double funcValue;
		for (int j = 0; j < DIM; j++) {
			file >> coords[j];
		}
		file >> funcValue;
		Vertex *v = new Vertex(coords, funcValue);
		int vPosition = this->addVertex(v);
		v->setVposition(vPosition);
	}
	cout << "\tDone." << endl;

	int numOfEdges;
	file >> numOfEdges;
	cout << "\tReading " << numOfEdges << "edges" << endl;
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
	cout << "\tDone." << endl;

	int numOfTris;
	file >> numOfTris;
	cout << "\tReading " << numOfTris << "triangles" << endl;
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
	sortVertices();
	cout << "\tDone." << endl;

	file.close();
}

void Simplicial2Complex::buildComplexFromFile2(string pathname) {
	ifstream file(pathname);
	int numOfVertices;
	file >> numOfVertices;
	cout << "\tReading " << numOfVertices << "vertices" << endl;

	for (int i = 0; i < numOfVertices; i++) {
		double coords[MAX_DIM];
		double funcValue;
		for (int j = 0; j < DIM; j++) {
			file >> coords[j];
		}
		file >> funcValue;
		Vertex *v = new Vertex(coords, funcValue);
		int vPosition = this->addVertex(v);
		v->setVposition(vPosition);
	}
	cout << "\tDone." << endl;

	int numOfEdges;
	file >> numOfEdges;
	cout << "\tReading " << numOfEdges << "edges" << endl;
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
	cout << "\tDone." << endl;

	int numOfTris;
	file >> numOfTris;
	cout << "\tReading " << numOfTris << "triangles" << endl;
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
	cout << "\tDone." << endl;
	file.close();
}


void Simplicial2Complex::buildComplexFromFile2_BIN(string pathname) {
	ifstream file(pathname, ios::binary);
	char* int_buffer = new char[sizeof(int)];
	int* int_reader = (int*) int_buffer;
	char* double_buffer = new char[sizeof(double)];
	double* double_reader = (double*) double_buffer;
	
	
	int numOfVertices;	
	file.read(int_buffer, sizeof(int));
	numOfVertices = *int_reader;	
	cout << "\tReading " << numOfVertices << "vertices" << endl;

	for (int i = 0; i < numOfVertices; i++) {
		double coords[MAX_DIM];
		double funcValue;
		for (int j = 0; j < DIM; j++) {
			file.read(double_buffer, sizeof(double));
			coords[j] = *double_reader;
		}
		file.read(double_buffer, sizeof(double));
		funcValue = *double_reader;
		funcValue = (int)(funcValue*1e5)/1.0e5;
		Vertex *v = new Vertex(coords, funcValue);
		int vPosition = this->addVertex(v);
		// this will be over written by sorted order later
		v->setVposition(vPosition);
		v->setoriposition(vPosition);
	}
	
	
	// Use flipped function --- maxma -> minima
	// So we can look at vertex-edge pair
	flipAndTranslateVertexFunction();
	
	
	cout << "\tSorting " << numOfVertices << "vertices" << endl;

	vector<Vertex*> sortedvert(vertexList);
	sort(sortedvert.begin(), sortedvert.end(), Simplex::simplexPointerCompare2);
	int counter = 0;
	for (auto i = sortedvert.begin(); i < sortedvert.end(); ++i){
		(*i)->setVposition(counter);
		counter++;
	}

	
	cout << "\tDone." << endl;

	int numOfEdges;
	file.read(int_buffer, sizeof(int));
	numOfEdges = *int_reader;
	cout << "\tReading " << numOfEdges << "edges" << endl;
	for (int i = 0; i < numOfEdges; i++) {
		int vIndex1, vIndex2;
		file.read(int_buffer, sizeof(int));
		vIndex1 = *int_reader;
		file.read(int_buffer, sizeof(int));
		vIndex2 = *int_reader;
		Vertex *v1 = this->getVertex(vIndex1);
		Vertex *v2 = this->getVertex(vIndex2);
		Edge *e = new Edge(v1, v2);
		v1->addEdge(e);
		v2->addEdge(e);
		int ePosition = this->addEdge(e);
		e->setEposition(ePosition);
	}
	cout << "\tDone." << endl;

	int numOfTris;
	file.read(int_buffer, sizeof(int));
	numOfTris = *int_reader;
	cout << "\tReading " << numOfTris << "triangles" << endl;
	for (int i = 0; i < numOfTris; i++) {
		int vIndex1, vIndex2, vIndex3;
		file.read(int_buffer, sizeof(int));
		vIndex1 = *int_reader;
		file.read(int_buffer, sizeof(int));
		vIndex2 = *int_reader;
		file.read(int_buffer, sizeof(int));
		vIndex3 = *int_reader;
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
	vertexList = sortedvert;
	cout << "\tDone." << endl;
	file.close();
}


void Simplicial2Complex::buildPsuedoMorseFunction(){
	cout << "\t Processing "<< edgeList.size() <<" edges\n";
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
		//e->output();
		//cout<<"E: "<< max->getFuncValue() << " " << total << endl;
	}

	cout << "\t Processing "<< triList.size() <<" triangles\n";
	for (unsigned int i = 0; i < this->triList.size(); i++){
		Triangle *t = this->triList.at(i);
		tuple<Edge*, Edge*, Edge*> edges = t->getEdges();
		Edge *max = get<0>(edges);
		Edge *e2 = get<1>(edges);
		Edge *e3 = get<2>(edges);
		tuple<Vertex*, Vertex*, Vertex*> vertices = t->getVertices();
		double total = get<0>(vertices)->getFuncValue() + get<1>(vertices)->getFuncValue() + get<2>(vertices)->getFuncValue();
		if (e2->getFuncValue() > max->getFuncValue()){
			max = e2;
		}
		if (e3->getFuncValue() > max->getFuncValue()){
			max = e3;
		}
		t->setFuncValue(max->getFuncValue());
		//t->output();
		//cout<<"T: "<<max->getFuncValue()<<" "<<max->getSymPerturb()<<" "<< total<<endl;
	}
}


set<Simplex*>* Simplicial2Complex::descendingManifold(Simplex* s){
	// Discrete Vector Field: V exist here
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
				Vertex *vert = (Vertex*)simplex;
				Edge* paired_edge = V->containsPair(vert);
				/*If it's a vertex, look at all incident edges and see where you can go. There'll be at most one direction*/
				if (paired_edge!=NULL){
					manifold->insert((Simplex*) paired_edge);
					st->push((Simplex*) paired_edge);
				}
			}
			else{
				/*An edge can only be entered from a vertex, so that uses up its arrow, and it can only leave via the other vertex*/
				/*if (V->CanJump(simplex)){
					vector<Simplex*>* path = new vector<Simplex*>();
					V->UnWarp(simplex, true, path);
					// insert all simplex to manifold;
					for(int s_index = 0; s_index < (path->size()); s_index++){
						manifold->insert((*path)[s_index]);
					}
				}*/
				Edge *edge = (Edge*)simplex;
				Vertex *vert1 = get<0>(edge->getVertices());
				Vertex *vert2 = get<1>(edge->getVertices());
				Vertex *exitVert = vert1;
				if (!(this->V->containsPair(vert2) == edge)){
					exitVert = vert2;
				}
				manifold->insert((Simplex*)exitVert);
				st->push((Simplex*)exitVert);
			}

		}
		delete st;
	}
	else{
		cout<<"we shouldn't be here now, SOMETHIING IS WRONG"<<endl;
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
					if (this->V->containsPair(vertex) == *it){
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
				if (!(this->V->containsPair(vert1) == edge)){
					manifold->insert((Simplex*)vert1);
					st->push((Simplex*)vert1);
				}
				if (!(this->V->containsPair(vert2) == edge)){
					manifold->insert((Simplex*)vert2);
					st->push((Simplex*)vert2);
				}

				for (vector<Triangle*>::iterator it = edge->begin(); it != edge->end(); it++){
					Triangle *tri = *it;
					if (this->V->containsPair(edge) == tri){
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
				if (this->V->containsPair(edge1) == tri){
					exitEdge1 = edge2;
					exitEdge2 = edge3;
				}
				else if (this->V->containsPair(edge2) == tri){
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
	if (simplexPointerCompare2(v2, v1))
		this->vertices = make_tuple(v1, v2);
	else
		this->vertices = make_tuple(v2, v1);
	this->dim = 1;
}

Edge::Edge(Vertex *v1, Vertex *v2, double funcValue){
	if (simplexPointerCompare2(v2, v1))
		this->vertices = make_tuple(v1, v2);
	else
		this->vertices = make_tuple(v2, v1);
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
	SortVertex();
}

Triangle::Triangle(Vertex *v1, Vertex *v2, Vertex *v3){
	this->vertices = make_tuple(v1, v2, v3);
	Edge *e1 = v1->findEdge(v2);
	Edge *e2 = v2->findEdge(v3);
	Edge *e3 = v3->findEdge(v1);
	this->edges = make_tuple(e1, e2, e3);
	this->dim = 2;
	SortVertex();
}

Triangle::Triangle(Edge *e1, Edge* e2, Edge *e3, double funcValue){
	this->edges = make_tuple(e1, e2, e3);
	Vertex *v1 = e1->findVertex(e2);
	Vertex *v2 = e2->findVertex(e3);
	Vertex *v3 = e3->findVertex(e1);
	this->vertices = make_tuple(v1, v2, v3);
	this->funcValue = funcValue;
	this->dim = 2;
	SortVertex();
}

Triangle::Triangle(Vertex *v1, Vertex *v2, Vertex *v3, double funcValue){
	this->vertices = make_tuple(v1, v2, v3);
	Edge *e1 = v1->findEdge(v2);
	Edge *e2 = v2->findEdge(v3);
	Edge *e3 = v3->findEdge(v1);
	this->edges = make_tuple(e1, e2, e3);
	this->funcValue = funcValue;
	this->dim = 2;
	SortVertex();
}

tuple<Edge*, Edge*, Edge*> Triangle::getEdges(){
	return this->edges;
}

tuple<Vertex*, Vertex*, Vertex*> Triangle::getVertices(){
	return this->vertices;
}


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
	this->Eflag = new unordered_set<Edge*>();
	this->jumpmap = new unordered_map<Simplex*, vector<Simplex*>*>();
}

bool DiscreteVField::CanJump(Simplex* s){
	if (this->jumpmap->count(s) == 1) return true;
	else return false;
}

void DiscreteVField::UnWarp(Simplex* s, bool flip, vector<Simplex*>* path){
	vector<Simplex*>* compressed_path = jumpmap->at(s);
	if (flip){
		for(vector<Simplex*>::reverse_iterator it = compressed_path->rbegin(); it != compressed_path->rend(); ++it){
			if (*it == s){
				// this should happen at the end.
				path->push_back(*it);
				continue;
			}
			if (this->CanJump(*it)){
				// order might not matter
				this->UnWarp(*it, !flip, path);
			}else{
				path->push_back(*it);
			}
		}
	}else{
		for(vector<Simplex*>::iterator it = compressed_path->begin(); it != compressed_path->end(); ++it){
			if (*it == s){
				// this should happen at the beginning.
				path->push_back(*it);
				continue;
			}
			if (this->CanJump(*it)){
				this->UnWarp(*it, !flip, path);
			}else{
				path->push_back(*it);
			}
		}
	}
}
void DiscreteVField::addJump(Simplex* s, vector<Simplex*>* jp){
	std::pair<Simplex*, vector<Simplex*>*> pair = std::make_pair(s, jp);
	this->jumpmap->insert(pair);
}

bool DiscreteVField::containsEdge(Edge *e){
	if (this->Eflag->count(e) == 1)
		return true;
	else
		return false;
}

Edge* DiscreteVField::containsPair(Vertex *v){
	if (this->VEmap->count(v) == 1){
		return this->VEmap->at(v);
	}
	return NULL;
}

Triangle* DiscreteVField::containsPair(Edge *e){
	if (this->ETmap->count(e) == 1){
			return this->ETmap->at(e);
	}
	return NULL;
}

void DiscreteVField::addPair(Vertex *v, Edge *e){
	std::pair<Vertex*, Edge*> pair = std::make_pair(v, e);
	this->VEmap->insert(pair);
	this->Eflag->insert(e);
}

void DiscreteVField::addPair(Edge *e, Triangle *t){
	std::pair<Edge*, Triangle*> pair = std::make_pair(e, t);
	this->ETmap->insert(pair);
}

void DiscreteVField::removePair(Vertex *v, Edge *e){
	if (this->containsPair(v) != e){
		cerr << "Removing unexisting pair" <<endl;
	}
	this->VEmap->erase(v);
	this->Eflag->erase(e);
}

void DiscreteVField::removePair(Edge *e, Triangle *t){
	if (this->containsPair(e) != t){
		cerr << "Removing unexisting pair" <<endl;
	}
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


bool Simplex::simplexPointerCompare2(const Simplex *s, const Simplex *t){
	//By function value
	double f1 = s->funcValue, f2 = t->funcValue;
	if (f1 < f2 - EPS_compare){
		return true;
	}else if(f1 > f2 + EPS_compare){
		return false;
	}
	else{
		//By dimension
		short int d1 = s->dim, d2 = t->dim;
		if (d1 < d2){
			return true;
		}else if(d1 > d2){
			return false;
		}
		else{
			//3 cases
			if (d1 == 0){
				//by vPosition
				return ((Vertex*)s)->getVPosition() < ((Vertex*)t)->getVPosition();
			}
			else if (d1 == 1){
				//by symbolic perturbation
				Edge* e1 = (Edge*)s;
				Edge* e2 = (Edge*)t;
				tuple<Vertex*, Vertex*> e1v = e1->getVertices();
				tuple<Vertex*, Vertex*> e2v = e2->getVertices();
				
				//for(int i = 0; i < 2; ++i){
				if (get<0>(e1v) != get<0>(e2v))
					return simplexPointerCompare2(get<0>(e1v), get<0>(e2v));
				if (get<1>(e1v) != get<1>(e2v))
					return simplexPointerCompare2(get<1>(e1v), get<1>(e2v));
				//}
				cout << "Caught edge with same set of vertices\n";
				return false;
			}
			else{
				//by first-order symbolic perturbation
				Triangle* t1 = (Triangle*) s;
				Triangle* t2 = (Triangle*) t;
				
				tuple<Vertex*, Vertex*, Vertex*> t1v = t1->getVertices();
				tuple<Vertex*, Vertex*, Vertex*> t2v = t2->getVertices();
				

				//by tPosition
				//for(int i = 0; i < 3; ++i){
				if (get<0>(t1v) != get<0>(t2v))
					return simplexPointerCompare2(get<0>(t1v), get<0>(t2v));
				if (get<1>(t1v) != get<1>(t2v))
					return simplexPointerCompare2(get<1>(t1v), get<1>(t2v));
				if (get<2>(t1v) != get<2>(t2v))
					return simplexPointerCompare2(get<2>(t1v), get<2>(t2v));
				//}
				cout << "Caught triangle with same set of vertices\n";
				return false;

			}
		}
	}
}
