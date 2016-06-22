#include "Simplicial2Complex.h"
#include "sparseZ2matrix.h"
#include <stack>
//#include <windows.h>
//#include <chrono>
#include <iostream>
using namespace std;

#define DEBUG 1

typedef struct {
	Vertex *min;
	Edge *saddle;
	double persistence;
	double symPerturb;
}persistencePair01;
typedef struct {
	Edge *saddle;
	Triangle *max;
	double persistence;
	double symPerturb1;
	double symPerturb2;
}persistencePair12;

class PersistencePairs{
	Simplicial2Complex *K;
	/*ms: min-saddle or 0-1
	  sm: saddle-max or 1-2*/
	vector<persistencePair01> msPersistencePairs;
	vector<persistencePair12> smPersistencePairs;

public:
	vector<Simplex*> filtration;
	static bool persistencePairCompare01(const persistencePair01& p, const persistencePair01& q){
		if(p.persistence < q.persistence){
			return true;
		}
		else if (p.persistence == q.persistence && p.symPerturb < q.symPerturb){
			return true;
		}else{
			return false;
		}
	}
	static bool persistencePairCompare(const persistencePair01& p, const persistencePair12& q){
		return p.persistence <= q.persistence;
	}
	static bool persistencePairCompare12(const persistencePair12& p, const persistencePair12& q){
		if(p.persistence < q.persistence){
			return true;
		} else if (p.persistence == q.persistence && p.symPerturb1 < q.symPerturb1){
			return true;
		} else if ( p.persistence == q.persistence && p.symPerturb1 == q.symPerturb1 && p.symPerturb2 < q.symPerturb2){
			return true;
		} else {
			return false;
		}
	}
	static bool persistencePairCompare(const persistencePair12& p, const persistencePair01& q){
		return p.persistence < q.persistence;
	}

	PersistencePairs(Simplicial2Complex *K){
		this->K = K;
	}
	void buildFiltration();
	
	void computePersistencePairs2();
	void computePersistencePairsWithClear();
	void computePersistencePairsSeparate();

	vector<Simplex*>* isCancellable(const persistencePair01&, ofstream&);
	vector <Simplex*>* isCancellable(const persistencePair12&, ofstream&);
	/*Requires that VPath was returned by isCancellable, and no calls to 
	cancelAlongPath have been made between then and this call*/
	void cancelAlongVPath(vector<Simplex*>* VPath);

	void cancelPersistencePairs(double delta);

	void outputPersistencePairs(string pathname){
		ofstream output(pathname, ios_base::trunc | ios_base::out);
		for (int i = 0; i < this->msPersistencePairs.size(); i++){
			output << "Vertex Edge " << this->msPersistencePairs[i].persistence<<"\n";
		}
		for (int i = 0; i < this->smPersistencePairs.size(); i++){
			output << "Edge Triangle " << this->smPersistencePairs[i].persistence << "\n";
		}
	}
};

void PersistencePairs::buildFiltration(){
	K = this->K;
	/*Reserve 3 times the vertices as there are in K. This is roughly how many simplices are usually in K, so filtration doesn't need to resize as often*/
	cout << "\t inserting simplicies\n";
	filtration.reserve(K->order() * 3);
	this->filtration.insert(filtration.end(), K->vBegin(), K->vEnd());
	this->filtration.insert(filtration.end(), K->eBegin(), K->eEnd());
	this->filtration.insert(filtration.end(), K->tBegin(), K->tEnd());

	/*Sort the simplices according to their function value (including symbolic perturbations)*/
	cout << "\t sorting simplicies...\n";
	sort(filtration.begin(), filtration.end(), Simplex::simplexPointerCompare);
	for (unsigned int i = 0; i < this->filtration.size(); i++){
		Simplex *s = filtration[i];
		s->filtrationPosition = i;
	}
} 

/*This computes (0,1)-persistence pairs including the ones with 0 persistence which are canceled out
in Simplex::cancel0PersistencePairs(). This is not a problem since they would be the first to be canceled anyway,
and we need to search paths between critical simplices to test for cancelability. By treating every simplex as critical,
we lose the ability to know right away whether a given pair has a path between them, but this is okay because we need
to search for multiple paths when we cancel them anyway.*/

//void PersistencePairs::compute01PersistencePairs(){
//	UnionFind uf;
//	for (vector<Simplex*>::iterator it = filtration.begin(); it != filtration.end(); it++){
//		Simplex *s = *it;
//		if (s->dim == 1){
//			Vertex *v = (Vertex*)s;
//			/*Every vertex creates a new connected component, so it is tagged positive*/
//			v->parity = Simplex::POSITIVE;
//			uf.makeSet(v);
//		}
//		else if (s->dim == 2){
//			Edge *e = (Edge*)s;
//			Vertex *v1 = get<0>(e->getVertices());
//			Vertex *v2 = get<1>(e->getVertices());
//			Vertex *rootV1 = uf.findSet(v1);
//			Vertex *rootV2 = uf.findSet(v2);
//
//			if (rootV1 != rootV2)
//			{
//				/*If rootV1 != rootV2, then adding e connects different connected components, so e is negative*/
//				e->parity = Simplex::NEGATIVE;
//				/*We need to pair e with one of the roots, namely the one that came into the filtration first. The "older" one*/
//				Vertex *min = rootV1;
//				if (rootV2->getFuncValue() < min->getFuncValue()){
//					min = rootV2;
//				}
//				persistencePair01 pair = { min, e, e->getFuncValue() - min->getFuncValue() };
//				/*If the persistence of pair is 0, then we don't bother adding it to the vector of (0,1)-persistence pairs, since we already canceled them*/
//				if (pair.persistence != 0){
//					this->msPersistencePairs.push_back(pair);
//				}
//			}
//
//
//		}
//	}
//}

void PersistencePairs::computePersistencePairs2(){
	SparseZ2Matrix *boundaryMatrix = new SparseZ2Matrix(this->filtration.size(), this->filtration.size());

	cout << "\tInitializing boundary matrix...\n";
	cout << "\t\t" << this->filtration.size() << " x " << this->filtration.size() << "\n";
	for (unsigned int i = 0; i < this->filtration.size(); i++) {
		Simplex *s = filtration[i];
		if (s->dim == 0) {
			//no boundary. no bits to set in matrix
		}
		else if (s->dim == 1) {
			Edge *e = (Edge*)s;
			Vertex *v1 = get<0>(e->getVertices());
			Vertex *v2 = get<1>(e->getVertices());
			if (v1->filtrationPosition < v2->filtrationPosition) {
				boundaryMatrix->set(v1->filtrationPosition, i);
				boundaryMatrix->set(v2->filtrationPosition, i);
			}
		}
		else if (s->dim == 2) {
			Triangle *t = (Triangle*)s;
			Edge *e1 = get<0>(t->getEdges());
			Edge *e2 = get<1>(t->getEdges());
			Edge *e3 = get<2>(t->getEdges());
			Edge *edges[3] = { e1,e2,e3 };
			Edge *min = e1;
			for (int j = 0; j < 3; j++)
			{
				int min = j;
				for (int k = j; k < 3; k++) {
					if (edges[k]->filtrationPosition < edges[min]->filtrationPosition) {
						min = k;
					}
				}
				Edge *tmp = edges[min];
				edges[min] = edges[j];
				edges[j] = tmp;
			}
			boundaryMatrix->set(edges[0]->filtrationPosition, i);
			boundaryMatrix->set(edges[1]->filtrationPosition, i);
			boundaryMatrix->set(edges[2]->filtrationPosition, i);
		}
	}

	boundaryMatrix->output("boundary.txt");

	cout << "\tInitializing pivot array...\n";
	/*The ith entry of pivots stores the number j such that the ith row is the pivot row of the jth column in the reduced matrix*/
	int *pivots = new int[this->filtration.size()];
	/*Initially, all entries are set to -1, meaning "unknown"*/
	for (unsigned int i = 0; i < this->filtration.size(); i++){
		pivots[i] = -1;
	}

	/*Go through the columns in filtration order*/
	for (unsigned int j = 0; j < this->filtration.size(); j++){
		if (j % 100000 == 0) cout << "\tHandling column " << j <<"/"<< this->filtration.size() << "...\n";
		/*If the column is 0, then there is nothing to do. Move on to the next column*/
		if (boundaryMatrix->isZeroColumn(j) == true) continue;
		/*If the column is not 0, find the pivot*/
		unsigned int pivotJ = boundaryMatrix->getPivotRow(j);
		/*As long as column j is not 0, and there is some previous column with the same pivot row,
		add the previous column to column j in order to cancel the pivot. Then find the new pivot 
		and repeat until column j is 0 or there is no previous column to cancel its pivot*/
		while (boundaryMatrix->isZeroColumn(j) == false && pivots[pivotJ] != -1){
			boundaryMatrix->add((unsigned int)pivots[pivotJ], j);
			if (boundaryMatrix->isZeroColumn(j) == false){
				pivotJ = boundaryMatrix->getPivotRow(j);
			}
		}

		/*If the column j gets to this point and is still not 0, then it must have a pivot that can't be canceled.
		In which case, the simplex corresponding to the pivot row and the the simplex corresponding to column j form
		a persistence pair*/
		if (boundaryMatrix->isZeroColumn(j) == false){
			/*Set column j as the column with row pivotJ as its pivot in the reduced matrix*/
			unsigned int i = pivotJ;
			pivots[i] = j;

			/*Add persistence pair*/
			Simplex *s1 = this->filtration[i];
			Simplex *s2 = this->filtration[j];
			if (s1->dim == 0){
				Vertex *v = (Vertex*)s1;
				Edge *e = (Edge*)s2;
				double persistence = e->funcValue - v->funcValue;

				persistencePair01 pp = { v, e, persistence };
				this->msPersistencePairs.push_back(pp);
			}
			else{
				Edge* e = (Edge*)s1;
				Triangle* t = (Triangle*)s2;
				double persistence = t->funcValue - e->funcValue;

				persistencePair12 pp = { e, t, persistence };
				this->smPersistencePairs.push_back(pp);
			}
		}
	}
	/*for (vector<Vertex*>::iterator it = this->K->vBegin(); it != this->K->vEnd(); it++) {
	K->addCriticalPoint((Simplex*)*it);
	}
	for (vector<Edge*>::iterator it = this->K->eBegin(); it != this->K->eEnd(); it++) {
	K->addCriticalPoint((Simplex*)*it);
	}
	for (vector<Triangle*>::iterator it = this->K->tBegin(); it != this->K->tEnd(); it++) {
	K->addCriticalPoint((Simplex*)*it);
	}*/

	delete boundaryMatrix;
	for (vector<persistencePair01>::iterator it = this->msPersistencePairs.begin(); it != this->msPersistencePairs.end(); it++) {
		persistencePair01 p = *it;
		this->K->addCriticalPoint(p.min);
		this->K->addCriticalPoint(p.saddle);
	}
	cout << "# of min-saddle pairs: " << this->msPersistencePairs.size() << "\n";
	cout << "# of saddle-max pairs: " << this->smPersistencePairs.size() << "\n";
	int count = 0;
	for (vector<Edge*>::iterator it = K->eBegin(); it != K->eEnd(); it++) {
		count++;
	}
	cout << "# of edges in complex: " << count << "\n";
}

//Compute persistence pairs
void PersistencePairs::computePersistencePairsWithClear(){
	SparseZ2Matrix *boundaryMatrix = new SparseZ2Matrix(this->filtration.size(), this->filtration.size());

	cout << "\tInitializing boundary matrix...\n";
	cout << "\t\t" << this->filtration.size() << " x " << this->filtration.size() << "\n";
	for (unsigned int i = 0; i < this->filtration.size(); i++){
		Simplex *s = filtration[i];
		if (s->dim == 0){
			//no boundary. no bits to set in matrix
		}
		else if (s->dim == 1){
			Edge *e = (Edge*)s;
			Vertex *v1 = get<0>(e->getVertices());
			Vertex *v2 = get<1>(e->getVertices());
			if (v1->filtrationPosition < v2->filtrationPosition) {
				boundaryMatrix->set(v1->filtrationPosition, i);
				boundaryMatrix->set(v2->filtrationPosition, i);
			}else{
				boundaryMatrix->set(v2->filtrationPosition, i);
				boundaryMatrix->set(v1->filtrationPosition, i);
			}
		}
		else if (s->dim == 2){
			Triangle *t = (Triangle*)s;
			Edge *e1 = get<0>(t->getEdges());
			Edge *e2 = get<1>(t->getEdges());
			Edge *e3 = get<2>(t->getEdges());
			Edge *edges[3] = { e1,e2,e3 };
			for(int j = 0; j < 3; j++){
				for(int k = j; k < 3; k++){
					if(edges[k]->filtrationPosition < edges[j]->filtrationPosition){
						Edge *tmp = edges[k];
						edges[k] = edges[j];
						edges[j] = tmp;
					}
				}
			}
			boundaryMatrix->set(edges[0]->filtrationPosition, i);
			boundaryMatrix->set(edges[1]->filtrationPosition, i);
			boundaryMatrix->set(edges[2]->filtrationPosition, i);
		}
	}

	boundaryMatrix->output("boundary.txt");

	cout << "\tInitializing pivot array...\n";
	/*The ith entry of pivots stores the number j such that the ith row is the pivot row of the jth column in the reduced matrix*/
	int *pivots = new int[this->filtration.size()];
	/*Initially, all entries are set to -1, meaning "unknown"*/
	for (unsigned int i = 0; i < this->filtration.size(); i++){
		pivots[i] = -1;
	}

	double average1 = 0, average2 = 0, average3 = 0;
	int a1Count = 0, a2Count = 0, a3Count = 0;


	for (int dim = 2; dim >= 0; dim--)
	{
		cout << "\tAt dimension " << dim << endl;
		/*Go through the columns in filtration order*/
		for (unsigned int j = 0; j < this->filtration.size(); j++){

			if (j % 100000 == 0) {
				cout << "\r";
				cout << "\tHandling column " << j << "/" << this->filtration.size() << "...";
			}

			/*If the column is 0, then there is nothing to do. Move on to the next column*/
			//chrono::time_point<chrono::system_clock> start, end;
			//start = chrono::system_clock::now();
			if (boundaryMatrix->isZeroColumn(j) == true) continue;
			//end = chrono::system_clock::now();
			//chrono::duration<double> elapsedDuration = end - start;
			//double elapsedSeconds = elapsedDuration.count();
			//average1 = (elapsedSeconds + a1Count * average1) / (a1Count + 1);
			//a1Count++;



			if (this->filtration[j]->dim == dim){
				/*If the column is not 0, find the pivot*/
				unsigned int pivotJ = boundaryMatrix->getPivotRow(j);

				/*As long as column j is not 0, and there is some previous column with the same pivot row,
				add the previous column to column j in order to cancel the pivot. Then find the new pivot
				and repeat until column j is 0 or there is no previous column to cancel its pivot*/
				while (boundaryMatrix->isZeroColumn(j) == false && pivots[pivotJ] != -1){
					//start = chrono::system_clock::now();
					boundaryMatrix->add((unsigned int)pivots[pivotJ], j);
					//end = chrono::system_clock::now();
					//elapsedDuration = end - start;
					//elapsedSeconds = elapsedDuration.count();
					//average2 = (elapsedSeconds + a2Count * average2) / (a2Count + 1);
					//a2Count++;
					if (boundaryMatrix->isZeroColumn(j) == false){
						pivotJ = boundaryMatrix->getPivotRow(j);
					}
				}
				
				/*If the column j gets to this point and is still not 0, then it must have a pivot that can't be canceled.
				In which case, the simplex corresponding to the pivot row and the the simplex corresponding to column j form
				a persistence pair*/
				//start = chrono::system_clock::now();
				if (boundaryMatrix->isZeroColumn(j) == false){
					/*Set column j as the column with row pivotJ as its pivot in the reduced matrix*/
					unsigned int i = pivotJ;
					pivots[i] = j;
					boundaryMatrix->clearColumn(i);

					/*Add persistence pair*/
					Simplex *s1 = this->filtration[i];
					Simplex *s2 = this->filtration[j];
					if (s1->dim == 0){
						Vertex *v = (Vertex*)s1;
						Edge *e = (Edge*)s2;
						double persistence = e->funcValue - v->funcValue;
						double symPerturb = e->filtrationPosition - v->filtrationPosition;

						persistencePair01 pp = { v, e, persistence, symPerturb };
						this->msPersistencePairs.push_back(pp);
					}
					else{
						Edge* e = (Edge*)s1;
						Triangle* t = (Triangle*)s2;
						double persistence = t->funcValue - e->funcValue;
						double symPerturb1 = t->filtrationPosition - e->filtrationPosition;
						double symPerturb2 = t->filtrationPosition;

						persistencePair12 pp = { e, t, persistence, symPerturb1, symPerturb2 };
						this->smPersistencePairs.push_back(pp);
					}
				}
				//end = chrono::system_clock::now();
				//elapsedDuration = end - start;
				//elapsedSeconds = elapsedDuration.count();
				//average3 = (elapsedSeconds + a3Count * average3) / (a3Count + 1);
				//a3Count++;
			}
		}
		cout<<"done\n";
	}
	for (vector<Vertex*>::iterator it = this->K->vBegin(); it != this->K->vEnd(); it++) {
		K->addCriticalPoint((Simplex*)*it);
	}
	for (vector<Edge*>::iterator it = this->K->eBegin(); it != this->K->eEnd(); it++) {
		K->addCriticalPoint((Simplex*)*it);
	}
	for (vector<Triangle*>::iterator it = this->K->tBegin(); it != this->K->tEnd(); it++) {
		K->addCriticalPoint((Simplex*)*it);
	}

	delete boundaryMatrix;
	/*for (vector<persistencePair01>::iterator it = this->msPersistencePairs.begin(); it != this->msPersistencePairs.end(); it++) {
		persistencePair01 p = *it;
		this->K->addCriticalPoint(p.min);
		this->K->addCriticalPoint(p.saddle);
	}*/
	cout << "# of min-saddle pairs: " << this->msPersistencePairs.size() << "\n";
	cout << "# of saddle-max pairs: " << this->smPersistencePairs.size() << "\n";
	int count = 0;
	for (vector<Edge*>::iterator it = K->eBegin(); it != K->eEnd(); it++) {
		count++;
	}
	cout << "# of edges in complex: " << count << "\n";
	cout << "average time to check zero column: " << average1 << "s\n";
	cout << "average time to add columns: " << average2 << "s\n";
	cout << "average time to add persistence pairs: " << average3 << "s\n";
}

//void PersistencePairs::computePersistencePairsSeparate() {
//	int count = 0;
//	for (vector<Simplex*>::iterator it = this->filtration.begin(); it != this->filtration.end(); it++) {
//		Simplex *s = *it;
//		if (s->dim == 1 || s->dim == 2) {
//			count++;
//		}
//	}
//
//	SparseZ2Matrix *boundaryMatrix1 = new SparseZ2Matrix(count, count);
//	int idx = 0;
//	for (unsigned int i = 0; i < this->filtration.size(); i++) {
//		Simplex *s = filtration[i];
//		if (s->dim == 1) {
//			Edge *e = (Edge*)s;
//
//		}
//		else if (s->dim == 2) {
//
//		}
//	}
//
//}

//test cancellability
vector<Simplex*>* PersistencePairs::isCancellable(const persistencePair01& pp, ofstream& cancelData){
	DiscreteVField *V = this->K->getDiscreteVField();
	Simplicial2Complex *K = this->K;

	Edge *e = pp.saddle;
	Vertex *v = pp.min;

	//if they have 0 persistence, we know they are cancellable and take this shortcut
	if (e->funcValue == v->funcValue && v->hasEdge(e)) {
		cancelData << "Index: 01, Persistence: " << pp.persistence << " + " << pp.symPerturb << "e" << " Cancellable: Yes (trivial)" << endl;
		return new vector<Simplex*>({ (Simplex*)e, (Simplex*)v });
	}

	/*Each vector<Simplex*> is a V-path starting with Edge e and ending on a minimum Vertex
	If, after performing a DFS, there are two V-paths from e to v, then we cannot cancel them.
	If there are no V-paths from e to v, then we cannot cancel them. The only time we would be able to cancel
	is if there is exactly one V-path from e to v*/
	vector<vector<Simplex*>*> paths;

	/*An edge-vertex V-path is not capable of branching beyond the first edge, so we start with two paths*/
	vector<Simplex*> *path1 = new vector<Simplex*>({ (Simplex*)e });
	vector<Simplex*> *path2 = new vector<Simplex*>({ (Simplex*)e });
	paths.push_back(path1);
	paths.push_back(path2);

	stack<Simplex*> *st = new stack<Simplex*>();
	Vertex *v1 = get<0>(e->getVertices());
	Vertex *v2 = get<1>(e->getVertices());
	paths[0]->push_back((Simplex*)v1);
	paths[1]->push_back((Simplex*)v2);
	st->push((Simplex*)v1);
	st->push((Simplex*)v2);

	while (!st->empty()){
		Simplex* s = st->top();
		st->pop();
		if (s->dim == 0){
			Vertex *vert = (Vertex*)s;
			Edge* paired_edge = V->containsPair(vert);
			if (paired_edge!=NULL){
				//if (pp.persistence == 0) cout << "This should never happen VE\n";
				/*If there is one found, add the edge to all paths ending with vert*/
				for (int i = 0;(unsigned) i < paths.size(); i++){
					vector<Simplex*> *path = paths[i];
					if (path->back() == s){
						path->push_back((Simplex*) paired_edge);
					}
				}
				/*Then push the edge onto the stack*/
				st->push((Simplex*) paired_edge);
			}
			/*If there is no way out of the Vertex, then we've reached a minimum*/
		}
		else{
			Edge *edge = (Edge*)s;
			/*Here, we only need to check two vertices*/
			Vertex *vert1 = get<0>(edge->getVertices());
			Vertex *vert2 = get<1>(edge->getVertices());

			/*We want the vertex such that the discrete gradient vector field does NOT contain
			(vertex, edge)*/
			Vertex *exitVert;
			
			if ((V->containsPair(vert1)) != edge){
				exitVert = vert1;
			}
			else{
				exitVert = vert2;
			}

			/*exitVert is the next on our path, so we add it to all paths ending with edge*/
			for (int i = 0;(unsigned) i < paths.size(); i++){
				vector<Simplex*> *path = paths[i];
				if (path->back() == s){
					path->push_back((Simplex*)exitVert);
				}
			}
			/*Then we need to push it onto the stack*/
			st->push((Simplex*)exitVert);
		}
	}

	/*Now if there are none or more than one paths ending at Vertex v, the persistence
	pair is not cancellable, and we return an empty vector*/
	bool hasPath = false;
	vector<Simplex*> *uniquePath = new vector<Simplex*>();
	for (int i = 0;(unsigned) i < paths.size(); i++){
		vector<Simplex*> *path = paths[i];
		if (path->back() == (Simplex*)v && hasPath == true){
			cancelData << "Index: 01, Persistence: " << pp.persistence << " + " << pp.symPerturb << "e" << ", Cancellable: No, Reason: path not unique\n";
			return new vector<Simplex*>();
		}
		else if (path->back() == (Simplex*)v){
			hasPath = true;
			uniquePath = path;
		}
	}
	/*If hasPath is true and the for loop didn't return as soon as it found a 2nd,
	Then uniquePath truly is unique. So we return it*/

	if (hasPath == true){
		cancelData << "Index: 01, Persistence: " << pp.persistence << " + " << pp.symPerturb << "e" << ", Cancellable: Yes\n";
		return uniquePath;
	}
	else{
		/*Otherwise, there is no path from e to v*/
		cancelData << "Index: 01, Persistence: " << pp.persistence << " + " << pp.symPerturb << "e" << ", Cancellable: No, Reason: no path found\n";
		return new vector<Simplex*>();
	}
}

//test cancellability
vector<Simplex*>* PersistencePairs::isCancellable(const persistencePair12& pp, ofstream& cancelData){
	DiscreteVField *V = this->K->getDiscreteVField();
	Simplicial2Complex* K = this->K;

	Triangle *max = pp.max;
	Edge* saddle = pp.saddle;

	if (max->funcValue == saddle->funcValue && saddle->hasTriangle(max)) {
		cancelData << "Index: 12, Persistence: "<< pp.persistence << " + " << pp.symPerturb1 << "e" << " Cancellable: Yes (trivial)" << endl;
		return new vector<Simplex*>({ (Simplex*)max,(Simplex*)saddle });
	}

	vector<vector<Simplex*>*> paths;
	paths.push_back(new vector<Simplex*>({ (Simplex*)max }));
	paths.push_back(new vector<Simplex*>({ (Simplex*)max }));
	paths.push_back(new vector<Simplex*>({ (Simplex*)max }));

	stack<Simplex*> *st = new stack<Simplex*>();
	Edge *edge1 = get<0>(max->getEdges());
	Edge *edge2 = get<1>(max->getEdges());
	Edge *edge3 = get<2>(max->getEdges());

	paths[0]->push_back((Simplex*)edge1);
	paths[1]->push_back((Simplex*)edge2);
	paths[2]->push_back((Simplex*)edge3);

	st->push((Simplex*)edge1);
	st->push((Simplex*)edge2);
	st->push((Simplex*)edge3);

	while (!st->empty()){
		Simplex* s = st->top();
		st->pop();
		if (s->dim == 1){
			Edge *e = (Edge*)s;
			Triangle* paired_triangle = V->containsPair(e);
			if (paired_triangle != NULL){
				for (int i = 0; i < paths.size(); i++){
					vector<Simplex*> *path = paths[i];
					if (path->back() == s){
						path->push_back((Simplex*) paired_triangle);
					}
				}
				st->push((Simplex*) paired_triangle);
			}
		}
		else{
			Triangle *t = (Triangle*)s;
			Edge *e1 = get<0>(t->getEdges());
			Edge *e2 = get<1>(t->getEdges());
			Edge *e3 = get<2>(t->getEdges());

			Edge *exitEdge1, *exitEdge2;
			/*if (V->containsPair(e1, t) && V->containsPair(e3, t)){
				cout << "Problem\n";
				Sleep(1000);
				exit(0);
			}*/
			if (V->containsPair(e1)==t){
				exitEdge1 = e2;
				exitEdge2 = e3;
			}
			else if (V->containsPair(e2)==t){
				exitEdge1 = e1;
				exitEdge2 = e3;
			}
			else{
				exitEdge1 = e1;
				exitEdge2 = e2;
			}
			for (int i = 0; i < paths.size(); i++){
				vector<Simplex*> *path = paths[i];
				/*The path branches into two paths, so we need to make a copy of the path for the second branch*/
				if (path->back() == s){
					/*Copy the path*/
					vector<Simplex*> *pathCopy = new vector<Simplex*>(path->begin(), path->end());
					/*Put exitEdge1 in the original path*/
					path->push_back((Simplex*)exitEdge1);
					/*Put exitEdge2 in the copy*/
					pathCopy->push_back((Simplex*)exitEdge2);
					
					paths.push_back(pathCopy);
				}
			}

			/*Now that the paths are updated, we push exitEdge1 and exitEdge2 onto the stack*/
			st->push(exitEdge1);
			st->push(exitEdge2);
		}
	}

	bool hasPath = false;
	vector<Simplex*> *uniquePath = new vector<Simplex*>();
	for (int i = 0; i < paths.size(); i++){
		vector<Simplex*> *path = paths[i];
		if (path->back() == (Simplex*)saddle && hasPath == true){
			cancelData << "Index: 12, Persistence: " << pp.persistence << " + " << pp.symPerturb1 << "e" << ", Cancellable: No, Reason: path not unique\n";
			return new vector<Simplex*>();
		}
		else if (path->back() == (Simplex*)saddle){
			hasPath = true;
			uniquePath = path;
		}
	}

	if (hasPath == true){
		cancelData << "Index: 12, Persistence: " << pp.persistence << " + " << pp.symPerturb1 << "e" << ", Cancellable: Yes\n";
		return uniquePath;
	}
	else{
		cancelData << "Index: 12, Persistence: " << pp.persistence << " + " << pp.symPerturb1 << "e" << ", Cancellable: No, Reason: No path found\n";
		return new vector<Simplex*>();
	}
}

//performs cancellation
void PersistencePairs::cancelAlongVPath(vector<Simplex*>* VPath){
	DiscreteVField *V = this->K->getDiscreteVField();
	/*As long as VPath is not empty, we may assume it has at least 2 entries*/
	/*Simplex* s = VPath->at(0);
	if (s->dim == 1){
		Simplex *first = s;
		for (int i = 1; i < VPath->size(); i++){
			Simplex* next = VPath->at(i);
			if (first->dim == 0){
				V->removePair((Vertex*)first, (Edge*)next);
			}
			else if (first->dim == 1){
				V->addPair((Vertex*)next, (Edge*)first);
			}
			first = next;
		}
	}
	else if(s->dim == 2){
		Simplex* first = s;
		for (int i = 1; i < VPath->size(); i++){
			Simplex* next = VPath->at(i);
			if (first->dim == 1){
				V->removePair((Edge*)first, (Triangle*)next);
			}
			else if (first->dim == 2){
				V->addPair((Edge*)next, (Triangle*)first);
			}
			first = next;
		}
	}*/

	if (VPath->at(0)->dim == 1){
		vector<Simplex*> * jp = new vector<Simplex*>();
		for (int i = 0; i < VPath->size(); i++){
			Simplex *s = VPath->at(i);
			jp->push_back(s);
			if (s->dim == 0){
				Vertex *v = (Vertex*)s;
				if (i < VPath->size() - 1){
					V->removePair(v, (Edge*)VPath->at(i + 1));
				}
				if (i > 0){
					V->addPair(v, (Edge*)VPath->at(i - 1));
				}
			}
		}
		//V->addPair((Vertex*)(jp->back()), (Edge*)(jp->front()));
		//V->addJump(jp->front(), jp);
	}
	else if (VPath->at(0)->dim == 2){
		vector<Simplex*> * jp = new vector<Simplex*>();
		for (int i = 0; i < VPath->size(); i++){
			Simplex *s = VPath->at(i);
			jp->push_back(s);
			if (s->dim == 1){
				Edge *e = (Edge*)s;
				if (i < VPath->size() - 1){
					V->removePair(e, (Triangle*)VPath->at(i + 1));
				}
				if (i > 0){
					V->addPair(e, (Triangle*)VPath->at(i - 1));
				}
			}
		}
	}
	else{
		cout << "This shouldn't happen\n";
	}
}


//cancels al pairs that can be cancelled
void PersistencePairs::cancelPersistencePairs(double delta){
	

	cout << "\tSorting persistence pairs...\n";
	sort(msPersistencePairs.begin(), msPersistencePairs.end(), PersistencePairs::persistencePairCompare01);
	sort(smPersistencePairs.begin(), smPersistencePairs.end(), PersistencePairs::persistencePairCompare12);
	cout << "\tDone\n";

	this->outputPersistencePairs("persistencePairs.txt");

	/*cout << "\tOutputting persistence pairs...\n";
	ofstream persistencePairs("persistencePairs.txt", ios_base::trunc | ios_base::out);
	int idx1 = 0, idx2 = 0;
	while (idx1 < this->msPersistencePairs.size() && idx2 < this->smPersistencePairs.size()){
		persistencePair01 pair1 = this->msPersistencePairs[idx1];
		persistencePair12 pair2 = this->smPersistencePairs[idx2];

		if (pair1.persistence <= pair2.persistence){
			persistencePairs << "Vertex Edge " << pair1.persistence<< "\n";
			idx1++;
		}
		else{
			persistencePairs << "Edge Triangle " << pair2.persistence << "\n";
			idx2++;
		}
	}

	if (idx1 < this->msPersistencePairs.size()){
		while (idx1 < this->msPersistencePairs.size()){
			persistencePair01 pair = this->msPersistencePairs[idx1];
			persistencePairs << "Vertex Edge " << pair.persistence << "\n";
			idx1++;
		}
	}
	else if (idx2 < this->smPersistencePairs.size()){
		while (idx2 < this->smPersistencePairs.size()){
			persistencePair12 pair = this->smPersistencePairs[idx2];
			persistencePairs << "Edge Triagle " << pair.persistence << "\n";
			idx2++;
		}
	}
	cout << "\tDone\n";*/

	/*cout << "\tCancelling min-saddle pairs...\n";
	int i = 0;
	while (i < this->msPersistencePairs.size()){
		//cout << "Handling pair " << i << "\n";
		persistencePair01 pair = msPersistencePairs[i];
		if (pair.persistence <= delta){
			vector<Simplex*>* VPath = this->isCancellable(pair);
			if (!VPath->empty()){
				this->cancelAlongVPath(VPath);
				this->K->removeCriticalPoint(pair.min);
				this->K->removeCriticalPoint(pair.saddle);
			}
			i++;
		}
		else{
			break;
		}
	}
	cout << "\tDone\n";

	cout << "\tCancelling saddle-max pairs...\n";
	int j = 0;
	while (j < this->smPersistencePairs.size()){
		persistencePair12 pair = this->smPersistencePairs[j];
		if (pair.persistence <= delta){
			vector<Simplex*>* VPath = this->isCancellable(pair);
			if (!VPath->empty()){
				this->cancelAlongVPath(VPath);
				this->K->removeCriticalPoint(pair.saddle);
				this->K->removeCriticalPoint(pair.max);
			}
			j++;
		}
		else{
			break;
		}
	}
	cout << "Done\n";*/

	DiscreteVField *V = this->K->getDiscreteVField();
	

	ofstream cancelData("cancellationData.txt", ios_base::trunc | ios_base::out);

	cout << "\tCancelling...\n";
	int i = 0, j = 0;
	double average4 = 0, average5 = 0;
	int a4Count = 0, a5Count = 0;
	int count = 0;
	//chrono::time_point<chrono::system_clock> start, end;
	cout << "msPair:" << msPersistencePairs.size() << "\tsmPair" << smPersistencePairs.size() << endl;
	int Pcounter = 0, verbose = 0;
	while (i < msPersistencePairs.size() && j < smPersistencePairs.size()){
		persistencePair01 pair1 = msPersistencePairs[i];
		persistencePair12 pair2 = smPersistencePairs[j];
		
		if (pair1.persistence < pair2.persistence || pair1.persistence == pair2.persistence && pair1.symPerturb < pair2.symPerturb1){
			if (pair1.persistence <= delta){
				//start = chrono::system_clock::now();
				if (verbose){
					cout<<"checking cancellable for pair1"<<endl;
				}
				vector<Simplex*> *VPath = this->isCancellable(pair1, cancelData);

				//end = chrono::system_clock::now();
				//chrono::duration<double> eDur = end - start;
				//double elapsedTime = eDur.count();
				//average4 = (elapsedTime + a4Count * average4) / (a4Count + 1);
				//a4Count++;
				if (verbose){
					cout<<"removing path"<<endl;
				}
				if (!VPath->empty()){
					count++;
					this->cancelAlongVPath(VPath);
					K->removeCriticalPoint(pair1.min);
					K->removeCriticalPoint(pair1.saddle);
				}
			}
			i++;
		}
		else{
			if (pair2.persistence <= delta){
				if (verbose){
					cout<<"checking cancellable for pair 2"<<endl;
				}
				//start = chrono::system_clock::now();
				vector<Simplex*> *VPath = this->isCancellable(pair2, cancelData);
				/*end = chrono::system_clock::now();
				chrono::duration<double> eDur = end - start;
				double elapsedTime = eDur.count();
				average4 = (elapsedTime + a4Count * average4) / (a4Count + 1);
				a4Count++;*/
				if (verbose){
					cout<<"removing path"<<endl;
				}
				if (!VPath->empty()){
					count++;
					this->cancelAlongVPath(VPath);
					K->removeCriticalPoint(pair2.saddle);
					K->removeCriticalPoint(pair2.max);
				}
			}
			j++;
		}
		
		
		Pcounter++;
		if (DEBUG && Pcounter%10000==0){
			cout << "\r";
			cout << "\t" << Pcounter << "/" << msPersistencePairs.size() + smPersistencePairs.size() << "...";
		}//else if (DEBUG && Pcounter == 802230) {
			//cout<< Pcounter << "/" << msPersistencePairs.size() + smPersistencePairs.size() << "...\n";
			//verbose = 1;
		//}
	}
	
	cout << "msPair: " << i << "/" << msPersistencePairs.size() <<endl;
	cout << "smPair: " << j << "/" << smPersistencePairs.size() <<endl;
	
	if (i < msPersistencePairs.size()){
		while (i < msPersistencePairs.size()){
			persistencePair01 pair = msPersistencePairs[i];
			if (pair.persistence <= delta){
				//start = chrono::system_clock::now();
				vector<Simplex*> *VPath = this->isCancellable(pair, cancelData);
				/*end = chrono::system_clock::now();
				chrono::duration<double> eDur = end - start;
				double elapsedTime = eDur.count();
				average4 = (elapsedTime + a4Count * average4) / (a4Count + 1);
				a4Count++;*/
				if (!VPath->empty()){
					count++;
					this->cancelAlongVPath(VPath);
					K->removeCriticalPoint(pair.min);
					K->removeCriticalPoint(pair.saddle);
				}
			}
			i++;
			if (DEBUG && i % 10 == 0){
				cout << "\r";
				cout << "\t" << i << "/" << msPersistencePairs.size() <<endl;
			}
		}
	}
	else if (j < smPersistencePairs.size()){
		while (j < smPersistencePairs.size()){
			persistencePair12 pair = smPersistencePairs[j];
			if (pair.persistence <= delta){
				//start = chrono::system_clock::now();
				vector<Simplex*> *VPath = this->isCancellable(pair, cancelData);
				/*end = chrono::system_clock::now();
				chrono::duration<double> eDur = end - start;
				double elapsedTime = eDur.count();
				average4 = (elapsedTime + a4Count * average4) / (a4Count + 1);
				a4Count++;*/
				if (!VPath->empty()){
					count++;
					this->cancelAlongVPath(VPath);
					K->removeCriticalPoint(pair.saddle);
					K->removeCriticalPoint(pair.max);
				}
			}
			
			j++;
			if (DEBUG && j % 10 == 0){
				cout << "\r";
				cout << "\t" << j << "/" << smPersistencePairs.size() <<endl;
			}
		}
	}
	cout << "Cancelled Pair: " <<count<<endl;

	/*for (vector<persistencePair01>::iterator it = this->msPersistencePairs.begin(); it != this->msPersistencePairs.end(); it++) {
		persistencePair01 p = *it;
		if (p.persistence <= delta) {
			vector<Simplex*> *VPath = this->isCancellable(p, cancelData);
			if (!VPath->empty()) {
				this->cancelAlongVPath(VPath);
				K->removeCriticalPoint(p.min);
				K->removeCriticalPoint(p.saddle);
			}
		}
	}*/

	cout << "\tDone\n";
	//cout << "average time to determine cancellability: " << average4 << "s\n";

	cout << "\tOutputting discrete gradient vector field...\n";
	// might need put this back
	// V->outputVEmap();
	// V->outputETmap();

	/*Edge *topEdge = msPersistencePairs[msPersistencePairs.size() * 2 / 3].saddle;
	set<Simplex*> *descMan = this->K->descendingManifold(topEdge);
	ofstream vertexIndices("vertices.txt");
	ofstream edgeIndices("MATLABedges.txt");
	vertexIndices << "index\n";
	edgeIndices << "index1 index2\n";
	for (set<Simplex*>::iterator it = descMan->begin(); it != descMan->end(); it++) {
		if ((*it)->dim == 0) {
			Vertex *v = (Vertex*)*it;
			vertexIndices << v->getVPosition() + 1 << "\n";
		}
		else if ((*it)->dim == 1) {
			Edge *e = (Edge*)*it;
			Vertex *v1 = get<0>(e->getVertices());
			Vertex *v2 = get<1>(e->getVertices());
			edgeIndices << v1->getVPosition() + 1 << " " << v2->getVPosition() + 1 << "\n";
		}
	}
	vertexIndices.close();
	edgeIndices.close();*/
	cout << "\tDone\n";
	cancelData.close();
}
