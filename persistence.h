#include "Simplicial2Complex.h"
#include "sparseZ2matrix.h"
#include <list>
#include <stack>
//#include <windows.h>
//#include <chrono>
#include <iostream>
#include <algorithm>

#include <phat/compute_persistence_pairs.h>

#include <phat/representations/vector_vector.h>
#include <phat/representations/vector_heap.h>
#include <phat/representations/vector_set.h>
#include <phat/representations/vector_list.h>
#include <phat/representations/sparse_pivot_column.h>
#include <phat/representations/heap_pivot_column.h>
#include <phat/representations/full_pivot_column.h>
#include <phat/representations/bit_tree_pivot_column.h>

#include <phat/algorithms/twist_reduction.h>
#include <phat/algorithms/standard_reduction.h>
#include <phat/algorithms/row_reduction.h>
#include <phat/algorithms/chunk_reduction.h>
#include <phat/algorithms/spectral_sequence_reduction.h>

#include <phat/helpers/dualize.h>

using namespace std;

typedef struct {
	Vertex *min;
	Edge *saddle;
	double persistence;
	double symPerturb;
	double loc_diff;
}persistencePair01;
typedef struct {
	Edge *saddle;
	Triangle *max;
	double persistence;
	double symPerturb1;
	double symPerturb2;
	double loc_diff;
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
		else if (fabs(p.persistence - q.persistence) < 1e-8 && p.symPerturb < q.symPerturb){
			return true;
		}else if(fabs(p.persistence - q.persistence) < 1e-8 && fabs(p.symPerturb - q.symPerturb) < 1e-8 && p.loc_diff < q.loc_diff){
			//cout<<"caught something"<<endl;
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
		} else if (fabs(p.persistence - q.persistence) < 1e-8 && p.symPerturb1 < q.symPerturb1){
			return true;
		} else if (fabs(p.persistence - q.persistence) < 1e-8 && fabs(p.symPerturb1 - q.symPerturb1) < 1e-8 && p.symPerturb2 < q.symPerturb2){
			return true;
		} else if (fabs(p.persistence - q.persistence) < 1e-8 && fabs(p.symPerturb1 - q.symPerturb1) < 1e-8 && fabs(p.symPerturb2 - q.symPerturb2) < 1e-8 && p.loc_diff < q.loc_diff){
			// cout << "caught something" <<endl;
			return true;
		} else{
			return false;
		}
	}
	static bool persistencePairCompare(const persistencePair12& p, const persistencePair01& q){
		return p.persistence < q.persistence;
	}

	PersistencePairs(Simplicial2Complex *K){
		this->K = K;
		msPersistencePairs.clear();
		smPersistencePairs.clear();
	}
	void buildFiltration();

	void computePersistencePairsWithClear();
	void PhatPersistence();

	vector<Simplex*>* isCancellable(const persistencePair01&, ofstream&);
	vector <Simplex*>* isCancellable(const persistencePair12&, ofstream&);
	/*Requires that VPath was returned by isCancellable, and no calls to
	cancelAlongPath have been made between then and this call*/
	void cancelAlongVPath(vector<Simplex*>* VPath);

	void cancelPersistencePairs(double delta);

	void outputPersistencePairs(string pathname){
		ofstream outst(pathname, ios_base::trunc | ios_base::out);
		for (int i = 0; i < this->msPersistencePairs.size(); i++){
			outst << "Vertex Edge " << this->msPersistencePairs[i].persistence << " " << this->msPersistencePairs[i].symPerturb << "\n";
			//this->msPersistencePairs[i].min->output();
			//this->msPersistencePairs[i].saddle->output();
			//cout << "\n";
		}
		for (int i = 0; i < this->smPersistencePairs.size(); i++){
			outst << "Edge Triangle " << this->smPersistencePairs[i].persistence << " " << this->smPersistencePairs[i].symPerturb1 << " " << this->smPersistencePairs[i].symPerturb2 << "\n";
			//this->smPersistencePairs[i].saddle->output();
            //this->smPersistencePairs[i].max->output();
			//cout <<"\n";
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


void PersistencePairs::PhatPersistence(){
	// generate boundary matrix
	cout << "\tInitializing boundary matrix...\n";
	cout << "\t\tMatrix size: " << this->filtration.size() << "\n";
	phat::boundary_matrix< phat::vector_vector > boundary_matrix;
	boundary_matrix.set_num_cols(this->filtration.size());
	
	std::vector< phat::index > temp_col;
	for (unsigned int i = 0; i < this->filtration.size(); i++){
		Simplex *s = filtration[i];
		temp_col.clear();
				
		if (s->dim == 0){
			// no boundary. no bits to set in matrix
			boundary_matrix.set_dim( i, 0 );
		}
		else if (s->dim == 1){
			// edge
			boundary_matrix.set_dim( i, 1 );
			Edge *e = (Edge*)s;
			Vertex *v1 = get<0>(e->getVertices());
			Vertex *v2 = get<1>(e->getVertices());
			if (v1->filtrationPosition < v2->filtrationPosition) {
				temp_col.push_back(v1->filtrationPosition);
				temp_col.push_back(v2->filtrationPosition);				
			}else{
				temp_col.push_back(v2->filtrationPosition);
				temp_col.push_back(v1->filtrationPosition);		
			}
		}
		else if (s->dim == 2){
			// triangle
			boundary_matrix.set_dim( i, 2 );
			Triangle *t = (Triangle*)s;
			Edge *e1 = get<0>(t->getEdges());
			Edge *e2 = get<1>(t->getEdges());
			Edge *e3 = get<2>(t->getEdges());
			Edge *edges[3] = { e1,e2,e3 };
			for(int j = 0; j < 2; j++){
				for(int k = j + 1; k < 3; k++){
					if(edges[k]->filtrationPosition \
					   < edges[j]->filtrationPosition){
						swap(edges[k], edges[j]);
					}
				}
			}
			temp_col.push_back(edges[0]->filtrationPosition);
			temp_col.push_back(edges[1]->filtrationPosition);
			temp_col.push_back(edges[2]->filtrationPosition);
		}
		boundary_matrix.set_col( i, temp_col );
	}
	cout << "\tInitialized!\n";
	
	// call Phat
	cout << "\tComputing persitence pairs...\n";
	phat::persistence_pairs pairs;
	phat::compute_persistence_pairs\
		< phat::twist_reduction >( pairs, boundary_matrix );
	pairs.sort();
	cout << "\tComputed and sorted!\n";
	
	// post processing: add sm-ms pairs, set critical points
	cout << "\tCounting ms-pairs and sm-pairs...";
	for( phat::index idx = 0; idx < pairs.get_num_pairs(); idx++ ){
		Simplex *s1 = this->filtration[pairs.get_pair( idx ).first];
		Simplex *s2 = this->filtration[pairs.get_pair( idx ).second];
		if (s1->dim == 0){
			Vertex *v = (Vertex*)s1;
			Edge *e = (Edge*)s2;
			double persistence = e->funcValue - v->funcValue;
			double symPerturb = e->getSymPerturb();
			double loc_diff = e->filtrationPosition - v->filtrationPosition;
			int c_type = 1;
			e->set_type(c_type);
			persistencePair01 pp = { v, e, persistence, symPerturb, loc_diff };
			this->msPersistencePairs.push_back(pp);
		}
		else{
			Edge* e = (Edge*)s1;
			Triangle* t = (Triangle*)s2;
			double persistence = t->funcValue - e->funcValue;
			double symPerturb1 = get<0>(t->getSymPerturb())\
							     - e->getSymPerturb();
			double symPerturb2 = get<1>(t->getSymPerturb());
			double loc_diff = t->filtrationPosition - e->filtrationPosition;
			int c_type = 2;
			e->set_type(c_type);
			persistencePair12 pp = { e, t, persistence, symPerturb1,\
									symPerturb2, loc_diff };
			this->smPersistencePairs.push_back(pp);
		}
	}
	cout << "done!\n";
	for (vector<Vertex*>::iterator it = this->K->vBegin();\
		 it != this->K->vEnd(); it++) {
		K->addCriticalPoint((Simplex*)*it);
	}
	for (vector<Edge*>::iterator it = this->K->eBegin();\
		 it != this->K->eEnd(); it++) {
		K->addCriticalPoint((Simplex*)*it);
	}
	for (vector<Triangle*>::iterator it = this->K->tBegin();\
		 it != this->K->tEnd(); it++) {
		K->addCriticalPoint((Simplex*)*it);
	}
	cout << "\tCritical points initialized\n";
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
				cout.flush();
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
						double symPerturb = e->getSymPerturb();
						double loc_diff = e->filtrationPosition - v->filtrationPosition;
						int c_type = 1;
						e->set_type(c_type);
						persistencePair01 pp = { v, e, persistence, symPerturb, loc_diff };
						this->msPersistencePairs.push_back(pp);
					}
					else{
						Edge* e = (Edge*)s1;
						Triangle* t = (Triangle*)s2;
						double persistence = t->funcValue - e->funcValue;
						double symPerturb1 = get<0>(t->getSymPerturb()) - e->getSymPerturb();
						double symPerturb2 = get<1>(t->getSymPerturb());
						double loc_diff = t->filtrationPosition - e->filtrationPosition;
						int c_type = 2;
						e->set_type(c_type);
						persistencePair12 pp = { e, t, persistence, symPerturb1, symPerturb2, loc_diff };
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

//test cancellability
vector<Simplex*>* PersistencePairs::isCancellable(const persistencePair01& pp, ofstream& cancelData){
	DiscreteVField *V = this->K->getDiscreteVField();
	Simplicial2Complex *K = this->K;

	Edge *e = pp.saddle;
	Vertex *v = pp.min;

	//if they have 0 persistence, we know they are cancellable and take this shortcut
	if (fabs(e->funcValue - v->funcValue)<1e-8 && v->hasEdge(e)) {
        if (DEBUG) {
            cancelData << "Index: 01, Persistence: " << pp.persistence << " + " << pp.symPerturb << "e" << " Cancellable: Yes (trivial)" << endl;
		}
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
		cancelData << "Index: 01, Persistence: " << pp.persistence << " + " << pp.symPerturb << "e" << ", Cancellable: No, Reason: no path found";
		return new vector<Simplex*>();
	}
}

//test cancellability

vector<Simplex*>* PersistencePairs::isCancellable(const persistencePair12& pp, ofstream& cancelData){
	DiscreteVField *V = this->K->getDiscreteVField();
	Simplicial2Complex* K = this->K;
	unordered_set<Simplex*> visited_once;
	unordered_set<Simplex*> visited_twice;
	//mark visited or not
	visited_once.clear();visited_twice.clear();

	Triangle *max = pp.max;
	Edge* saddle = pp.saddle;

	if (fabs(max->funcValue - saddle->funcValue) < 1e-8 && saddle->hasTriangle(max)) {
        if (DEBUG){
			cancelData << "Index: 12, Persistence: "<< pp.persistence << " + " << pp.symPerturb1 << "e " << pp.symPerturb2 << " Cancellable: Yes (trivial)" << endl;
        }
		
		return new vector<Simplex*>({ (Simplex*)max,(Simplex*)saddle });
	}
	
	list<Simplex*> BFSqueue;
	BFSqueue.clear();
	Simplex* curr = (Simplex*)max;
	curr->prev = NULL;
	BFSqueue.push_back(curr);
	visited_once.insert(curr);

	while (!BFSqueue.empty()){
		curr = BFSqueue.front();
		BFSqueue.pop_front();
		if (curr->dim == 1){
			Edge *e = (Edge*)curr;
			Triangle* paired_triangle = V->containsPair(e);
			if (paired_triangle != NULL){
                Simplex* simp_triangle = (Simplex*) paired_triangle;
				if (visited_twice.count(simp_triangle)>0){
					// if (DEBUG) cout << "Twice+";
				}else if (visited_once.count(simp_triangle) > 0){
					// if (DEBUG) cout << "Twice!";
					visited_twice.insert(simp_triangle);
				}else{
					// if (DEBUG) cout << "Once";
				    simp_triangle->prev = curr;
					visited_once.insert(simp_triangle);
					BFSqueue.push_back(simp_triangle);
				}
			}
			// if (DEBUG) cout << endl;
		}
		else{
			Triangle *t = (Triangle*)curr;
			Edge *e1 = get<0>(t->getEdges());
			Edge *e2 = get<1>(t->getEdges());
			Edge *e3 = get<2>(t->getEdges());

            Simplex* simp_edge = (Simplex*) e1;
            if (visited_once.count(simp_edge) == 0){
                simp_edge->prev = curr;
                visited_once.insert(simp_edge);
                BFSqueue.push_back(simp_edge);
				//if (DEBUG) cout << "@found: \tT" << simp_edge->funcValue << "Once";
            }else if (visited_twice.count(simp_edge) == 0){
                if (curr->prev != simp_edge){
					visited_twice.insert(simp_edge);
					//if (DEBUG) cout << "@found: \tT" << simp_edge->funcValue << "Twice!";
				}
            }

            simp_edge = (Simplex*) e2;
            if (visited_once.count(simp_edge) == 0){
                simp_edge->prev = curr;
                visited_once.insert(simp_edge);
                BFSqueue.push_back(simp_edge);
				//if (DEBUG) cout << "@found: \tT" << simp_edge->funcValue << "Once";
            }else if (visited_twice.count(simp_edge) == 0){
                if (curr->prev != simp_edge){
					visited_twice.insert(simp_edge);
					//if (DEBUG) cout << "@found: \tT" << simp_edge->funcValue << "twice!";
				}
            }

            simp_edge = (Simplex*) e3;
            if (visited_once.count(simp_edge) == 0){
                simp_edge->prev = curr;
                visited_once.insert(simp_edge);
                BFSqueue.push_back(simp_edge);
				//if (DEBUG) cout << "@found: \tT" << simp_edge->funcValue << "Once";
            }else if (visited_twice.count(simp_edge) == 0){
                if (curr->prev != simp_edge){
					visited_twice.insert(simp_edge);
					//if (DEBUG) cout << "@found: \tT" << simp_edge->funcValue << "twice!";
				}
            }

			//if (DEBUG) cout << endl;
		}
	}

	vector<Simplex*> *uniquePath = new vector<Simplex*>();
	list<Simplex*> path;
	path.clear();
    curr = (Simplex*)saddle;
    if (visited_twice.count(curr) > 0){
        if (DEBUG) cancelData << "Index: 12, Persistence: " << pp.persistence << " + " << pp.symPerturb1 << "e" << ", Cancellable: No, Reason: path not unique1\n";
        return new vector<Simplex*>();
    }else if (visited_once.count(curr) > 0){
        while (curr!=NULL){
            if (visited_twice.count(curr) > 0){
                if (DEBUG) cancelData << "Index: 12, Persistence: " << pp.persistence << " + " << pp.symPerturb1 << "e" << ", Cancellable: No, Reason: path not unique2\n";
                return new vector<Simplex*>();
            }
            path.push_front(curr);
            curr = curr->prev;
        }
		while(!path.empty()){
			uniquePath->push_back(path.front());
			path.pop_front();
		}
		if (DEBUG) cancelData << "Index: 12, Persistence: " << pp.persistence << " + " << pp.symPerturb1 << "e" << ", Cancellable: Yes.\n";
        return uniquePath;
    }
	if (DEBUG) cancelData << "Index: 12, Persistence: " << pp.persistence << " + " << pp.symPerturb1 << "e" << pp.symPerturb2 <<", Cancellable: No, Reason: No path.\n";
    return new vector<Simplex*>();
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
				//if (DEBUG) {e->output(); cout<<" ";}
			}else{
			    //if (DEBUG) {((Triangle*)s)->output(); cout<<" ";}
			}
		}
		//if (DEBUG) cout<<endl;
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

	if (DEBUG) this->outputPersistencePairs("persistencePairs.txt");

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
	int Pcounter = 0;
	while (i < msPersistencePairs.size() && j < smPersistencePairs.size()){
		persistencePair01 pair1 = msPersistencePairs[i];
		persistencePair12 pair2 = smPersistencePairs[j];

		if (pair1.persistence < pair2.persistence || (fabs(pair1.persistence - pair2.persistence) < 1e-8) && pair1.symPerturb < pair2.symPerturb1 + 1e-8){
			if (pair1.persistence < delta + 1e-8){
				//start = chrono::system_clock::now();
				vector<Simplex*> *VPath = this->isCancellable(pair1, cancelData);

				//end = chrono::system_clock::now();
				//chrono::duration<double> eDur = end - start;
				//double elapsedTime = eDur.count();
				//average4 = (elapsedTime + a4Count * average4) / (a4Count + 1);
				//a4Count++;
				// if (DEBUG) cout<< pair1.persistence << endl;
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
			if (pair2.persistence < delta + 1e-8){
				//start = chrono::system_clock::now();
				vector<Simplex*> *VPath = this->isCancellable(pair2, cancelData);
				/*end = chrono::system_clock::now();
				chrono::duration<double> eDur = end - start;
				double elapsedTime = eDur.count();
				average4 = (elapsedTime + a4Count * average4) / (a4Count + 1);
				a4Count++;*/
				// if (DEBUG) cout<< pair2.persistence << endl;
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
		if (Pcounter%10000==0){
			cout << "\r";
			cout << "\t" << Pcounter << "/" << msPersistencePairs.size() + smPersistencePairs.size() << "...";
			cout.flush();
		}
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
			if (i % 100 == 0){
				cout << "\r";
				cout << "\t" << i << "/" << msPersistencePairs.size() <<"...";
				cout.flush();
			}
		}
		cout<<endl;
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
			if (j % 100 == 0){
				cout << "\r";
				cout << "\t" << j << "/" << smPersistencePairs.size();
				cout.flush();
			}
		}
		cout<<endl;
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
