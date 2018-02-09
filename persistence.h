#include <algorithm>
#include <list>
#include <stack>

// PHAT dependencies
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


//  This stores a vertex edge pair
typedef struct {
	// here, use index in sorted_vertex!!!
	// because we need the sorted vertex info for persistence pair sorting.
	int min;
	int saddle;
	double persistence;
	int loc_diff;
}persistencePair01;



//  This stores an edge triangle pair
typedef struct {
	int saddle;
	int max;
	double persistence;
	int loc_diff;
}persistencePair12;


//  Container of all pairs
class PersistencePairs{
	/*ms: min-saddle or 0-1
	  sm: saddle-max or 1-2*/
	vector<persistencePair01> msPersistencePairs;
	// we don't sort smPersitencePairs now, so no need for pointers.
	vector<persistencePair12> smPersistencePairs;

public:
	PersistencePairs(){
		msPersistencePairs.clear();
		smPersistencePairs.clear();
	}
	
	
	//  Compares v-e pair
	static bool persistencePairCompare01(const persistencePair01& p, const persistencePair01& q){
		if(p.persistence < q.persistence - EPS_compare){
			return true;
		}
		else if (p.persistence > q.persistence + EPS_compare){
			return false;
		}else{
			if (p.loc_diff < q.loc_diff)
				return true;
			else if (p.loc_diff > q.loc_diff)
				return false;
			else {
				return p.min < q.min;
			}
		}
	}
	
	//  Compares e-t pair
	static bool nzf_compare(const persistencePair12 &a, const persistencePair12 &b){
		if(a.persistence < b.persistence){
			return true;
		}else return false;
	}
	
	
	//  data r/w interface
	static persistencePair01 read_ve_pair(ifstream &presave);
	static persistencePair12 read_et_pair(ifstream &presave);
	static void write_ve_pair(persistencePair01, ofstream& presave);
	static void write_et_pair(persistencePair12, ofstream& presave);
	static void write_ve_pair_debug(persistencePair01, int, ofstream&);
	static void write_et_pair_debug(persistencePair12, ofstream&);
	
	
	//  iterators
	vector<persistencePair01>::iterator msBegin(){
		return msPersistencePairs.begin();
	}
	vector<persistencePair01>::iterator msEnd(){
		return msPersistencePairs.end();
	}
	vector<persistencePair12>::iterator smBegin(){
		return smPersistencePairs.begin();
	}
	vector<persistencePair12>::iterator smEnd(){
		return smPersistencePairs.end();
	}
	
	void msinsert(persistencePair01 ms){
		msPersistencePairs.push_back(ms);
	}
	void sminsert(persistencePair12 sm){
		smPersistencePairs.push_back(sm);
	}
	void sortmspair(){
		sort(msPersistencePairs.begin(), msPersistencePairs.end(), persistencePairCompare01);
	}
	
	int mssize(){
		return msPersistencePairs.size();
	}
	int smsize(){
		return smPersistencePairs.size();
	}
	
	void output_sm_pair(ofstream &et_stream){
		vector<persistencePair12> nz_sm_pair;
		nz_sm_pair.clear();
		for(auto pp = smBegin(); pp!=smEnd(); ++pp){
			if (pp->persistence > EPS_compare)
				nz_sm_pair.push_back(*pp);
		}
		cout << "\t" << nz_sm_pair.size() << " non-zero persistence pairs\n";
		sort(nz_sm_pair.begin(), nz_sm_pair.end(), nzf_compare);
		for(auto pp = nz_sm_pair.begin(); pp != nz_sm_pair.end(); ++pp){
			et_stream << pp->persistence << endl;
		}
	}
	// deprecated
	
};

persistencePair01 PersistencePairs::read_ve_pair(ifstream &presave){
	char* int_buffer = new char[sizeof(int)];
	int* int_reader = (int*) int_buffer;
	char* double_buffer = new char[sizeof(double)];
	double* double_reader = (double*) double_buffer;
	
	persistencePair01 pp;
	presave.read(int_buffer, sizeof(int));
	pp.min = *int_reader;
	presave.read(int_buffer, sizeof(int));
	pp.saddle = *int_reader;
	presave.read(double_buffer, sizeof(double));
	pp.persistence = *double_reader;
	presave.read(int_buffer, sizeof(int));
	pp.loc_diff = *int_reader;
	
	delete int_buffer;
	delete double_buffer;
	return pp;
}
persistencePair12 PersistencePairs::read_et_pair(ifstream &presave){
	char* int_buffer = new char[sizeof(int)];
	int* int_reader = (int*) int_buffer;
	char* double_buffer = new char[sizeof(double)];
	double* double_reader = (double*) double_buffer;
	
	persistencePair12 pp;
	presave.read(int_buffer, sizeof(int));
	pp.saddle = *int_reader;
	presave.read(int_buffer, sizeof(int));
	pp.max = *int_reader;
	presave.read(double_buffer, sizeof(double));
	pp.persistence = *double_reader;
	presave.read(int_buffer, sizeof(int));
	pp.loc_diff = *int_reader;
	
	delete int_buffer;
	delete double_buffer;
	return pp;
}
void PersistencePairs::write_ve_pair(persistencePair01 pp, ofstream& presave){
	char* int_buffer = new char[sizeof(int)];
	int* int_writer = (int*) int_buffer;
	char* double_buffer = new char[sizeof(double)];
	double* double_writer = (double*) double_buffer;
	
	*int_writer = pp.min;
	presave.write(int_buffer, sizeof(int));
	*int_writer = pp.saddle;
	presave.write(int_buffer, sizeof(int));
	*double_writer = pp.persistence;
	presave.write(double_buffer, sizeof(double));
	*int_writer = pp.loc_diff;
	presave.write(int_buffer, sizeof(int));
	
	delete int_buffer;
	delete double_buffer;
}
void PersistencePairs::write_et_pair(persistencePair12 pp, ofstream& presave){
	char* int_buffer = new char[sizeof(int)];
	int* int_writer = (int*) int_buffer;
	char* double_buffer = new char[sizeof(double)];
	double* double_writer = (double*) double_buffer;
	
	*int_writer = pp.saddle;
	presave.write(int_buffer, sizeof(int));
	*int_writer = pp.max;
	presave.write(int_buffer, sizeof(int));
	*double_writer = pp.persistence;
	presave.write(double_buffer, sizeof(double));
	*int_writer = pp.loc_diff;
	presave.write(int_buffer, sizeof(int));
	
	delete int_buffer;
	delete double_buffer;
}

void PersistencePairs::write_ve_pair_debug(persistencePair01 pp, int ori, ofstream& ppair){
	ppair << "1 " << pp.persistence << " " 
		  << ori << " "
		  << pp.saddle << "\n";
}

void PersistencePairs::write_et_pair_debug(persistencePair12 pp, ofstream& ppair){
	ppair << "2 " << pp.persistence << " "
		  << pp.saddle << " "
		  << pp.max << "\n";
}
