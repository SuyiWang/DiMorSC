#include <set>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iterator>
using namespace std;

class SparseZ2Matrix{
private:
	/*Each list<int> is a column in the matrix. The entries in the set are
	the indices of the rows of the matrix with nonzero entries. e.g. the column
	[000000010000000100000101]^T would have the set {7,15,21,23}*/
	vector<vector<unsigned int>*>* matrix;
	unsigned int numRows;

public:
	/*constructs new SparseZ2Matrix with numCols columns. The number of rows*/
	SparseZ2Matrix(unsigned int numRows,unsigned int numCols);

	/*Sets the entry at (rowIndex, colIndex) to value. value = true by default*/
	void set(unsigned int rowIndex, unsigned int colIndex);

	/*Returns the entry at (rowIndex, colIndex)*/
	bool get(unsigned int rowIndex, unsigned int colIndex);

	/*Returns the highest index of a row with nonzero entry in column colIndex
	If the column colIndex is all zero, it returns -1*/
	unsigned int getPivotRow(unsigned int colIndex);

	/*Returns true if column colIndex has all zero entries. false otherwise*/
	bool isZeroColumn(unsigned int colIndex);

	/*Replaces the column rightCol in the matrix with (leftCol + rightCol), where
	'+' is Z2 addition (i.e. xor)*/
	void add(unsigned int leftCol, unsigned int rightCol);

	void clearColumn(unsigned int colIndex);

	void output(string pathname);

	~SparseZ2Matrix(){
		for (vector<vector<unsigned int>*>::iterator it = this->matrix->begin(); it != this->matrix->end(); it++) {
			delete *it;
		}
		delete this->matrix;
	}
};

SparseZ2Matrix::SparseZ2Matrix(unsigned int numRows, unsigned int numCols){
	this->matrix = new vector<vector<unsigned int>*>();
	this->matrix->reserve(numCols);
	/*Initialize new lists in every entry of the vector*/
	for (unsigned int i = 0; i < numCols; i++){
		matrix->push_back(new vector<unsigned int>());
	}
	this->numRows = numRows;
	this->matrix->shrink_to_fit();
}

void SparseZ2Matrix::set(unsigned int rowIndex, unsigned int colIndex){
	vector<unsigned int> *column = this->matrix->at(colIndex);
	column->push_back(rowIndex);
}

bool SparseZ2Matrix::get(unsigned int rowIndex, unsigned int colIndex){
	vector<unsigned int> *column = this->matrix->at(colIndex);
	bool containsIndex = binary_search(column->begin(), column->end(), rowIndex);
	return containsIndex;
}

unsigned int SparseZ2Matrix::getPivotRow(unsigned int colIndex){
	vector<unsigned int> *column = this->matrix->at(colIndex);
	/*The values in a set are sorted, so end() returns and iterator to just past the greatest element*/
	return column->at(column->size() - 1);
}

bool SparseZ2Matrix::isZeroColumn(unsigned int colIndex){
	vector<unsigned int> *column = this->matrix->at(colIndex);
	return column->empty();
}

void SparseZ2Matrix::add(unsigned int leftCol, unsigned int rightCol){

	/*std::set<unsigned int> *sum = new std::set<unsigned int>(this->matrix->at(leftCol)->begin(), this->matrix->at(leftCol)->end());

	for (std::set<unsigned int>::iterator it = this->matrix->at(rightCol)->begin(); it != this->matrix->at(rightCol)->end(); it++){
		unsigned int entry = *it;
		if (sum->count(entry) == 1){
			sum->erase(entry);
		}
		else{
			sum->insert(entry);
		}
	}

	std::set<unsigned int> *rightColumnSet = this->matrix->at(rightCol);
	*rightColumnSet = *sum;
	delete sum;*/

	/*std::set<unsigned int> *rightColumnSet = this->matrix->at(rightCol);
	std::set<unsigned int> *leftColumnSet = this->matrix->at(leftCol);
	for (std::set<unsigned int>::iterator it = leftColumnSet->begin(); it != leftColumnSet->end(); it++) {
		unsigned int entry = *it;
		if (rightColumnSet->count(entry) == 1) {
			rightColumnSet->erase(entry);
		}
		else {
			rightColumnSet->insert(entry);
		}
	}*/

	/*std::set<unsigned int> *rightColumnSet = this->matrix->at(rightCol);
	std::set<unsigned int> *leftColumnSet = this->matrix->at(leftCol);
	std::set<unsigned int> *sum = new std::set<unsigned int>();
	insert_iterator<std::set<unsigned int>> sum_in(*sum, sum->begin());
	set_symmetric_difference(rightColumnSet->begin(), rightColumnSet->end(), leftColumnSet->begin(), leftColumnSet->end(), sum_in);
	std::swap(*rightColumnSet, *sum);
	delete sum;*/

	vector<unsigned int> *rightColumn = this->matrix->at(rightCol);
	vector<unsigned int> *leftColumn = this->matrix->at(leftCol);
	vector<unsigned int> *sum = new vector<unsigned int>();
	set_symmetric_difference(rightColumn->begin(), rightColumn->end(), leftColumn->begin(), leftColumn->end(), back_inserter(*sum));
	std::swap(*rightColumn, *sum);
	delete sum;

}

void SparseZ2Matrix::clearColumn(unsigned int colIndex){
	vector<unsigned int> *column = this->matrix->at(colIndex);
	column->clear();
}

void SparseZ2Matrix::output(string pathname){
	ofstream out(pathname);
	for (int i = 0; i < this->matrix->size(); i++){
		vector<unsigned int> *s = this->matrix->at(i);
		out << "column " << i << " "  << ":";
		for (vector<unsigned int>::iterator it = s->begin(); it != s->end(); it++){
			out << " " << *it;
		}
		out << "\n";
	}
}