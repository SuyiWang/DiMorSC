class DiscreteVField;
class DiscreteVField{
private:
	// use original index.
	unordered_map<int, int> VEmap;
	unordered_map<int, int> ETmap;
	
	// gonna deprecate the following
	// unordered_set<Edge*> *Eflag;
	// unordered_map<Simplex*, vector<Simplex*>* > *jumpmap;

public:
	DiscreteVField();
	int containsVE(int v);
	int containsET(int e);
	
	void addVE(int v, int e);
	void addET(int e, int t);
	
	void removeVE(int v, int e);
	void removeET(int e, int t);
	// deprecated
	/*
	void outputVEmap();
	void outputETmap();
	
	bool containsEdge(Edge *e);
	void UnWarp(Simplex* s, bool flip, vector<Simplex*>* path);
	bool CanJump(Simplex*);
	void addJump(Simplex*, vector<Simplex*>*);
	*/
};

DiscreteVField::DiscreteVField(){
	VEmap.clear();
	ETmap.clear();
}

int DiscreteVField::containsVE(int v){
	if (VEmap.count(v) > 0){
		return VEmap.at(v);
	}
	return -1;
}

int DiscreteVField::containsET(int e){
	if (ETmap.count(e) > 0){
		return ETmap.at(e);
	}
	return -1;
}

void DiscreteVField::addVE(int v, int e){
	std::pair<int, int> pair = std::make_pair(v, e);
	VEmap.insert(pair);
}

void DiscreteVField::addET(int e, int t){
	std::pair<int, int> pair = std::make_pair(e, t);
	ETmap.insert(pair);
}

void DiscreteVField::removeVE(int v, int e){
	if (containsVE(v) != e){
		cerr << "Removing unexisting pair" <<endl;
	}
	VEmap.erase(v);
}

void DiscreteVField::removeET(int e, int t){
	if (containsET(e) != t){
		cerr << "Removing unexisting pair" <<endl;
	}
	ETmap.erase(e);
}
