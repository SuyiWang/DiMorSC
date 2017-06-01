class Simplex;
class Triangle;
class Edge;
class Vertex;

class Simplex{
public:
	/*Simplex may be cast to Vertex, Edge, or Triangle based on whether
	dim = 0, 1, or 2, respectively*/
	// It contains only local info: A Simplex knows only its vertices
	// Interconnection between Simplex is in Simplicial Complex
	int dim;
	double funcValue;
	unsigned int filtrationPosition;
};


class Vertex:public Simplex{
	// will deprecate coords, since it is uncessary to have them.
	double coords[MAX_DIM];
	int vPosition;
	int oriPosition;
	// vector<Edge*> incidenceList;
	// moved to Class simplicialcomplex
	
	
public:
	Vertex(double *coords, double funcValue);
	
	double* getCoords();
	double getFuncValue();
	void setFuncValue(double v){
		funcValue = v;
	}
	int getVPosition();
	void setVposition(int p);
	int getoriPosition(){
		return oriPosition;
	}
	void setoriposition(int p){
		oriPosition = p;
	}
	void output(ofstream& ofs){
		for(int i = 0; i< DIM; i++){
			ofs << coords[i] <<" ";
		}
	    ofs << funcValue << " " << vPosition << "/" << oriPosition << " " << dim <<"\n";
	}
	// deprecated functions
};


class Edge:public Simplex{
	// all use original position
	int vertices[2];		// vertices should be sorted
	int ePosition;
	Vertex* vp[2];
	double eval = 0;		// used to mark supporting saddle
	// moved to Class Simplicialcomplex
	// vector<Triangle*> incidenceList;

public:
	double persistence;
	int critical_type;
	Edge(int v1, int v2){
		vertices[0] = v1; vertices[1] = v2;
		dim = 1;
	}
	
	int* getVertices();
	double getFuncValue();
	void setFuncValue(double funcValue);
	int getEPosition();
	void setEposition(int p);
	void set_vp(Vertex** v){
		vp[0] = v[0]; vp[1] = v[1];
	}
	Vertex** get_vp(){
		return vp;
	}
	void output(ofstream& ofs){
	    ofs << vertices[0] << "~" << vertices[1] << ":" << ePosition << "@" << dim <<"\n";
	}
	double getEval() {return eval;}
	void setEval(double v) {eval = v;}
	// Deprecated functions
};


class Triangle:public Simplex{
	int edges[3];
	int vertices[3];		// vertices should be sorted
	int tPosition;
	Vertex* vp[3];

public:
	Triangle(vector<int> v, vector<int> e){
		for(int i = 0; i < 3; i++){
			edges[i] = e[i];
			vertices[i] = v[i];
		}
		dim = 2;
	}
	
	int* getEdges();
	int* getVertices();
	double getFuncValue();
	void setFuncValue(double funcValue);
	int getTPosition();
	void setTposition(int p);
	void set_vp(Vertex** v){
		vp[0] = v[0]; vp[1] = v[1]; vp[2] = v[2];
	}
	Vertex** get_vp(){
		return vp;
	}
	
	void output(ofstream& ofs){
		for(int i = 0; i < 3; i++){
			ofs << vertices[i] << " ";
		}
	    for(int i = 0; i < 3; i++){
			ofs << edges[i] << " ";
		}
		ofs << tPosition << "@"<< dim << endl;
		ofs <<endl;
	}
	
	// Deprecated
};

Vertex::Vertex(double *coords, double funcValue){
	for (int i = 0; i < DIM; i++) {
		this->coords[i] = coords[i];
	}
	this->funcValue = funcValue;
	this->dim = 0;
}

double* Vertex::getCoords(){
	return coords;
}

double Vertex::getFuncValue(){
	return funcValue;
}

int Vertex::getVPosition(){
	return vPosition;
}

void Vertex::setVposition(int p){
	vPosition = p;
}

int* Edge::getVertices(){
	return vertices;
}

double Edge::getFuncValue(){
	return funcValue;
}

void Edge::setFuncValue(double v){
	funcValue = v;
}

int Edge::getEPosition(){
	return ePosition;
}

void Edge::setEposition(int p){
	ePosition = p;
}

int* Triangle::getEdges(){
	return edges;
}

int* Triangle::getVertices(){
	return vertices;
}

double Triangle::getFuncValue(){
	return funcValue;
}
void Triangle::setFuncValue(double v){
	funcValue = v;
}
int Triangle::getTPosition(){
	return tPosition;
}
void Triangle::setTposition(int p){
	tPosition = p;
}
