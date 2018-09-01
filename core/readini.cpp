#include"readini.h"

void default_ini(){
	// if output exist
	// ofstream outstream("output/0.ini");
	/*
		outstream << "# vertex filename\n";
		outstream << "output/0_vert.txt\n";
		outstream << "# edge filename\n";
		outstream << "output/0_edge.txt\n";
		outstream << "# output prefix: program writes to output/0.swc\n";
		outstream << "output/0\n";
		outstream << "# root\n";
		outstream << "0 0 0\n";
		outstream << "# threshold for simplification\n";
		outstream << "5\n";
	*/
	// else create file in current directory
	ofstream outstream("output/0.ini");
	outstream << "# vertex filename\n";
	outstream << "0_vert.txt\n";
	outstream << "# edge filename\n";
	outstream << "0_edge.txt\n";
	outstream << "# output prefix: program writes to output/0.swc\n";
	outstream << "0\n";
	outstream << "# root\n";
	outstream << "0 0 0\n";
	outstream << "# threshold for simplification\n";
	outstream << "5\n";
}


string getnext(ifstream & filein, int & stat){
	string txt;
	getline(filein, txt);
	while(!filein.eof() && !filein.fail() && txt[0] == '#'){
		getline(filein, txt);
	}

	if (!filein.eof() && !filein.fail())
		stat = 0; 
	else
		stat = -1;

	if (stat != -1) cout << txt << endl;
	return txt;
}


vector<int> tokenize(const string &line, int & stat){
	/*
	// work with cstring
	pch = strtok (str," ,.-");
	while (pch != NULL){
		printf ("%s\n",pch);
		pch = strtok (NULL, " ,.-");
	}*/
	vector<int> rtn;

    // Skip delimiters at beginning.
    string::size_type lastPos = line.find_first_not_of(" ", 0);
    // Find first "non-delimiter".
    string::size_type pos     = line.find_first_of(" ", lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        rtn.push_back(stoi(line.substr(lastPos, pos - lastPos)));
        // Skip delimiters.  Note the "not_of"
        lastPos = line.find_first_not_of(" ", pos);
        // Find next "non-delimiter"
        pos = line.find_first_of(" ", lastPos);
    }
    return rtn;
}


string getpath(const string &filename){
	int pathstop = filename.find_last_of("/\\");
	string path = filename.substr(0, pathstop);
	if (path.size() > 0) path = path + "/";
	return path;
}


int loadini(const string &filename, parameters &p){
	ifstream filein(filename);
	

	int stat = 0;
	string path = getpath(filename);

	p.vertfile = path + getnext(filein, stat);
	if (stat == -1) {cout << "cannot get vert file\n"; return -1;}
	
	p.edgefile = path + getnext(filein, stat);
	if (stat == -1) {cout << "cannot get edge file\n"; return -1;}
	
	p.outputprefix = getnext(filein, stat);
	if (stat == -1) {cout << "cannot get output path\n"; return -1;}
	
	string rootstr = getnext(filein, stat);
	if (stat == -1) {cout << "cannot read root value\n"; return -1;}
	p.root = tokenize(rootstr, stat);
	if (stat == -1) {cout << "cannot parse root value\n"; return -1;}

	p.thd = atof(getnext(filein, stat).c_str());
	if (stat == -1) {cout << "cannot get threshold\n"; return -1;}

	
	// optional
	int optstat = 0;
	p.comp = atoi(getnext(filein, optstat).c_str());
	if (optstat == -1) {
		// use only one component
		p.comp = -1;
	}

	string shiftstr = getnext(filein, optstat);
	if (optstat == -1) {cout << "shift [0 0 0]\n";}
	else{
		p.shift = tokenize(shiftstr, optstat);
		if (optstat == -1) {cout << "Cannot parse shift value\n"; return -1;}
	}

	return stat;
}