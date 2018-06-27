#include<vector>
#include<string>
#include<fstream>
#include<iostream>
#include<string.h>
using namespace std;

struct parameters{
	string ininame;
	string vertfile;
	string edgefile;
	string outputprefix;
	double thd;
	int comp=0;
	vector<int> shift;
	vector<int> root;
};

void default_ini();
string getnext(ifstream & filein, int & stat);
vector<int> tokenize(const string &line, int & stat);
int loadini(const string &filename, parameters &p);
string getpath(const string &filename);