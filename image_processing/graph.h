#include<vector>
#include<string>
#include<cstdio>
#include<iostream>
#include<algorithm>
#include<unordered_map>
#include<fstream>
#include<cmath>
#include<limits>
#include<queue>
using namespace std;

struct point{
	// We are processing pixels so the position is always integer
	vector<int> pos;
	double f;

	point(int x, int y, int z, int func){
		pos.push_back(x);pos.push_back(y);pos.push_back(z);
		f = func;
	}
	point(const point &p){
		pos = p.pos;
		f = p.f;
	}
	point(const vector<int> &p){
		pos = p;
		f = 0;
	}
	point(){return;}

	// Maybe there is better solution
	friend ostream& operator<<(ostream& os, const point & p)  
	{
		for(auto x : p.pos)
			os << x << ' ';
		os << p.f << " 0";
		return os;
	}
};


struct vertex{
	int idx;
	int prev;
	double f;

	vertex(int id, int from_idx, double dist){
		idx = id; f = dist; prev = from_idx;
	}
	vertex(){return;}
};


class vertex_cmp
{
public:
    bool operator() (const vertex &a, const vertex &b)
    {
    	// Min Heap
        return a.f > b.f;
    }
};


// C++ template to print pair<>
// class by using template
template <typename T, typename S>
ostream& operator<<(ostream& os, const pair<T, S>& v)
{
    os << v.first << " "
       << v.second;
    return os;
}

/*
ostream& operator<<(ostream& os, const pair<int, int>& edge)
{
    os << edge.first << " " << edge.second << "-1 1";
    return os;
}
*/


class graph{
private:
	vector<point> v;
	vector<vector<int> > e;
	int dist_created=0;
	
public:
	graph();
	graph(string filename);
	graph(string vert, string edge);
	
	int fromfile(const string & filename);
	int loadvert(const string & filename);
	int loadedge(const string & filename);

	int size();
	
	int threshold_saddle(double thd);
	
	int shiftvert(vector<double> pos);
	int find_vert(vector<int> pos);
	double get_dist(point a, point b);
	
	int check_redundancy();
	
	graph traverse(const int start, const vector<point> & v, 
					const vector<vector<int> > & e, vector<bool> & touched);
	vector<pair<int, int> > components();
	
	vector<vector<int> > dijkstra(const vector<int> &pos);
	
	vector<graph> split(int n);

	int set_edge(const vector<vector<int> > & tree_edge);
	int to_swc(string filename);
	int to_file(string filename);
	int to_simplex(string filename);
};