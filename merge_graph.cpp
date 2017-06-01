/*
Merge result
Author: Suyi Wang

Input: xxx_vert.txt xxx_edge.txt
Output: binary file

Comments: Vertex index start from 1. All edges and triangles uses vertex index.
Output segments from DM are lines connecting neighbouring cells along axis
so diffusing is simplified.
*/


// g++ merge_graph.cpp -std=c++11 -I /home/wangsu/Develop/Neuro/DiademMetric/persistence/boost_1_59_0/ -o merge_graph_exe


#include<stdio.h>
#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<algorithm>
#include<cmath>


#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>

using namespace std;


#define DEBUG 0

// Do not fill inner part of a cube - 12
// Do     fill inner part of a cube - 16
int nb = 12;


// *********** begin vertex pair***********
struct point{
    int x,y,z;
    double v;
    bool in_bbox;
};

std::size_t hash_value(const point &pt){
  std::size_t seed = 0;
  boost::hash_combine(seed, pt.x);
  boost::hash_combine(seed, pt.y);
  boost::hash_combine(seed, pt.z);
  return seed;
}

bool operator==(const point &a, const point &b)
{
  return a.x == b.x && a.y == b.y && a.z == b.z;
}

class VertexHash{
	private:
		boost::unordered_map<point, int> v_hash;
	public:
		VertexHash(){
			v_hash.clear();
		}
		int GetIndex(point p){
			if (v_hash.count(p) > 0)
				return v_hash[p];
			else
				return -1;
		}
		void InsertVertex(point p, int n){
			v_hash.insert(make_pair(p, n));
		}
		int size(){return v_hash.size();}
		boost::unordered_map<point, int>::iterator begin(){
			return v_hash.begin();
		}
		boost::unordered_map<point, int>::iterator end(){
			return v_hash.end();
		}
};
// *********** End of vertex pair ***********


// *********** begin Edge pair & it's hashing***********
class cp{
	public:
		int p1,p2;  // must be vertices after merge.
		void Reorder(){
			if (p1>p2) swap(p1, p2);
		}
};

std::size_t hash_value(const cp &e){
  std::size_t seed = 0;
  boost::hash_combine(seed, e.p1);
  boost::hash_combine(seed, e.p2);
  return seed;
}

bool operator==(const cp &a, const cp &b)
{
  return a.p1 == b.p1 && a.p2 == b.p2;
}

class EdgeHash{
	// Assuming vertex are given in ascending order
	private:
		boost::unordered_set<cp> e_hash;
	
	public:
		EdgeHash(){
			e_hash.clear();
		}
		bool HasEdge(cp edge){
			return e_hash.count(edge) > 0;
		}
		void InsertEdge(cp edge){
			e_hash.insert(edge);
		}
		int size(){return e_hash.size();}
		boost::unordered_set<cp>::iterator begin(){return e_hash.begin();}
		boost::unordered_set<cp>::iterator end(){return e_hash.end();}
};
// *********** End of Edge pair***********


// *********** begin Triangle pair***********
struct tp{
	public:
		int p1,p2,p3;
		void Reorder(){
			if (p1>p2) swap(p1, p2);
			if (p1>p3) swap(p1, p3);
			if (p2>p3) swap(p2, p3);
		}
};

std::size_t hash_value(const tp &t){
  std::size_t seed = 0;
  boost::hash_combine(seed, t.p1);
  boost::hash_combine(seed, t.p2);
  boost::hash_combine(seed, t.p3);
  return seed;
}

bool operator==(const tp &a, const tp &b)
{
  return a.p1 == b.p1 && a.p2 == b.p2 && a.p3 == b.p3;
}

class TriangleHash{
	// Assuming vertex are given in ascending order
	private:
		boost::unordered_set<tp> t_hash;
	public:
		TriangleHash(){
			t_hash.clear();
		}
		bool HasTriangle(tp triangle){
			return t_hash.count(triangle) > 0;
		}
		void InsertTriangle(tp triangle){
			t_hash.insert(triangle);
		}
		int size(){return t_hash.size();}
		boost::unordered_set<tp>::iterator begin(){return t_hash.begin();}
		boost::unordered_set<tp>::iterator end(){return t_hash.end();}
};
// *********** End of Triangle pair ***********


int vcount = 0; //for vector vertex
vector<point> vertex;
vector<cp> edge;
vector<tp> triangle;


VertexHash vh;
EdgeHash eh;
TriangleHash th;

vector<int> trans_info;
vector<point> graph_vert;


int triangle_cube(int i, int j, int k, int AB){
    if (AB == 0){ // triangles of Type A
        int TypeAtri[3*16][3] = {{0,0,0}, {1,0,1}, {0,0,1},
                          {0,0,0}, {1,0,0}, {1,0,1},
                          {0,0,0}, {0,1,1}, {0,1,0},
                          {0,0,0}, {0,0,1}, {0,1,1},
                          {0,0,0}, {1,1,0}, {1,0,0},
                          {0,0,0}, {0,1,0}, {1,1,0},
						  {1,0,0}, {1,1,0}, {1,0,1},
						  {1,0,1}, {1,1,0}, {1,1,1},
						  {0,1,0}, {0,1,1}, {1,1,0},
						  {0,1,1}, {1,1,0}, {1,1,1},
						  {0,0,1}, {0,1,1}, {1,0,1},
						  {0,1,1}, {1,1,1}, {1,0,1},
                          {0,0,0}, {0,1,1}, {1,0,1},// fill
                          {0,0,0}, {1,0,1}, {1,1,0},
                          {1,1,0}, {1,0,1}, {0,1,1},
                          {0,0,0}, {1,1,0}, {0,1,1}
                   };
        for (int cnt = 0; cnt < nb; cnt++){
            int sub1, sub2, sub3;
            point p;
            p.x = i + TypeAtri[cnt*3][0]; p.y = j + TypeAtri[cnt*3][1]; p.z = k + TypeAtri[cnt*3][2];
            sub1 = vh.GetIndex(p);
            if (sub1 < 0){
            	p.v = 1e-6;
                vh.InsertVertex(p, vcount);
                sub1 = vcount;
                vcount++;
                vertex.push_back(p);
            }
			
			p.x = i + TypeAtri[cnt*3+1][0]; p.y = j + TypeAtri[cnt*3+1][1];
			p.z = k + TypeAtri[cnt*3+1][2];
			sub2 = vh.GetIndex(p);
            if (sub2 < 0){
                p.v = 1e-6;
                vh.InsertVertex(p, vcount);
                sub2 = vcount;
              	vcount++;
                vertex.push_back(p);
            }

			p.x = i + TypeAtri[cnt*3+2][0]; p.y = j + TypeAtri[cnt*3+2][1];
			p.z = k + TypeAtri[cnt*3+2][2];
			sub3 = vh.GetIndex(p);
            if (sub3 < 0){
                p.v = 1e-6;
                vh.InsertVertex(p, vcount);
                sub3 = vcount;
                vcount++;
                vertex.push_back(p);
            }
            
            tp new_triangle;
            new_triangle.p1 = sub1; new_triangle.p2 = sub2; new_triangle.p3 = sub3;
			new_triangle.Reorder();
			if (!th.HasTriangle(new_triangle)){
				triangle.push_back(new_triangle);
				th.InsertTriangle(new_triangle);
			}

            cp new_edge1;
            new_edge1.p1 = sub1; new_edge1.p2 = sub2;
			new_edge1.Reorder();
            cp new_edge2;
            new_edge2.p1 = sub1; new_edge2.p2 = sub3;
			new_edge2.Reorder();
            cp new_edge3;
            new_edge3.p1 = sub2; new_edge3.p2 = sub3;
			new_edge3.Reorder();

            if (DEBUG){
                printf("%d %d %d\n", sub1, sub2, sub3);
            }

            if (!eh.HasEdge(new_edge1)){
				edge.push_back(new_edge1);
				eh.InsertEdge(new_edge1);
			}
			if (!eh.HasEdge(new_edge2)){
				edge.push_back(new_edge2);
				eh.InsertEdge(new_edge2);
			}
			if (!eh.HasEdge(new_edge3)){
				edge.push_back(new_edge3);
				eh.InsertEdge(new_edge3);
			}
        }
    }
    else{// triangles of Type B
        int TypeAtri[3*16][3] = {{0,0,0}, {0,0,1}, {0,1,0},
                               {0,0,0}, {0,1,0}, {1,0,0},
                               {0,0,0}, {1,0,0}, {0,0,1},
                               {0,1,1}, {0,1,0}, {0,0,1},
                               {0,1,0}, {1,1,0}, {1,0,0},
                               {0,0,1}, {1,0,0}, {1,0,1},
							   {0,1,1}, {1,1,1}, {0,0,1},
							   {0,0,1}, {1,0,1}, {1,1,1},
							   {1,0,1}, {1,1,1}, {1,0,0},
							   {1,0,0}, {1,1,0}, {1,1,1},
							   {0,1,1}, {1,1,1}, {0,1,0},
							   {0,1,0}, {1,1,1}, {1,1,0},
                               {0,0,1}, {0,1,0}, {1,0,0},// fill
                               {0,1,0}, {1,0,0}, {1,1,1},
                               {0,0,1}, {1,1,1}, {1,0,0},
                               {0,0,1}, {1,1,1}, {0,1,0}
                               };
		for (int cnt = 0; cnt < nb; cnt++){
            int sub1, sub2, sub3;
            point p;
            p.x = i + TypeAtri[cnt*3][0]; p.y = j + TypeAtri[cnt*3][1]; p.z = k + TypeAtri[cnt*3][2];
            sub1 = vh.GetIndex(p);
            if (sub1 < 0){
            	p.v = 1e-6;
                vh.InsertVertex(p, vcount);
                sub1 = vcount;
                vcount++;
                vertex.push_back(p);
            }
			
			p.x = i + TypeAtri[cnt*3+1][0]; p.y = j + TypeAtri[cnt*3+1][1];
			p.z = k + TypeAtri[cnt*3+1][2];
			sub2 = vh.GetIndex(p);
            if (sub2 < 0){
                p.v = 1e-6;
                vh.InsertVertex(p, vcount);
                sub2 = vcount;
              	vcount++;
                vertex.push_back(p);
            }

			p.x = i + TypeAtri[cnt*3+2][0]; p.y = j + TypeAtri[cnt*3+2][1];
			p.z = k + TypeAtri[cnt*3+2][2];
			sub3 = vh.GetIndex(p);
            if (sub3 < 0){
                p.v = 1e-6;
                vh.InsertVertex(p, vcount);
                sub3 = vcount;
                vcount++;
                vertex.push_back(p);
            }
            
            tp new_triangle;
            new_triangle.p1 = sub1; new_triangle.p2 = sub2; new_triangle.p3 = sub3;
			new_triangle.Reorder();
			if (!th.HasTriangle(new_triangle)){
				triangle.push_back(new_triangle);
				th.InsertTriangle(new_triangle);
			}

            cp new_edge1;
            new_edge1.p1 = sub1; new_edge1.p2 = sub2;
			new_edge1.Reorder();
            cp new_edge2;
            new_edge2.p1 = sub1; new_edge2.p2 = sub3;
			new_edge2.Reorder();
            cp new_edge3;
            new_edge3.p1 = sub2; new_edge3.p2 = sub3;
			new_edge3.Reorder();

            if (DEBUG){
                printf("%d %d %d\n", sub1, sub2, sub3);
            }

            if (!eh.HasEdge(new_edge1)){
				edge.push_back(new_edge1);
				eh.InsertEdge(new_edge1);
			}
			if (!eh.HasEdge(new_edge2)){
				edge.push_back(new_edge2);
				eh.InsertEdge(new_edge2);
			}
			if (!eh.HasEdge(new_edge3)){
				edge.push_back(new_edge3);
				eh.InsertEdge(new_edge3);
			}
        }
    }
    return 0;
}


double norm_const = 1;
int diffuse(point p){
	double sigma = 1;
	int sum = 0;
	for (int i = -1; i <= 1; ++i)
		for(int j = -1; j <= 1; ++j)
			for(int k = -1; k <= 1; ++k){
				point dif_p;
				
				dif_p.x = p.x + i;
				dif_p.y = p.y + j;
				dif_p.z = p.z + k;
				dif_p.v = p.v * norm_const * exp(- (i*i + j*j +k*k)/(2* sigma * sigma));
				dif_p.in_bbox = false;
				/*if (dif_p.v < 1e-8) {
					continue;
				}*/
				
				int idx = vh.GetIndex(dif_p);
				if (idx < 0){
					dif_p.in_bbox = false;
					vh.InsertVertex(dif_p, vcount);
					vertex.push_back(dif_p);
					vcount++;
				}else{
					vertex[idx].v += dif_p.v;
					vertex[idx].in_bbox = false;
				}
    }
    return sum;
}


bool in_range(point p, vector<int> &bbox){
	if (p.x >= bbox[0] && p.x <= bbox[1] &&
		p.y >= bbox[2] && p.y <= bbox[3] &&
		p.z >= bbox[4] && p.z <= bbox[5])
	{
	   return true;
	}
	else{
		return false;
	}
}


void ProcessGraph(string filename, vector<int> bbox){
	string vert_name = filename + "_vert.txt";
	string edge_name = filename + "_edge.txt";
	
	FILE* fp = fopen(vert_name.c_str(), "r");
	if (fp==NULL){
		cout << "Cannot open "<<vert_name << endl;
		return;
	}
	double x, y, z, v;
	int nl;
	graph_vert.clear();

	while(fscanf(fp, "%lf%lf%lf%lf%d", &x, &y, &z, &v, &nl) != EOF){
		point p;
		p.x = x; p.y = y; p.z = z; p.v = v;
		graph_vert.push_back(p);  // + 1 or not
	}
	fclose(fp);
	printf("Read %d vertices\n", graph_vert.size());
	
	
	fp = fopen(edge_name.c_str(), "r");
	int e1, e2;
	int counter = 0;
	int c_interior = 0;
	while(fscanf(fp, "%d%d%d", &e1, &e2, &nl) != EOF){
		e1--;e2--;
		if (in_range(graph_vert[e1], bbox) && in_range(graph_vert[e2], bbox)){
			// only insert a point.
			int idx = vh.GetIndex(graph_vert[e1]);
			int idx1 = -1, idx2 = -1;
			if (idx < 0){
				graph_vert[e1].in_bbox = true;
				vh.InsertVertex(graph_vert[e1], vcount);
				vertex.push_back(graph_vert[e1]);
				idx1 = vcount;
				vcount++;
			}else{
				vertex[idx].v += graph_vert[e1].v;
				idx1 = idx;
			}
			
			idx = vh.GetIndex(graph_vert[e2]);
			if (idx < 0){
				graph_vert[e2].in_bbox = true;
				vh.InsertVertex(graph_vert[e2], vcount);
				vertex.push_back(graph_vert[e2]);
				idx2 = vcount;
				vcount++;
			}else{
				vertex[idx].v += graph_vert[e2].v;
				idx2 = idx;
			}

			cp new_edge;
            new_edge.p1 = idx1; new_edge.p2 = idx2;
			new_edge.Reorder();
			if (!eh.HasEdge(new_edge)){
				edge.push_back(new_edge);
				eh.InsertEdge(new_edge);
			}
			c_interior++;
		}	
		else{
			if (!in_range(graph_vert[e1], bbox)){
				diffuse(graph_vert[e1]);
			}
			if (!in_range(graph_vert[e2], bbox)){
				diffuse(graph_vert[e2]);
			}
			counter ++;
		}
	}
	fclose(fp);
	printf("processed %d interior edges, %d diffused points\n", c_interior, counter, vh.size());
	
    printf("done\n");
}


void bin_output(){
    ofstream ofs("vert.bin",ios::binary);
    printf("writing %d vertex\n", vertex.size());
    char* vert = new char[sizeof(double) * 4];
    double* vert_buffer = (double*) vert;
    for (int i = 0; i < vertex.size(); i++){
    	vert_buffer[0] = vertex[i].x; vert_buffer[1] = vertex[i].y;
    	vert_buffer[2] = vertex[i].z; vert_buffer[3] = vertex[i].v;
		ofs.write(vert, sizeof(double) * 4);
    }
    ofs.close();

	ofs.open("edge.bin", ios::binary);
    printf("writing %d edge\n", edge.size());
    char* edgechar = new char[sizeof(int) * 2];
    int* edge_buffer = (int*) edgechar;
    for (int i = 0; i < edge.size(); i++){
    	edge_buffer[0] = edge[i].p1; edge_buffer[1] = edge[i].p2;
        ofs.write(edgechar, sizeof(int) * 2);
    }
    ofs.close();

    ofs.open("triangle.bin", ios::binary);
    printf("writing %d triangle\n", triangle.size());
    char* trianglechar = new char[sizeof(int) * 3];
    int* triangle_buffer = (int*) trianglechar;
    for (int i = 0; i < triangle.size(); i++){
    	triangle_buffer[0] = triangle[i].p1;
    	triangle_buffer[1] = triangle[i].p2;
    	triangle_buffer[2] = triangle[i].p3;
        ofs.write(trianglechar, sizeof(int) * 3);
    }
    ofs.close();
}


double kernel_init(double sigma){
	double sum = 0;
	for (int i = -1; i <= 1; ++i)
		for(int j = -1; j <= 1; ++j){
			for(int k = -1; k <= 1; ++k){
				sum += exp(- (i*i + j*j +k*k)/(2* sigma * sigma));
				// cout << exp(- (i*i + j*j +k*k)/(2* sigma * sigma)) << " ";
			}
		// cout << "\n";
    }
    return 1.0/sum;
}


int Triangulate(){
	int original_total = vertex.size();
	
	double THD = 1e-6;
	int skip_count = 0;
    for(int v = 0; v < original_total; ++v){
		int i, j, k, val;
		
		if (v%10000==0){
			cout << '\r';
			cout << v << '/' << original_total;
			cout.flush();
		}
		
		i = vertex[v].x; j = vertex[v].y; k = vertex[v].z; val = vertex[v].v;
		if (val < THD ||vertex[v].in_bbox) {
			// cout << "skipped something\n";
			skip_count++;
			continue;
		}

		if ((i+j+k)%2==1){
			triangle_cube(i, j, k, 0);
		}
		else{
			triangle_cube(i, j, k, 1);
		}
	}
    
    printf("Skipped %d points\n", skip_count);
    return 0;
}


void init_trans(string trans_file){
	FILE* fp = fopen(trans_file.c_str(), "r");
	int xlen;
	while(fscanf(fp, "%d", &xlen)!=EOF){
		trans_info.push_back(xlen);
	}
	// should be dividable by 4
	fclose(fp);
}


int main(int argc, char* argv[])
{
	// input contains prefix S for dataset & number # for max number.
	// format: s_#_vert.txt
	// merge_graph_exe <prefix> <trans_matrix> <max_file_num>
	string prefix;
	string trans;
	int max_file;
	if (argc == 4){
		prefix = string(argv[1]);
		trans = string(argv[2]);
		max_file = atoi(argv[3]);
		nb = 12;
	}
	else {
		cout << "usage: merge_graph <dataset_name> <trans_matrix> <max_filecount>\n";
		return 0;
	}
	
	trans_info.clear();
	init_trans(trans);
	
	vertex.clear(); edge.clear(); triangle.clear();
	norm_const = kernel_init(0.5);
	cout << "Normalize factor: " << norm_const << endl;

	vector<int> bbox;
	for(int i = 0; i < max_file; i++){
		string filename = prefix + to_string(i);
		printf("Processing Graph %s...\n", filename.c_str());
		
		bbox.clear();
		int lz;
		int j = i-1;
		lz = trans_info[j*4 + 3];
		

		// for NeuroMuscular and Neocortical
		// This part is still hard coded.
		
		int xmin = trans_info[j*4] + 512/50.0;		bbox.push_back(xmin);
		int xmax = trans_info[j*4] + 512*0.98;		bbox.push_back(xmax);
		int ymin = trans_info[j*4 + 1] + 512/50.0;	bbox.push_back(ymin);
		int ymax = trans_info[j*4 + 1] + 512*0.98;	bbox.push_back(ymax);
		int zmin = trans_info[j*4 + 2] + 0;		bbox.push_back(zmin);
		int zmax = trans_info[j*4 + 2] + 300;		bbox.push_back(zmax);
		

		// for OP
		/*
		int xmin = trans_info[j*4] + 5;		bbox.push_back(xmin);
		int xmax = trans_info[j*4] + 261 - 5;		bbox.push_back(xmax);
		int ymin = trans_info[j*4 + 1] + 5;	bbox.push_back(ymin);
		int ymax = trans_info[j*4 + 1] + 261 - 5;	bbox.push_back(ymax);
		int zmin = trans_info[j*4 + 2] + 0;		bbox.push_back(zmin);
		int zmax = trans_info[j*4 + 2] + 200;		bbox.push_back(zmax);
		*/

		
			for(auto &x : bbox){
				cout << x << " ";
			}
			cout << endl;
		
		ProcessGraph(filename, bbox);
	}
	
	Triangulate();
	
	bin_output();

    printf("Done\n");
    return 0;
}


