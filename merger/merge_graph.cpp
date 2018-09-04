/*
Merge multiple graph in global coordinate.
Author: Suyi Wang

Input: merger_config
Output: Simplicial complex (to be pipelined to DiMorSC)
*/

// g++ merger/merge_graph.cpp core/hash.cpp core/readini.cpp -O3 -std=c++11 -I./extern/boost -I./core -o bin/merge_graph -w 2>error

/*
Merger config has the following format:
Actual config should not contain comments

# output filename
merged
# overlap for x,y,z dimension
3 3 3
# filename bounding box (including the overlapping region)
0 0 0 0 71 50 50
1 0 0 47 71 50 97
2 0 0 94 71 50 117
3 0 47 0 71 97 50
4 0 47 47 71 97 97
5 0 47 94 71 97 117
6 0 94 0 71 101 50
7 0 94 47 71 101 97
8 0 94 94 71 101 117
*/


#include<stdio.h>
#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<algorithm>
#include<cmath>

#include"hash.h"
#include"readini.h"


using namespace std;


#define DEBUG 0

// Do not fill inner part of a cube - 12
// Do     fill inner part of a cube - 16
int nb = 12;

int vcount = 0; //for vector vertex
vector<point> vertex;
vector<cp> edge;
vector<tp> triangle;

VertexHash vh;
EdgeHash eh;
TriangleHash th;

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
			p.x = i + TypeAtri[cnt*3][0];
			p.y = j + TypeAtri[cnt*3][1];
			p.z = k + TypeAtri[cnt*3][2];
			sub1 = vh.GetIndex(p);
			if (sub1 < 0){
				p.v = 1e-6;
				vh.InsertVertex(p, vcount);
				sub1 = vcount;
				vcount++;
				vertex.push_back(p);
			}
			
			p.x = i + TypeAtri[cnt*3+1][0];
			p.y = j + TypeAtri[cnt*3+1][1];
			p.z = k + TypeAtri[cnt*3+1][2];
			sub2 = vh.GetIndex(p);
			if (sub2 < 0){
				p.v = 1e-6;
				vh.InsertVertex(p, vcount);
				sub2 = vcount;
				vcount++;
				vertex.push_back(p);
			}

			p.x = i + TypeAtri[cnt*3+2][0];
			p.y = j + TypeAtri[cnt*3+2][1];
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
				printf("\t%d %d %d\n", sub1, sub2, sub3);
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
			p.x = i + TypeAtri[cnt*3][0];
			p.y = j + TypeAtri[cnt*3][1];
			p.z = k + TypeAtri[cnt*3][2];
			sub1 = vh.GetIndex(p);
			if (sub1 < 0){
				p.v = 1e-6;
				vh.InsertVertex(p, vcount);
				sub1 = vcount;
				vcount++;
				vertex.push_back(p);
			}
			
			p.x = i + TypeAtri[cnt*3+1][0];
			p.y = j + TypeAtri[cnt*3+1][1];
			p.z = k + TypeAtri[cnt*3+1][2];
			sub2 = vh.GetIndex(p);
			if (sub2 < 0){
				p.v = 1e-6;
				vh.InsertVertex(p, vcount);
				sub2 = vcount;
				vcount++;
				vertex.push_back(p);
			}

			p.x = i + TypeAtri[cnt*3+2][0];
			p.y = j + TypeAtri[cnt*3+2][1];
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
				printf("\t%d %d %d\n", sub1, sub2, sub3);
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


bool in_range(point p, const vector<int> &bbox){
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

void update_bbox(vector<double> &vb, bool &first, point p){
	if (first){
		first = 0;
		vb[0] = p.x; vb[1] = p.x;
		vb[2] = p.y; vb[3] = p.y;
		vb[4] = p.z; vb[5] = p.z;
	}else{
		vb[0] = p.x < vb[0]?p.x:vb[0]; vb[1] = p.x > vb[1]?p.x:vb[1];
				vb[2] = p.y < vb[2]?p.y:vb[2]; vb[3] = p.y > vb[3]?p.y:vb[3];
				vb[4] = p.z < vb[4]?p.z:vb[4]; vb[5] = p.z > vb[5]?p.z:vb[5];
	}
}


void simplex_output(string fname){
    string binname = fname + ".sc";
    ofstream ofs(binname,ios::binary);
    printf("writing vertex\n");

    char* intwriter = new char[sizeof(int)];
    int* intbuffer = (int*) intwriter;

    char* vert = new char[sizeof(double) * 4];
    double* vert_buffer = (double*) vert;

    intbuffer[0] = vertex.size();
    ofs.write(intwriter, sizeof(int));
    for (int i = 0; i < vertex.size(); i++){
        vert_buffer[0] = vertex[i].x; vert_buffer[1] = vertex[i].y;
        vert_buffer[2] = vertex[i].z; vert_buffer[3] = vertex[i].v;
        ofs.write(vert, sizeof(double) * 4);
    }
    
    printf("writing edge\n");
    char* edgechar = new char[sizeof(int) * 2];
    int* edge_buffer = (int*) edgechar;
    intbuffer[0] = edge.size();
    ofs.write(intwriter, sizeof(int));
    for (int i = 0; i < edge.size(); i++){
        edge_buffer[0] = edge[i].p1; edge_buffer[1] = edge[i].p2;
        ofs.write(edgechar, sizeof(int) * 2);
    }
    
    printf("writing triangle\n");
    char* trianglechar = new char[sizeof(int) * 3];
    int* triangle_buffer = (int*) trianglechar;
    intbuffer[0] = triangle.size();
    ofs.write(intwriter, sizeof(int));
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
		if (vertex[v].in_bbox) {
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
	
	printf("\tSkipped %d points\n", skip_count);
	return 0;
}

// For MinGW only
/*
std::string to_string(int i)
{
	std::stringstream ss;
	ss << i;
	return ss.str();
}*/

const int overlapdim = 3;
struct parameter{
	string out_prefix;
	int overlap[overlapdim];
};
struct fileinfo{
	string name;
	// pixel coordinate
	// 2 pairs of 3D coordinate defines the bounding box
	// offset: xmin ymin zmin xmax ymax zmax
	int offset[6];
};

ostream& operator<<(ostream& os, const parameter& p)  
{
	os << "[outputname] " << p.out_prefix << " || [overlap]";
	for(int i = 0; i < overlapdim; i++)
		 os << " " << p.overlap[i];
	return os;
}

ostream& operator<<(ostream& os, const vector<int> &vec)  
{
	os << "[innerbox]";
	for(int i = 0; i < vec.size(); i++)
		 os << " " << vec[i];
	return os;
}


int init_file(const string &inputname,
			  parameter &para,
			  vector<fileinfo> &blocks){
	fstream ifs(inputname.c_str(), ios::in);
	// get filename
	ifs >> para.out_prefix;
	
	// get overlap
	for(int i = 0; i < 3; i++){
		ifs >> para.overlap[i];
	}
	
	// get block info
	blocks.clear();
	while(!ifs.eof()){
		fileinfo block_file;
		ifs >> block_file.name;
		if (ifs.eof()) break;

		for(int i = 0; i < 6; i++)
			ifs >> block_file.offset[i];
		blocks.push_back(block_file);
	}
	ifs.close();
	return 0;
}

int adjust_bbox(vector<fileinfo> &blocks, parameter para){
	for(int i = 0; i < blocks.size(); ++i){
		// make bbox in order: xmin, xmax, ymin, ymax, zmin, zmax
		swap(blocks[i].offset[1], blocks[i].offset[3]);
		swap(blocks[i].offset[2], blocks[i].offset[4]);
		swap(blocks[i].offset[2], blocks[i].offset[3]);
		for(int j = 0; j < overlapdim; ++j){
			// shrink axis-min
			blocks[i].offset[j*2] += para.overlap[j];
			// shrink axis-max
			blocks[i].offset[j*2+1] -= para.overlap[j];
		}
	}
	return 0;
}


void ProcessGraph(fileinfo blk){
	string vert_name = blk.name + "_vert.txt";
	string edge_name = blk.name + "_edge.txt";
	
	fstream fp(vert_name.c_str(), ios::in);
	if (fp.fail()){
		cout << "Cannot open "<<vert_name << endl;
		return;
	}
	
	int nl;
	graph_vert.clear();
	bool first = 1;
	vector<double> vertexbound(6, 0);
	
	string input_str;
	getline(fp, input_str);
	while(!fp.eof()){
		point p;
		sscanf(input_str.c_str(), 
			   "%d%d%d%lf",
			   &p.x, &p.y, &p.z, &p.v);
		//cout << input_str << endl;
		//cout << vector<int>({p.x, p.y, p.z, p.v}) << endl;
		graph_vert.push_back(p);  // + 1 or not
		update_bbox(vertexbound, first, p);
		getline(fp, input_str);
	}
	fp.close();
	printf("\tRead %d vertices\n \tbounded in:", graph_vert.size());
	for(int i =0; i< 6; i++) printf(" %.0f ", vertexbound[i]);
	printf("\n");
	
	
	fp.open(edge_name.c_str(), ios::in);
	int e1, e2;
	int counter = 0;
	int c_interior = 0;
	double persist;

	vector<int> bbox(blk.offset, blk.offset+ 6);
	getline(fp, input_str);
	while(!fp.eof()){
		sscanf(input_str.c_str(), "%d%d%d%lf", &e1, &e2, &nl,&persist);
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
				counter ++;
			}
			if (!in_range(graph_vert[e2], bbox)){
				diffuse(graph_vert[e2]);
				counter ++;
			}
		}
		getline(fp, input_str);
	}
	
	printf("\tprocessed %d interior edges, %d diffused points\n", c_interior, counter, vh.size());
	fp.close();
	printf("\tdone\n");
}


int main(int argc, char* argv[])
{
	// input contains prefix S for dataset & number # for max number.
	// format: s_#_vert.txt
	// merge_graph_exe <prefix> <trans_matrix> <max_file_num>
	string search_path;
	string setting_name;
	string trans;
	bool degreed;
	vector<double> bboundary;
	vector<int> filelist;

	parameter para;
	vector<fileinfo> blocks;
	if (argc == 2){
		/*
		prefix = string(argv[1]);
		trans = string(argv[2]);
		setting_name = string(argv[3]);
		nb = 12;
		*/
		string configname(argv[1]);
		cout << "Reading parameters from: " << configname << endl;
		init_file(configname, para, blocks);
		cout << para << endl;
		cout << "block count: " << blocks.size() << endl;

		search_path = getpath(argv[1]);
		cout << "work folder: " << search_path << endl;
	}
	else {
		cout << "usage: merge_graph <config_file>";
		return 0;
	}

	// adjust bounding box for each block according to overlap
	adjust_bbox(blocks, para);

	vertex.clear(); edge.clear(); triangle.clear();
	norm_const = kernel_init(0.5);
	cout << "Normalize factor: " << norm_const << endl;
	
	for(int i = 0; i < blocks.size(); i++){
		string filename = search_path + blocks[i].name;
		blocks[i].name = filename;
		printf("Processing Graph %s...\n", filename.c_str());
		cout << "\t" << vector<int>(blocks[i].offset, blocks[i].offset+ 6) << endl;
		
		ProcessGraph(blocks[i]);
	}
	
	Triangulate();
	
	simplex_output(search_path+para.out_prefix);

	printf("Done\n");
	
	return 0;
}
