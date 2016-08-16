/*
Compute Triangulation on a 3D domain.
Author: Suyi Wang

Input: mapinput.txt
Output: vert.txt edge.txt triangle.txt

Comments: Vertex index start from 1. All edges and triangles uses vertex index.
*/

#include<stdio.h>
#include<fstream>
#include<vector>
#include<string>
#include<algorithm>
#include<unordered_map>
#include<unordered_set>
using namespace std;
#define DEBUG 0

// Do not fill inner part of a cube - 12
// Do     fill inner part of a cube - 16
const int nb = 12;

struct point{
    double x,y,z,v;
};

// Edge pair
class cp{
	public:
		int p1,p2;
		void Reorder(){
			if (p1>p2) swap(p1, p2);
		}
};


// Triangle pair
struct tp{
	public:
		int p1,p2,p3;
		void Reorder(){
			if (p1>p2) swap(p1, p2);
			if (p1>p3) swap(p1, p3);
			if (p2>p3) swap(p2, p3);
		}
};


//	Tetrahedron pair: if need this, uncomment all lines with tetra
struct tetra{
	vector<int> p{-1, -1, -1, -1};
};

vector<point> vertex;
vector<cp> edge;
vector<tp> triangle;
// vector<tetra> tetrahedron;

int vertcount = 1;
vector<vector<vector<int> > > rev_idx;

// vector<vector<vector<double> > > MAP;
// Density map stored in 3D array
int HEIGHT, WIDTH, DEPTH, LENGTH;
	
class EdgeHash{
	// Assume vertex are given in ascending order
	private:
		unordered_map<int, unordered_set<int>* >* v1;
		
	public:
		EdgeHash(){
			v1 = new unordered_map<int, unordered_set<int>* >();
			v1->clear();
		}
		bool HasEdge(cp edge){
			if(v1->count(edge.p1)==0){
				return false;
			}else{
				unordered_set<int>* v2 = v1->at(edge.p1);
				if (v2->count(edge.p2) == 0){
					return false;
				}else{
					return true;
				}
			}
		}
		void InsertEdge(cp edge){
			if (v1->count(edge.p1) > 0){
				unordered_set<int>* v2 = v1->at(edge.p1);
				v2->insert(edge.p2);
			}else{
				unordered_set<int>* v2 = new unordered_set<int>();
				v2->insert(edge.p2);
				std::pair<int, unordered_set<int>*> pair = std::make_pair(edge.p1, v2);
				v1->insert(pair);
			}
		}
};

class TriangleHash{
	// Assume vertex are given in ascending order
	private:
		unordered_map<int, EdgeHash* >* v1;
		
	public:
		TriangleHash(){
			v1 = new unordered_map<int, EdgeHash*>();
			v1->clear();
		}
		bool HasTriangle(tp triangle){
			if(v1->count(triangle.p1)==0){
				return false;
			}else{
				EdgeHash* v2 = v1->at(triangle.p1);
				cp edge;
				edge.p1 = triangle.p2;
				edge.p2 = triangle.p3;
				if (v2->HasEdge(edge)){
					return true;
				}else{
					return false;
				}
			}
		}
		void InsertTriangle(tp triangle){
			if (v1->count(triangle.p1) > 0){
				EdgeHash* v2 = v1->at(triangle.p1);
				cp edge;
				edge.p1 = triangle.p2;
				edge.p2 = triangle.p3;
				
				v2->InsertEdge(edge);
			}else{
				EdgeHash* v2 = new EdgeHash();
				cp edge;
				edge.p1 = triangle.p2;
				edge.p2 = triangle.p3;
				v2->InsertEdge(edge);
				std::pair<int, EdgeHash*> pair = std::make_pair(triangle.p1, v2);
				v1->insert(pair);
			}
		}
};


void bin_init(string filename){
	vertex.clear(); edge.clear(); triangle.clear(); 
	// tetrahedron.clear();
    ifstream binaryIO;
    binaryIO.open("mapinput.bin", ios::binary);
    if (!binaryIO.is_open()){
    	printf("Error opening file\n");
    	return;
    }
    
    char* header_data = new char[sizeof(int) * 4];
    
    binaryIO.read(header_data, sizeof(int) * 4);
	int* header_value = (int*) header_data;
	HEIGHT = header_value[0];
	WIDTH = header_value[1];
	DEPTH = header_value[2];
	LENGTH = header_value[3];

    printf("Preparing rev-index matrix %d-%d-%d\n", HEIGHT, WIDTH, DEPTH);
	
    rev_idx.clear();	
    for (int i = 0; i < HEIGHT; ++i) {
        // DIM 2
        vector<vector<int> > tmp2di;
        rev_idx.push_back(tmp2di);
        for (int j = 0; j < WIDTH; ++j){
            // DIM 3
            rev_idx[i].push_back(vector<int>(DEPTH, 0));
        }
    }

    printf("Reading density matrix: total %d Lines\n", LENGTH);
    char* density_data = new char[sizeof(double) * 4];
    double* density_value = (double*) density_data;
    for (int len = 0; len < LENGTH; ++len){
        if (DEBUG&&len%10000 == 0)
            printf("%d\n", len);
        int i, j, k;
        double v;
        
        binaryIO.read(density_data, sizeof(double) * 4);
		i = floor(density_value[0] + 0.5);
		j = floor(density_value[1] + 0.5);
		k = floor(density_value[2] + 0.5);
        v = density_value[3];

        i--;j--;k--;
        rev_idx[i][j][k] = vertcount++;
        point p;
        p.x = i; p.y = j; p.z = k; p.v = v;
        vertex.push_back(p);
    }
    binaryIO.close();
    printf("done\n");
}

int triangle_cube(int i, int j, int k, int AB, TriangleHash &th, EdgeHash &eh){
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
		int flag = 1;
        for (int cnt = 0; cnt<nb; cnt++){
            if ((i + TypeAtri[cnt*3][0]>=HEIGHT)||
                (i + TypeAtri[cnt*3+1][0]>=HEIGHT)||
                (i + TypeAtri[cnt*3+2][0]>=HEIGHT)||
                (j + TypeAtri[cnt*3][1]>=WIDTH)||
                (j + TypeAtri[cnt*3+1][1]>=WIDTH)||
                (j + TypeAtri[cnt*3+2][1]>=WIDTH)||
                (k + TypeAtri[cnt*3][2]>=DEPTH)||
                (k + TypeAtri[cnt*3+1][2]>=DEPTH)||
                (k + TypeAtri[cnt*3+2][2]>=DEPTH)){
                    flag = 0;
                    continue;
               }

            int sub1, sub2, sub3;
            if (rev_idx[i + TypeAtri[cnt*3][0]][j + TypeAtri[cnt*3][1]][k + TypeAtri[cnt*3][2]] <= 0){
                rev_idx[i + TypeAtri[cnt*3][0]][j + TypeAtri[cnt*3][1]][k + TypeAtri[cnt*3][2]] = vertcount++;
                point p;
                p.x = i+TypeAtri[cnt*3][0]; p.y = j+TypeAtri[cnt*3][1]; p.z = k+TypeAtri[cnt*3][2]; p.v = 1e-6;
                vertex.push_back(p);
                sub1 = rev_idx[i + TypeAtri[cnt*3][0]][j + TypeAtri[cnt*3][1]][k + TypeAtri[cnt*3][2]];
            }
            else{
                sub1 = rev_idx[i + TypeAtri[cnt*3][0]][j + TypeAtri[cnt*3][1]][k + TypeAtri[cnt*3][2]];
            }

            if (rev_idx[i + TypeAtri[cnt*3+1][0]][j + TypeAtri[cnt*3+1][1]][k + TypeAtri[cnt*3+1][2]] <= 0){
                rev_idx[i + TypeAtri[cnt*3+1][0]][j + TypeAtri[cnt*3+1][1]][k + TypeAtri[cnt*3+1][2]] = vertcount++;
                point p;
                p.x = i+TypeAtri[cnt*3+1][0]; p.y = j+TypeAtri[cnt*3+1][1]; p.z = k+TypeAtri[cnt*3+1][2]; p.v = 1e-6;
                vertex.push_back(p);
                sub2 = rev_idx[i + TypeAtri[cnt*3+1][0]][j + TypeAtri[cnt*3+1][1]][k + TypeAtri[cnt*3+1][2]];
            }
            else{
                sub2 = rev_idx[i + TypeAtri[cnt*3+1][0]][j + TypeAtri[cnt*3+1][1]][k + TypeAtri[cnt*3+1][2]];
            }

            if (rev_idx[i + TypeAtri[cnt*3+2][0]][j + TypeAtri[cnt*3+2][1]][k + TypeAtri[cnt*3+2][2]]<=0){
                rev_idx[i + TypeAtri[cnt*3+2][0]][j + TypeAtri[cnt*3+2][1]][k + TypeAtri[cnt*3+2][2]] = vertcount++;
                point p;
                p.x = i+TypeAtri[cnt*3+2][0]; p.y = j+TypeAtri[cnt*3+2][1]; p.z = k+TypeAtri[cnt*3+2][2]; p.v = 1e-6;
                vertex.push_back(p);
                sub3 = rev_idx[i + TypeAtri[cnt*3+2][0]][j + TypeAtri[cnt*3+2][1]][k + TypeAtri[cnt*3+2][2]];
            }
            else{
                sub3 = rev_idx[i + TypeAtri[cnt*3+2][0]][j + TypeAtri[cnt*3+2][1]][k + TypeAtri[cnt*3+2][2]];
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
        
        /*
        int Type4D[4 * 5][3] = {{0,0,0}, {0,0,1}, {1,0,1}, {0,1,1},
        						{0,0,0}, {1,0,0}, {1,0,1}, {1,1,0},
        						{0,1,1}, {0,0,0}, {1,1,0}, {0,1,0},
        						{0,1,1}, {1,0,1}, {1,1,0}, {1,1,1},
        						{0,1,1}, {1,0,1}, {1,1,0}, {0,0,0}
        						};
        for (int cnt = 0; cnt < 5; ++cnt){
        	if (!flag) break;
        	tetra tmp;
        	for (int vert_id = 0; vert_id < 4; ++vert_id){
        		tmp.p[vert_id] = rev_idx[i + Type4D[cnt * 4 + vert_id][0]]
        								[j + Type4D[cnt * 4 + vert_id][1]]
        								[k + Type4D[cnt * 4 + vert_id][2]];
        	}
        	tetrahedron.push_back(tmp);
        }*/
    }
    else{// 10 triangles of Type B
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
		int flag = 1;
        for (int cnt = 0; cnt<nb; cnt++){
            if ((i + TypeAtri[cnt*3][0]>=HEIGHT)||
                (i + TypeAtri[cnt*3+1][0]>=HEIGHT)||
                (i + TypeAtri[cnt*3+2][0]>=HEIGHT)||
                (j + TypeAtri[cnt*3][1]>=WIDTH)||
                (j + TypeAtri[cnt*3+1][1]>=WIDTH)||
                (j + TypeAtri[cnt*3+2][1]>=WIDTH)||
                (k + TypeAtri[cnt*3][2]>=DEPTH)||
                (k + TypeAtri[cnt*3+1][2]>=DEPTH)||
                (k + TypeAtri[cnt*3+2][2]>=DEPTH)){
                    flag = 0;
               	    continue;
               }
                
                
            int sub1, sub2, sub3;
            if (rev_idx[i + TypeAtri[cnt*3][0]][j + TypeAtri[cnt*3][1]][k + TypeAtri[cnt*3][2]] <= 0){
                rev_idx[i + TypeAtri[cnt*3][0]][j + TypeAtri[cnt*3][1]][k + TypeAtri[cnt*3][2]] = vertcount++;
                point p;
                p.x = i+TypeAtri[cnt*3][0]; p.y = j+TypeAtri[cnt*3][1]; p.z = k+TypeAtri[cnt*3][2]; p.v = 1e-6;
                vertex.push_back(p);
                sub1 = rev_idx[i + TypeAtri[cnt*3][0]][j + TypeAtri[cnt*3][1]][k + TypeAtri[cnt*3][2]];
            }
            else{
                sub1 = rev_idx[i + TypeAtri[cnt*3][0]][j + TypeAtri[cnt*3][1]][k + TypeAtri[cnt*3][2]];
            }

            if (rev_idx[i + TypeAtri[cnt*3+1][0]][j + TypeAtri[cnt*3+1][1]][k + TypeAtri[cnt*3+1][2]] <= 0){
                rev_idx[i + TypeAtri[cnt*3+1][0]][j + TypeAtri[cnt*3+1][1]][k + TypeAtri[cnt*3+1][2]] = vertcount++;
                point p;
                p.x = i+TypeAtri[cnt*3+1][0]; p.y = j+TypeAtri[cnt*3+1][1]; p.z = k+TypeAtri[cnt*3+1][2]; p.v = 1e-6;
                vertex.push_back(p);
                sub2 = rev_idx[i + TypeAtri[cnt*3+1][0]][j + TypeAtri[cnt*3+1][1]][k + TypeAtri[cnt*3+1][2]];
            }
            else{
                sub2 = rev_idx[i + TypeAtri[cnt*3+1][0]][j + TypeAtri[cnt*3+1][1]][k + TypeAtri[cnt*3+1][2]];
            }

            if (rev_idx[i + TypeAtri[cnt*3+2][0]][j + TypeAtri[cnt*3+2][1]][k + TypeAtri[cnt*3+2][2]]<=0){
                rev_idx[i + TypeAtri[cnt*3+2][0]][j + TypeAtri[cnt*3+2][1]][k + TypeAtri[cnt*3+2][2]] = vertcount++;
                point p;
                p.x = i+TypeAtri[cnt*3+2][0]; p.y = j+TypeAtri[cnt*3+2][1]; p.z = k+TypeAtri[cnt*3+2][2]; p.v = 1e-6;
                vertex.push_back(p);
                sub3 = rev_idx[i + TypeAtri[cnt*3+2][0]][j + TypeAtri[cnt*3+2][1]][k + TypeAtri[cnt*3+2][2]];
            }
            else{
                sub3 = rev_idx[i + TypeAtri[cnt*3+2][0]][j + TypeAtri[cnt*3+2][1]][k + TypeAtri[cnt*3+2][2]];
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
        
        /*
        int Type4D[4 * 5][3] = {{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1},
        						{0,1,0}, {1,0,0}, {1,1,0}, {1,1,1},
        						{0,1,0}, {1,1,1}, {0,0,1}, {0,1,1},
        						{0,0,1}, {1,0,1}, {1,1,1}, {1,0,0},
        						{0,0,1}, {1,1,1}, {0,1,0}, {1,0,0}
        						};
        for (int cnt = 0; cnt < 5; ++cnt){
        	if (!flag) break;
        	tetra tmp;
        	for (int vert_id = 0; vert_id < 4; ++vert_id){
        		tmp.p[vert_id] = rev_idx[i + Type4D[cnt * 4 + vert_id][0]]
        								[j + Type4D[cnt * 4 + vert_id][1]]
        								[k + Type4D[cnt * 4 + vert_id][2]];
        	}
        	tetrahedron.push_back(tmp);
        }*/
    }
    return 0;
}


void bin_output(){
    ofstream ofs("vert.bin",ios::binary);
    printf("writing vertex\n");
    char* vert = new char[sizeof(double) * 4];
    double* vert_buffer = (double*) vert;
    for (int i = 0; i < vertex.size(); i++){
    	vert_buffer[0] = vertex[i].x; vert_buffer[1] = vertex[i].y;
    	vert_buffer[2] = vertex[i].z; vert_buffer[3] = vertex[i].v;
		ofs.write(vert, sizeof(double) * 4);
    }
    ofs.close();

	ofs.open("edge.bin", ios::binary);
    printf("writing edge\n");
    char* edgechar = new char[sizeof(int) * 2];
    int* edge_buffer = (int*) edgechar;
    for (int i = 0; i < edge.size(); i++){
    	edge_buffer[0] = edge[i].p1; edge_buffer[1] = edge[i].p2;
        ofs.write(edgechar, sizeof(int) * 2);
    }
    ofs.close();

    ofs.open("triangle.bin", ios::binary);
    printf("writing triangle\n");
    char* trianglechar = new char[sizeof(int) * 3];
    int* triangle_buffer = (int*) trianglechar;
    for (int i = 0; i < triangle.size(); i++){
    	triangle_buffer[0] = triangle[i].p1;
    	triangle_buffer[1] = triangle[i].p2;
    	triangle_buffer[2] = triangle[i].p3;
        ofs.write(trianglechar, sizeof(int) * 3);
    }
    ofs.close();
    
    
    /*
    ofs.open("tetrahedron.bin", ios::binary);
    printf("writing tetrahedron\n");
    char* tetchar = new char[sizeof(int) * 4];
    int* tet_buffer = (int*) tetchar;
    for (int i = 0; i < tetrahedron.size(); i++){
		for (int j = 0; j < 4; ++j)
			tet_buffer[j] = tetrahedron[i].p[j];
        ofs.write(tetchar, sizeof(int) * 4);
    }
    ofs.close();*/
}


int triangulation_with_vertex(){
	TriangleHash th;
	EdgeHash eh;

    double THD = -1e-6;
	int original_total = vertex.size();
//    int counter = 0;
    for(int v = 0; v < original_total; ++v){
		int i, j, k, val;
		i = vertex[v].x; j = vertex[v].y; k = vertex[v].z; val = vertex[v].v;
		if (val < THD) continue;
		if ((i+j+k)%2==1){
			triangle_cube(i, j, k, 0, th, eh);
		}
		else{
			triangle_cube(i, j, k, 1, th, eh);
		}
    }
    printf("Done\n");
    return 0;
}


int main()
{
    string filename = "mapinput.txt";

    printf("Initializing input\n");
	bin_init(filename);

    printf("Computing triangulation\n");
	triangulation_with_vertex();

    printf("Writing output\n");
	bin_output();

    printf("Done\n");
    return 0;
}
