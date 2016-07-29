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

vector<point> vertex;
vector<cp> edge;
vector<tp> triangle;

int vertcount = 1;
vector<vector<vector<int> > > rev_idx;

vector<vector<vector<double> > > MAP;
	// Density map stored in 3D array
	
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


void init(vector<vector<vector<double> > > &MAP, string filename){
    FILE* fp = fopen(filename.c_str(), "r");
    int HEIGHT, WIDTH, DEPTH, LENGTH;

    printf("Preparing density matrix\n");
    fscanf(fp, "%d %d %d %d", &HEIGHT, &WIDTH, &DEPTH, &LENGTH);

	// HEIGHT = HEIGHT + 1; WIDTH = WIDTH + 1; LENGTH = LENGTH + 1;
	
    MAP.clear();
    rev_idx.clear();

    for (int i = 0; i < HEIGHT; ++i) {
        // DIM 2
        vector<vector<double> > tmp2dd;
        vector<vector<int> > tmp2di;
        MAP.push_back(tmp2dd);
        rev_idx.push_back(tmp2di);
        for (int j = 0; j < WIDTH; ++j){
            // DIM 3
            MAP[i].push_back(vector<double>(DEPTH, 0));
            rev_idx[i].push_back(vector<int>(DEPTH, 0));
        }
    }

    printf("Reading density matrix: total %d Lines\n", LENGTH);
    for (int len = 0; len < LENGTH; len++){
        if (DEBUG&&len%10000 == 0)
            printf("%d\n", len);
        int i, j, k;
        double v;
        fscanf(fp, "%d %d %d %lf", &i, &j, &k, &v);
        i--;j--;k--;
        MAP[i][j][k] = v;
        rev_idx[i][j][k] = vertcount++;
        point p;
        p.x = i; p.y = j; p.z = k; p.v = v;
        vertex.push_back(p);
    }
    fclose(fp);
    printf("done\n");
}

void bin_init(vector<vector<vector<double> > > &MAP, string filename){
    ifstream binaryIO;
    binaryIO.open("mapinput.bin", ios::binary);
    if (!binaryIO.is_open()){
    	printf("Error opening file\n");
    	return;
    }
    
    char* header_data = new char[sizeof(int) * 4];
    int HEIGHT, WIDTH, DEPTH, LENGTH;
    binaryIO.read(header_data, sizeof(int) * 4);
	int* header_value = (int*) header_data;
	HEIGHT = header_value[0];
	WIDTH = header_value[1];
	DEPTH = header_value[2];
	LENGTH = header_value[3];

    printf("Preparing density matrix %d-%d-%d\n", HEIGHT, WIDTH, DEPTH);
	
    MAP.clear();
    rev_idx.clear();	
    for (int i = 0; i < HEIGHT; ++i) {
        // DIM 2
        vector<vector<double> > tmp2dd;
        vector<vector<int> > tmp2di;
        MAP.push_back(tmp2dd);
        rev_idx.push_back(tmp2di);
        for (int j = 0; j < WIDTH; ++j){
            // DIM 3
            MAP[i].push_back(vector<double>(DEPTH, 0));
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
        MAP[i][j][k] = v;
        rev_idx[i][j][k] = vertcount++;
        point p;
        p.x = i; p.y = j; p.z = k; p.v = v;
        vertex.push_back(p);
    }
    binaryIO.close();
    printf("done\n");
}

int triangle_cube(int i, int j, int k, int AB, TriangleHash &th, EdgeHash &eh){
    int HEIGHT, WIDTH, DEPTH;
	
    HEIGHT = MAP.size();
    WIDTH  = MAP[0].size();
    DEPTH  = MAP[0][0].size();
	
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
                   continue;
               }
//
//            dens1 = MAP[i + TypeAtri[cnt*3][0]][j + TypeAtri[cnt*3][1]][k + TypeAtri[cnt*3][2]];
//            dens2 = MAP[i + TypeAtri[cnt*3+1][0]][j + TypeAtri[cnt*3+1][1]][k + TypeAtri[cnt*3+1][2]];
//            dens3 = MAP[i + TypeAtri[cnt*3+2][0]][j + TypeAtri[cnt*3+2][1]][k + TypeAtri[cnt*3+2][2]];
//            if (dens1>THD && dens2 > THD && dens3 >THD){
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
//                nonempty(cnt) = 1;
//            }
        }
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
                   continue;
               }
                int sub1, sub2, sub3;
//
//            dens1 = MAP[i + TypeAtri[cnt*3][0]][j + TypeAtri[cnt*3][1]][k + TypeAtri[cnt*3][2]];
//            dens2 = MAP[i + TypeAtri[cnt*3+1][0]][j + TypeAtri[cnt*3+1][1]][k + TypeAtri[cnt*3+1][2]];
//            dens3 = MAP[i + TypeAtri[cnt*3+2][0]][j + TypeAtri[cnt*3+2][1]][k + TypeAtri[cnt*3+2][2]];
//            if (dens1>THD && dens2 > THD && dens3 >THD){
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
//                nonempty(cnt) = 1;
//            }
        }
    }
    return 0;
}

void output(){
    FILE* fp = fopen("vert.txt","w");
    printf("writing vertex\n");
    for (int i = 0; i < vertex.size(); i++)
        fprintf(fp, "%f %f %f %f\n", vertex[i].x, vertex[i].y, vertex[i].z, vertex[i].v);
    fclose(fp);

    printf("writing edge\n");
    fp = fopen("edge.txt","w");
    for (int i = 0; i < edge.size(); i++)
        fprintf(fp, "%d %d\n", edge[i].p1, edge[i].p2);
    fclose(fp);

    printf("writing triangle\n");
    fp = fopen("triangle.txt","w");
    for (int i = 0; i < triangle.size(); i++)
        fprintf(fp, "%d %d %d\n", triangle[i].p1, triangle[i].p2, triangle[i].p3);
    fclose(fp);
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
}


int triangulation(vector<vector<vector<double> > > &MAP){
    int HEIGHT, WIDTH, DEPTH;
    HEIGHT = MAP.size();
    WIDTH  = MAP[0].size();
    DEPTH  = MAP[0][0].size();
	
	TriangleHash th;
	EdgeHash eh;

    double THD = 1e-6;
//    int counter = 0;
    for(int i = 0; i<HEIGHT; i++)
        for(int j = 0; j<WIDTH; j++)
            for(int k = 0; k<DEPTH; k++){
                if (MAP[i][j][k] < THD) continue;
                //printf("%d %f\n", rev_idx[i][j][k], MAP[i][j][k]);
                if ((i+j+k)%2==1){
                    triangle_cube(i, j, k, 0, th, eh);
                }
                else{
                    triangle_cube(i, j, k, 1, th, eh);
                }
//                counter ++;
//                if (counter % 1000 ==0) printf("%d\n", counter);
            }
    printf("Done\n");
    return 0;
}

int triangulation_with_vertex(vector<vector<vector<double> > > &MAP){
	TriangleHash th;
	EdgeHash eh;

    double THD = 1e-6;
	int original_total = vertex.size();
//    int counter = 0;
    for(int v = 0; v < original_total; ++v){
		int i, j, k;
		i = vertex[v].x; j = vertex[v].y; k = vertex[v].z;
		if (MAP[i][j][k] < THD) continue;
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
//    init(MAP, filename);
	bin_init(MAP, filename);

    printf("Computing triangulation\n");
//  triangulation(MAP);
	triangulation_with_vertex(MAP);

    printf("Writing output\n");
//    output();
	bin_output();

    printf("Done\n");
    return 0;
}
