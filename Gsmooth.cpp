/*
Smooth point cloud 2D/3D.
Point cloud must have integer coordinates.
Author: Suyi Wang

Input: output/<id>_dens.bin
Output: <id>_dens.bin
Complie: g++ Gsmooth.cpp -static-libstdc++ -I boost_1_64_0/ -std=c++11 -o Gsmooth -w
*/


#include<stdio.h>
#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<algorithm>
#include<time.h>
#include<cstdlib>
// #include<unordered_map>
// #include<unordered_set>

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>

using namespace std;

#define DEBUG 0

struct point{
	int x,y,z;
	double v;
	bool in_bbox;
};

int vertcount = 0;
vector<point> vertex;


// *********** begin vertex pair***********
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

VertexHash vh;

void init_3D(double selTHD, string filename){
	vertex.clear();
    ifstream binaryIO;
    
    if (DEBUG) cout << "Reading from " << filename << endl;
    binaryIO.open(filename.c_str(), ios::binary);
    if (!binaryIO.is_open()){
    	printf("Error opening file\n");
    	return;
    }
    
    char* header_data = new char[sizeof(int) * 1];
    
    binaryIO.read(header_data, sizeof(int) * 1);
	int* header_value = (int*) header_data;
	int LENGTH = header_value[0];
    if (DEBUG) cout << "Reading " << LENGTH << " points." << endl;

	cout << *header_value << endl;
    printf("Reading 3D density matrix: total %d Lines\n", LENGTH);
    char* density_data = new char[sizeof(double) * 4];
    double* density_value = (double*) density_data;
    for (int len = 0; len < LENGTH; ++len){
        if (DEBUG&&len%10000 == 0)
            printf("%d\n", len);
        int i, j, k;
        double v;
        
        binaryIO.read(density_data, sizeof(double) * 4);
		i = round(density_value[0]);
		j = round(density_value[1]);
        k = round(density_value[2]);
        v = density_value[3];
        if (v < selTHD) continue;

        point p;
        p.x = i; p.y = j; p.z = k; p.v = v;
        int idx = -1;
        idx = vh.GetIndex(p);
        if (idx < 0){
	        vertex.push_back(p);
	        vh.InsertVertex(p, vertcount);
	        vertcount++;
	    }
    }
    binaryIO.close();
    printf("done\n");
}


void bin_output(string id, double selTHD){
    string vname = id + "_dens.bin";
    ofstream ofs(vname,ios::binary);
    printf("writing vertex\n");
    
    
    char* header_data = new char[sizeof(int) * 1];
    int* header_buffer = (int*) header_data;
    header_buffer[0] = vertex.size();
    ofs.write(header_data, sizeof(int) * 1);
    if (DEBUG) cout << "Writing " << header_buffer[0] << " points" << endl;
    
    
    char* vert = new char[sizeof(double) * 4];
    double* vert_buffer = (double*) vert;
    for (int i = 0; i < vertex.size(); i++){
        if (vertex[i].v < selTHD) continue;
    	vert_buffer[0] = vertex[i].x; vert_buffer[1] = vertex[i].y;
    	vert_buffer[2] = vertex[i].z; vert_buffer[3] = vertex[i].v;
		ofs.write(vert, sizeof(double) * 4);
    }
    ofs.close();
}

//  TO DO
void init_2D(string filename){
	vertex.clear();
    ifstream binaryIO;
    binaryIO.open(filename.c_str(), ios::binary);
    if (!binaryIO.is_open()){
    	printf("Error opening file\n");
    	return;
    }
    
    char* header_data = new char[sizeof(int) * 1];
    
    binaryIO.read(header_data, sizeof(int) * 1);
	int* header_value = (int*) header_data;
	int LENGTH = header_value[0];

	cout << *header_value << endl;
    printf("Reading 2D density matrix: total %d Lines\n", LENGTH);
    char* density_data = new char[sizeof(double) * 3];
    double* density_value = (double*) density_data;
    for (int len = 0; len < LENGTH; ++len){
        if (DEBUG&&len%10000 == 0)
            printf("%d\n", len);
        int i, j, k;
        double v;
        
        binaryIO.read(density_data, sizeof(double) * 3);
		i = floor(density_value[0] + 0.5);
		j = floor(density_value[1] + 0.5);
        v = density_value[2];

        point p;
        p.x = i; p.y = j; p.z = 0; p.v = v;
        int idx = -1;
        idx = vh.GetIndex(p);
        if (idx < 0){
	        vertex.push_back(p);
	        vh.InsertVertex(p, vertcount);
	        vertcount++;
	    }
    }
    binaryIO.close();
    printf("done\n");
}


vector<double> read_ini(string ini_filename){
    vector<double> rtn;
    rtn.clear();
    
    FILE* fp = fopen(ini_filename.c_str(), "r");
    if (fp == NULL){
        cerr << "Failed to open " << ini_filename << endl;
        return rtn;
    }
    double rd = 0;
    for(int i = 0; i < 6; ++i){
        if(fscanf(fp, "%lf", &rd) != EOF){
            rtn.push_back(rd);
            if (DEBUG) cout << rd << endl;
        }
        else
            cerr << "Reached EOF, ini file error\n";
    }
    fclose(fp);
    return rtn;
}

double norm_const = 1;
int diffuse(point p, vector<double> stepsize){
	double sigma = 3;
	int sum = 0;
	for (int i = -2; i <= 2; ++i)
		for(int j = -2; j <= 2; ++j)
			for(int k = -1; k <= 1; ++k){
				point dif_p;
				
				dif_p.x = p.x + i;
				dif_p.y = p.y + j;
				dif_p.z = p.z + k;
                double x = i*stepsize[0];
                double y = j*stepsize[1];
                double z = k*stepsize[2];
				dif_p.v = p.v * norm_const * exp(-(x*x + y*y + z*z)/(2* sigma * sigma));
				
				/*if (dif_p.v < 1e-8) {
					continue;
				}*/
				
				int idx = vh.GetIndex(dif_p);
				if (idx < 0){
					vh.InsertVertex(dif_p, vertcount);
					vertex.push_back(dif_p);
					vertcount++;
				}else{
					vertex[idx].v += dif_p.v;
				}
    }
    return sum;
}

double kernel_init(double sigma, vector<double> stepsize){
	double sum = 0;
    for(int k = -1; k <= 1; ++k){
        for (int i = -2; i <= 2; ++i)
            for(int j = -2; j <= 2; ++j){
                double x = i*stepsize[0];
                double y = j*stepsize[1];
                double z = k*stepsize[2];
				sum += exp(-(x*x + y*y + z*z)/(2* sigma * sigma));
				if (DEBUG) cout << exp(-(x*x + y*y + z*z)/(2* sigma * sigma)) 
                                << " ";
			}
		if (DEBUG) cout << "\n\n";
    }
    return 1.0/sum;
}


void smooth3D(vector<double> steps){
    int LEN = vertex.size();
    for(int i = 0; i < LEN; ++i){
        if (DEBUG && i%1000 == 0) cout << i << endl;
        diffuse(vertex[i], steps);
    }
}


int main(int argc, char* argv[])
{
	if (argc != 2){
		cout << "usage: Gsmooth <ini_file_name>\n";
		return 0;
	}
    string inifilename(argv[1]);
    vector<double> ini = read_ini(inifilename);
    
    // Cutting threshold - usually very low.
    double selTHD = round(ini[0]);
    // Folder ID
    char idc[20];
    int idi = round(ini[1]);

    sprintf(idc, "%d", idi);
    string id(idc);
    
    // Dimension 2 or 3
	int dimension = round(ini[2]);
    
    // Init step size;
    vector<double> steps;
    steps.clear();
	
    // Input filename
    string filename = "output/" + id + "_dens.bin";

    printf("Initializing input\n");
    if (dimension ==3){
        steps.push_back(ini[3]);
        steps.push_back(ini[4]);
        steps.push_back(ini[5]);
		init_3D(selTHD, filename);
    }else{
        steps.push_back(ini[3]);
        steps.push_back(ini[4]);
		init_2D(filename);
    }
    norm_const = kernel_init(3.0, steps);
    printf("Applying Gaussian Kernel\n");
    if (dimension == 3)
		smooth3D(steps);
	else
		smooth3D(steps);// should be 2D

    printf("Writing output\n");
    bin_output(id, selTHD);
    printf("Done\n");
    return 0;
}
