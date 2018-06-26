#include<iostream>
#include<fstream>
#include"graph.h"
#include"readini.h"

using namespace std;
#define DEBUG 1

/*
	Enforce tree structure using shortest path tree
	Input: Graph file (output from DiMorSC)
	Output: Tree represented by adjacency list

	
	saddle_threshold: Removes saddles whose density is lower than threshold
	pos: Root location
*/

/*
	usage: ./graph2tree <graph_file.ini>
	graph_file.ini: if does not exist, create one
	format:
		<input vertex filename>
		<input edge filename>
		<output folder+prefix>
		<root>
		<threshold>
		[Component #]
		[shift]
		
	default:
		# vertex filename
		0_vert.txt
		# edge filename
		0_edge.txt
		# output prefix: program writes to output/0.swc
		output/0
		# root
		0 0 0
		# threshold for simplification
		5
		[not specified: use maximum component]
		[not specified: use default]
	
	outputs to 0.swc or 0_0.swc ~ 0_x.swc if n connected components are available
*/

int main(int argc, char* argv[]){
	parameters para;
	string ininame;
	
	//  Resolving parameters
    if (argc != 2){
		cout << "Usage: ./graph2tree graphfile.ini"	<<endl;
		cout << "default output/0.ini is generated" << endl;
		default_ini();
		return 0;
    }else{
		ininame = string(argv[1]);
		
		int stat = loadini(ininame, para);
		if (stat == -1){
			cout << "file not correctly read, exit\n";
			return 0;
		}
	}
    
	//	Processing Graph
	//  graph G(para.ininame);
	cout << "Loading Graph from " + para.vertfile + " " + para.edgefile << endl;
	graph G(para.vertfile, para.edgefile);
	cout << "checking vertex and edge redundancy\n";
	G.check_redundancy();
	cout << "Counting Components: ";
	vector<graph> subgraph = G.split(para.comp);
	//cout << subgraph.size() << endl;

	int counter = 0;
	for (auto subg : subgraph){
		if (subg.size() < para.comp) continue;
		vector<vector<int> > edge = subg.dijkstra(para.root);
		subg.set_edge(edge);
		//subg.to_swc(para.outputprefix + to_string(counter) + ".swc");
		subg.to_file(para.outputprefix + "_tree_" + to_string(counter++));
	}
	cout << "Components written: " << counter << endl;
}
