#include"graph.h"

graph::graph(){
	v.clear();
	e.clear();
	dist_created=0;
}

graph::graph(string vertfile, string edgefile){
	graph();
	loadvert(vertfile);
	loadedge(edgefile);
}


int graph::loadvert(const string & filename){
	FILE* vertinput = fopen(filename.c_str(), "r");
	int x,y,z,c;
	double f;
	while(fscanf(vertinput, "%d%d%d%lf%d",&x,&y,&z,&f,&c)!=EOF){
		//printf("%d %d %d %f %d\n", x,y,z,f,c);
		point p;
		p.pos.push_back(x);p.pos.push_back(y);p.pos.push_back(z);
		p.f = f;
		v.push_back(p);
		e.push_back(vector<int>());
	}
	fclose(vertinput);
	return 0;
}


int graph::loadedge(const string & filename){
	FILE* edgeinput = fopen(filename.c_str(), "r");
	int x,y,c;
	double persist;
	while(fscanf(edgeinput, "%d%d%d%lf",&x,&y,&c,&persist)!=EOF){
		//printf("%d %d %d %f\n", x,y,c,persist);
		// maybe need check redundancy.
		// the edge is bi-directional
		e[x-1].push_back(y-1);
		//cout << x-1 << " " << y-1 <<endl;
		e[y-1].push_back(x-1);
	}
	fclose(edgeinput);
	return 0;
}


int graph::check_redundancy(){
	// check zero connectivity
	vector<bool> rmvmark(v.size(), 0);
	bool zerovert = false;
	for(auto i = 0; i < e.size(); ++i)
		if (e[i].size() == 0) {
			rmvmark[i] = 1;
			zerovert = true;
		}

	if (zerovert){
		cout << "Found and shrunk isolated vert.\n";
		int realpointer = 0;
		for(auto i = 0; i < e.size(); ++i){
			// skip without moving data
			if (rmvmark[i]) continue;
				
			// copy data if i is in consistent with real
			if(i != realpointer){
				v[realpointer] = v[i];
				e[realpointer] = e[i];
			}
			realpointer++;
		}
	}

	// check redundant edge
	bool erased = false;
	for(auto i = 0; i < e.size(); ++i){
		// remove redundancy using unique
		// if elemenent is huge, use set
		sort( e[i].begin(), e[i].end() );
		auto realend = unique(e[i].begin(), e[i].end());
		if (realend!=e[i].end()){
			e[i].erase(realend, e[i].end());
			erased = true;
		}
	}
	if (erased) cout << "Duplicate edges removed\n";
	return 0;
}


graph graph::traverse(const int start, const vector<point> & v, const vector<vector<int> > & e, 
					  vector<bool> & touched){
	graph rtngraph;
	// uses original idx
	vector<int> queued_idx;
	int visit=0;
	queued_idx.push_back(start);
	touched[start] = true;
	// mapping original idx to new idx
	unordered_map<int, int> idxmap;
	idxmap[start] = 0;
	rtngraph.v.push_back(v[start]);
	rtngraph.e.push_back(vector<int>());


	// loop when there is unvisited vertex
	while(visit < queued_idx.size()){
		int now = queued_idx[visit++];

		// for every touched vert, put neighbours into queue
		for(auto vert:e[now]){
			if (!touched[vert]){
				touched[vert] = true;
				queued_idx.push_back(vert);
				// copy new vert
				idxmap[vert] = rtngraph.v.size();
				rtngraph.v.push_back(v[vert]);
				rtngraph.e.push_back(vector<int>());
			}
			// copy edges
			// all appeared vert should be touched
			rtngraph.e[idxmap[now]].push_back(idxmap[vert]);
			rtngraph.e[idxmap[vert]].push_back(idxmap[now]);
		}
	}

	return rtngraph;
}


vector<graph> graph::split(int n){
	vector<graph> rtn;

	vector<bool> touched(v.size(), false);
	int vnum;
	for(vnum = 0; vnum < (int)v.size(); vnum++){
		if (!touched[vnum]){
			graph g = traverse(vnum, v, e, touched);
			rtn.push_back(g);
		}
	}
	cout << "# of components: " << to_string(rtn.size()) << endl;
	return rtn;
}

int graph::size(){
	return v.size();
}


double graph::get_dist(point a, point b){
	double sum = 0;
	// assert a.pos.size() == b.pos.size()
	for(auto i = 0; i < a.pos.size(); i++)
		sum += (a.pos[i]-b.pos[i]) * (a.pos[i]-b.pos[i]);
	return sqrt(sum);
}


int graph::find_vert(vector<int> pos){
	if (pos.size() != 3){
		cout << "position size mismatch\n";
		return -1;
	}
	int rtn = 0;
	double d = get_dist(v[0], point(pos));

	for(auto idx = 1; idx < v.size(); idx++){
		double newd = get_dist(v[rtn], v[idx]);
		if (newd < d){
			d = newd; rtn = idx;
		}
	}
	return rtn;
}


vector<vector<int> > graph::dijkstra(const vector<int> &pos){
	int root = find_vert(pos);
	cout << "Root idx: " << root << endl;
	
	vector<double> dist(v.size(), std::numeric_limits<double>::infinity());
	vector<vector<int> > tree_edge(v.size(), vector<int>());

	dist[root] = 0;
	// using priority queue is not optimal algorithm, but 
	// it is still near-linear and is much easier to develop.
	// Optimal: Fibonacci heap.
	priority_queue<vertex, std::vector<vertex>, vertex_cmp> pq;
	pq.push(vertex(root, -1, 0));

	while(!pq.empty()){
		vertex dij_v = pq.top();
		pq.pop();
		double cur_d = dist[dij_v.idx];

		if (dij_v.idx != root){
			tree_edge[dij_v.prev].push_back(dij_v.idx);
			tree_edge[dij_v.idx].push_back(dij_v.prev);
		}
		//cout << "At " << dij_v.idx << " found " << endl;

		for(auto adj_v : e[dij_v.idx]){
			double d = get_dist(v[dij_v.idx], v[adj_v]);
			//cout << adj_v << " : " << d << "<>" << dist[adj_v] << endl;
			if (cur_d + d < dist[adj_v]){
				dist[adj_v] = cur_d + d;
				pq.push(vertex(adj_v, dij_v.idx, cur_d + d));
			}
		}

		while(!pq.empty() && dist[pq.top().idx] < pq.top().f) pq.pop();
	}

	for(auto i = 0; i < v.size(); ++i)
		v[i].f = dist[i];

	dist_created = 1;
	return tree_edge;
}


int graph::set_edge(const vector<vector<int> > &tree_edge){
	e = tree_edge;
	return 0;
}

int graph::to_swc(string filename){
	return 0;
}

int graph::to_file(string filename){
	ofstream vertout;
	vertout.open(filename+"_vert.txt");
	for(auto vert : v){
		vertout << vert << endl;
	}
	vertout.close();

	ofstream edgeout;
	edgeout.open(filename+"_edge.txt");
	for(int v0 = 0; v0 < (int)e.size(); v0++){
		for(auto v : e[v0]){
			if (v0 < v)
				// write only one pair of edge
				edgeout << make_pair(v0+1, v+1) << endl;
		}
	}
	edgeout.close();
	return 0;
}

