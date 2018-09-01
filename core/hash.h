#include "geometry.h"

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>

// *********** begin vertex pair***************
class VertexHash{
	private:
		boost::unordered_map<point, int> v_hash;
	public:
		VertexHash(){
			v_hash.clear();
		}
		int GetIndex(point p);
		void InsertVertex(point p, int n);
		int size(){return v_hash.size();}
		boost::unordered_map<point, int>::iterator begin(){
			return v_hash.begin();
		}
		boost::unordered_map<point, int>::iterator end(){
			return v_hash.end();
		}
};
// *********** End of vertex pair ***************


// *********** begin Edge pair ******************
class EdgeHash{
	// Assuming vertex are given in ascending order
	private:
		boost::unordered_set<cp> e_hash;
	
	public:
		EdgeHash(){
			e_hash.clear();
		}
		bool HasEdge(cp edge);
		void InsertEdge(cp edge);
		int size(){return e_hash.size();}
		boost::unordered_set<cp>::iterator begin(){return e_hash.begin();}
		boost::unordered_set<cp>::iterator end(){return e_hash.end();}
};
// *********** End of Edge pair******************

// *********** begin Triangle pair***************
class TriangleHash{
	// Assuming vertex are given in ascending order
	private:
		boost::unordered_set<tp> t_hash;
	public:
		TriangleHash(){
			t_hash.clear();
		}
		bool HasTriangle(tp triangle);
		void InsertTriangle(tp triangle);
		int size(){return t_hash.size();}
		boost::unordered_set<tp>::iterator begin(){return t_hash.begin();}
		boost::unordered_set<tp>::iterator end(){return t_hash.end();}
};
// *********** End of Triangle pair *************
