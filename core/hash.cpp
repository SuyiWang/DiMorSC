#include"hash.h"


//**********vertex hasing implementation***************
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

int VertexHash::GetIndex(point p){
	if (v_hash.count(p) > 0)
		return v_hash[p];
	else
		return -1;
}

void VertexHash::InsertVertex(point p, int n){
	v_hash.insert(std::make_pair(p, n));
}
//**********vertex hasing implementation***************


//**********edge hashing implementation****************
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

void EdgeHash::InsertEdge(cp edge){
	e_hash.insert(edge);
}

bool EdgeHash::HasEdge(cp edge){
	return e_hash.count(edge) > 0;
}
//**********edge hashing implementation****************


//**********triangle hashing implementation************
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

bool TriangleHash::HasTriangle(tp triangle){
	return t_hash.count(triangle) > 0;
}

void TriangleHash::InsertTriangle(tp triangle){
	t_hash.insert(triangle);
}
//**********triangle hashing implementation************
