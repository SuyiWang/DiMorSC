#include<algorithm>
// for vertex
class point{
	public:
    int x,y,z;
    double v;
    bool in_bbox;
};

// for edge
class cp{
	public:
		int p1,p2;  // must be vertices after merge.
		void Reorder(){
			if (p1>p2) std::swap(p1, p2);
		}
};

// for triangle
class tp{
	public:
		int p1,p2,p3;
		void Reorder(){
			if (p1>p2) std::swap(p1, p2);
			if (p1>p3) std::swap(p1, p3);
			if (p2>p3) std::swap(p2, p3);
		}
};
