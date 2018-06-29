# the compiler
CXX = g++

# flags:
# -static-libstdc++ might be needed if running in matlab
CXXFLAGS = -std=c++11 -w

# includes
COREINCLUDES = -I./extern/phat/include
TRI_INCLUDES = -I./extern/boost

# target
CORE = core/DiMorSC.cpp core/DiscreteVField.h core/persistence.h core/Simplex.h core/Simplicial2Complex.h
TRI = Triangulate
TREE = graph2tree

all: DiMorSC $(TRI) $(TREE)
	
DiMorSC: $(CORE)
	mkdir -p bin
	mkdir -p output
	$(CXX) $(CXXFLAGS) $(COREINCLUDES) -o bin/DiMorSC core/DiMorSC.cpp

$(TRI): pointcloud/$(TRI).cpp
	$(CXX) $(CXXFLAGS) $(TRI_INCLUDES) -o bin/$(TRI) pointcloud/$(TRI).cpp

$(TREE): tree_simplification/$(TREE).cpp
	$(CXX) $(CXXFLAGS) -o bin/$(TREE) tree_simplification/$(TREE).cpp tree_simplification/graph.cpp tree_simplification/readini.cpp
#clean:
	
	
# Test Command
# ./bin/DiMorSC data/OP_7_trunc.bin output/OP7 5 3
