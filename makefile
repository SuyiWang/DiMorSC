# the compiler
CC = g++

# flags:
# -static-libstdc++ might be needed if running in matlab
CFLAGS = -std=c++11 -w

# includes
INCLUDES = -I./core/phat/include
TRI_INCLUDES = -I./core/boost

# target
TARGET = DiMorSC
TRI = Triangulation
TREE = graph2tree

all: $(TARGET) $(TRI) $(TREE)

$(TARGET): core/$(TARGET).cpp
	mkdir -p bin
	mkdir -p output
	$(CC) $(CFLAGS) $(INCLUDES) -o bin/$(TARGET) core/$(TARGET).cpp

$(TRI): core/$(TRI).cpp
	$(CC) $(CFLAGS) $(TRI_INCLUDES) -o bin/Triangulate core/$(TRI).cpp

$(TREE): image_processing/$(TREE).cpp
	$(CC) $(CFLAGS) -o bin/$(TREE) image_processing/$(TREE).cpp image_processing/graph.cpp image_processing/readini.cpp
#clean:
	
	
# Test Command
# ./bin/DiMorSC data/OP_7_trunc.bin output/OP7 5 3
