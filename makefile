# the compiler
CC = g++

# flags:
CFLAGS = -static-libstdc++ -std=c++11 -w

# includes
INCLUDES = -I./core/phat/include
# target
TARGET = DiMorSC

all: $(TARGET)

$(TARGET): core/$(TARGET).cpp
	mkdir -p bin
	mkdir -p output
	$(CC) $(CFLAGS) $(INCLUDES) -o bin/$(TARGET) core/$(TARGET).cpp

#clean:
	
	
# Test Command
# ./bin/DiMorSC data/OP_7_trunc.bin output/OP7 5 3
