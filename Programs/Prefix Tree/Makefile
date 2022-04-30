C = g++
CFLAGS = -Wall -g -std=c++17

prefixTreeConstructs: prefixTree.o main.o
	$(C) $(CFLAGS) -o prefixTreeConstructs prefixTree.o main.o
prefixTree.o: prefixTree.cpp prefixTree.h
	$(C) $(CFLAGS) -c prefixTree.cpp
main.o: main.cpp
	$(C) $(CFLAGS) -c main.cpp
clean:	
	rm prefixTreeConstructs prefixTree.o main.o