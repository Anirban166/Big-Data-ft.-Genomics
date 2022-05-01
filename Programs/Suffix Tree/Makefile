C = g++
CFLAGS = -Wall -g -std=c++17

suffixTreeConstructs: suffixTree.o main.o
	$(C) $(CFLAGS) -o suffixTreeConstructs suffixTree.o main.o
suffixTree.o: suffixTree.cpp suffixTree.h
	$(C) $(CFLAGS) -c suffixTree.cpp
main.o: main.cpp
	$(C) $(CFLAGS) -c main.cpp
clean:	
	rm suffixTreeConstructs suffixTree.o main.o