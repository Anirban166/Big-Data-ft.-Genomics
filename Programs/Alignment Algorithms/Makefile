C = g++
CFLAGS = -Wall -g -std=c++17

alignmentConstructs: BLAST.o main.o
	$(C) $(CFLAGS) -o alignmentConstructs BLAST.o main.o
BLAST.o: BLAST.cpp BLAST.h
	$(C) $(CFLAGS) -c BLAST.cpp
main.o: main.cpp
	$(C) $(CFLAGS) -c main.cpp
clean:	
	rm alignmentConstructs BLAST.o main.o