C = g++
CFLAGS = -Wall -g -std=c++17

hashTableConstructs: FASTAreadset_HT.o main.o
	$(C) $(CFLAGS) -o hashTableConstructs FASTAreadset_HT.o main.o
FASTAreadset_HT.o: FASTAreadset_HT.cpp FASTAreadset_HT.h
	$(C) $(CFLAGS) -c FASTAreadset_HT.cpp
main.o: main.cpp
	$(C) $(CFLAGS) -c main.cpp
clean:	
	rm hashTableConstructs FASTAreadset_HT.o main.o