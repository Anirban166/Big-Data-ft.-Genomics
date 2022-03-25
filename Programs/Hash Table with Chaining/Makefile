CXX = g++
CXXFLAGS = -Wall -g -std=c++17

hashTableConstructs: FASTAreadset_HT.o main.o
	$(CXX) $(CXXFLAGS) -o hashTableConstructs FASTAreadset_HT.o main.o
FASTAreadset_HT.o: FASTAreadset_HT.cpp FASTAreadset_HT.h
	$(CXX) $(CXXFLAGS) -c FASTAreadset_HT.cpp
main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c main.cpp
clean:	
	rm hashTableConstructs FASTAreadset_HT.o main.o