CXX = g++
CXXFLAGS = -Wall -g -std=c++17

linkedListConstructs: FASTAreadset_LL.o main.o
	$(CXX) $(CXXFLAGS) -o linkedListConstructs FASTAreadset_LL.o main.o
FASTAreadset_LL.o: FASTAreadset_LL.cpp FASTAreadset_LL.h
	$(CXX) $(CXXFLAGS) -c FASTAreadset_LL.cpp
main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c main.cpp
clean:	
	rm linkedListConstructs FASTAreadset_LL.o main.o
