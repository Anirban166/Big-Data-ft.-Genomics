C = g++
CFLAGS = -Wall -g -std=c++17

linkedListConstructs: FASTAreadset_LL.o main.o
	$(C) $(CFLAGS) -o linkedListConstructs FASTAreadset_LL.o main.o
FASTAreadset_LL.o: FASTAreadset_LL.cpp FASTAreadset_LL.h
	$(C) $(CFLAGS) -c FASTAreadset_LL.cpp
main.o: main.cpp
	$(C) $(CFLAGS) -c main.cpp
clean:	
	rm linkedListConstructs FASTAreadset_LL.o main.o
