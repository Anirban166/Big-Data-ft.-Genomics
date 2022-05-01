/*-------------------------
  Author: Anirban166/Ani
  Email:  ac4743@nau.edu
--------------------------*/
#pragma once
#include <cmath>
#include <ctime>
#include <stack>
#include <chrono>
#include <random>
#include <string>
#include <cstdlib>
#include <cstring>
#include <cstddef>
#include <fstream>
#include <iostream>
using namespace std::chrono;
// Stringify based macro for traditional variable debugging (prints):
#define varDebug(v) std::cout << #v << " = " << v << "\n"
// Turn off synchronization with standard streams before I/O:
#define fastio std::ios::sync_with_stdio(0); std::cin.tie(0); std::cout.tie(0)
// Template to print multiple variables (comma-separated arguments) of any type each:
template<typename ...T>
void print(T&&... args) 
{
    ((std::cout << args), ...);
}

// Suffix tree node:
struct Node 
{
	Node* A; Node* C; Node* G; Node* T; Node* $;
};

// Information to push onto the stack for traversal:
struct TraversalInformation 
{
	Node* address;
	int mismatchCount;
	int remainingLength;
	char* remainingQuery;	
};

class suffixTree
{
	private:
		Node* root;		
		char* genome;
		int genomeLength, treeNodeCount;
		std::ifstream genomeFile;
		std::string genomeFilePath;		

		// Functions for saving the genome:
		void saveGenome();
		int getLengthOfGenome();		

	public:
	    // Constructors and destructor:
		suffixTree();
		suffixTree(std::string inputGenome);
		~suffixTree();

		// Functions for my suffix tree:
		void createGenomeSuffixTree();
		TraversalInformation* saveTraversalInformation(TraversalInformation* currentTraversalInformation, Node* currentNode, int currentRemainingLength, char* currentQuery, char whichPtr);
		bool searchSuffixTree(char* inputQuery, int inputQueryLength);
		void searchSuffixTreeForRead(char* sequence, int length);
		char** select_N_RandomFragmentsFromGenome(int readCount, int readLength);
		int searchSuffixTrieFor_N_RandomReadsFromGenome(int readCount, int readLength);
};
