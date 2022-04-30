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

// Prefix tree node:
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

// Required class:
class prefixTree
{
	private:
		Node* root;		
		char* genome;
		char** readset;		
		char** genomeKmerArray;
		std::ifstream genomeFile;
		std::string genomeFilePath;		
		int genomeLength, readLength, rowCountForReads, rowCountForGenomeKmerArray, trieNodeCount; // Counter for total number of nodes in the trie (excluding root)
		void splitGenomeIntoKmer(int kmer);		
		void saveNRandomFragmentsFromGenomeAsReadset();

	public:
		// Constructors and destructor:
		prefixTree();                                               
		prefixTree(std::string genomeFilePath, int readCount, int readLength);
		prefixTree(std::string genomeFilePath, int readCount, int readLength, int percentError);
		prefixTree(prefixTree& object);
		~prefixTree();

		// Functions to save the genome and the readset:
		int getLengthOfGenome();
		void saveGenome();		
		void createTrieReadset();
		void introducePercentError(int percentError);

		// Functions to search and traverse the trie:
		int searchTrieForGenomeKmers();		
		int fuzzySearchTrie(char* inputGenomeQuery);		
		TraversalInformation* saveTraversalInformation(TraversalInformation* currentTraversalData, Node* currentNode, int currentMismatchCount, int currentRemainingLength, char* currentQuery, char whichPtr);
};
