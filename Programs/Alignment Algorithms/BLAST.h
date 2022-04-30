/*-------------------------
  Author: Anirban166/Ani
  Email:  ac4743@nau.edu
--------------------------*/
#pragma once
#include <cmath>
#include <ctime>
#include <chrono>
#include <random>
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

// Struct for the hash table: (with collision-based entries chained in a link list)
struct Node 
{
	int data;
	int position;
	Node* next;
};

// Required class:
class BLAST
{
	public:
		int readLength, matchScore, gapScore, mismatchScore;
		BLAST(int scoreForMatch, int penalityForMismatch, int penaltyForGap);

		// Functions dealing with the readset:
		int countNumberOfReads(std::string readFilePath);
		void saveReads(int readCount, int readLength, char** readset, std::string readFilePath);
		void createRandomReads(char** readset, int readCount, int readLength);
		char** selectNRandomFragmentsFromGenome(int readCount, int readLength, char* genome, int genomeLength);
		void introducePercentError(int readCount, int readLength, char** readset, int percentError);
		char** splitFragmentIntoKmer(char* input, int inputLength, int kmer, int* rowCount);		

		// Functions dealing with the genome:
		int getLengthOfGenome(std::string genomeFilePath);
		void saveGenome(std::string genomeFilePath, char* genome, int genomeLength);

		// Functions for alignment:
		void computeArray(int** D, char** traceback, int matchScore, int gapScore, int mismatchScore, char* sequenceOne, char* sequenceTwo, int sequenceOneLength, int sequenceTwoLength);
		textAlignment* BLASTAlignment(char* sequenceOne, char* sequenceTwo, int sequenceOneLength, int sequenceTwoLength, int matchScore, int gapScore, int mismatchScore, int* max);
		int findMaxInArray(int** D, int* indexOne, int* indexTwo, int sequenceOneLength, int sequenceTwoLength);
		int performTraceback(int** D, char** traceback, int* indexOne, int* indexTwo, char* tempSequenceOne, char* tempSequenceTwo, char* tempCode, char* sequenceOne, char* sequenceTwo);
		int getSubGenomeLength(int* lowBound, int* upBound, int* genomeIndex, int currentReadKmer, int kmer, int currentReverseIndex, int genomeLength);
		void expandSeed(int* lowBound, int* upBound, char* genome, char* subGenome);
		char* searchAndExpandSeed(int readNumber, char** readset, int readLength, int kmer, int hashTableSize, int* genomeIndex, Node** T, int genomeLength, char* genome, int* subGenomeLength, bool* seedFound, char* subGenome, bool printFlag);
		void performBLASTForNRandomReads(int readCount, int readLength, int kmer, int hashTableSize, Node** T, char* genome, int genomeLength, int* genomeIndex, int matchScore, int gapScore, int mismatchScore, int* max, bool printFlag);

		// Functions for my hash table:
		Node** createHashTable(int hashTableSize);
		void insertFragmentInTable(char* fragment, int kmer, int hashTableSize, Node** T, int* genomeIndex);
		int convertFragmentToRadixNumber(char* inputRead, int kmer);
		int* characterToIntegerConversion(char* readCharacter, int* readInteger, int kmer);
		int radixBaseConversion(int* readInteger, int kmer);
		int divisionMethodHash(int radixNumber, int hashTableSize);
		void insertInTable(int address, int radixNumber, Node** T, int hashTableSize, int* genomeIndex);
		bool searchTableForRadix(Node* tableAddress, int radixNumber, int* genomeIndex);
		bool searchTableForFragment(char* fragment, int kmer, int hashTableSize, Node** T, int* genomeIndex);
		void deleteHashTable(int hashTableSize, Node** T);
};

// Class for text alignment:
class textAlignment
{
	public:
		int length;
		char* code;
		char* sequenceOneText;
		char* sequenceTwoText;
		textAlignment(int len);
		void printAlignment(int alignmentLength);
		void saveAlignment(int alignmentLength, char* sequenceOne, char* sequenceTwo, char* alignmentCode);
};