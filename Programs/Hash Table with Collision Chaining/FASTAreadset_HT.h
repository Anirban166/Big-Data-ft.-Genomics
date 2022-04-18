/*-------------------------
  Author: Anirban166/Ani
  Email:  ac4743@nau.edu
--------------------------*/
#pragma once
#include <cmath>
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

// Required node structure, composed of a positive integer for the data and the usual pointer to the next node:
struct Node 
{
	unsigned int data;
	Node* next;
};

// Required class:
class FASTAreadset_HT
{
    	private:
	        unsigned int m;                       // Variable to denote the size of the hash table
	        Node** T;                            // Pointer to array of pointers to linked lists
		unsigned long long size;            // Largest possible radix number value
		unsigned long long elementCount;   // Number of elements added to hash table
		unsigned long long collisions;    // Number of sequences with the same address but different radix number
		unsigned long long duplicates;   // Number of sequences with same radix number
		std::string genomeFilePath;     // Variable to store the genome file location
		std::ifstream genomeFile;      // Variable to store the genome file itself
		int base, sequenceLength;     // Variables to store the length of the base and the sequence resepectively (4, 16)
	
	public:
		int genomeRowCount;        // Variable to account for the number of rows in the genome array (i.e., number of k-mer fragments)
		FASTAreadset_HT(std::string filePath, unsigned int hashTableSize);                    // Constructor
		~FASTAreadset_HT();                                                                  // Destructor
		void insertSequenceInTable(char* sequence);                                         // Function to insert a given n-mer sequence into the hash table 		
		bool searchTableForSequence(char* sequence);                                       // Function to search the hash table for a given n-mer sequence
		unsigned int radixBaseConversion(int* readInt);                                   // Function to convert the integer sequence array into a radix number
		unsigned int divisionMethodHash(unsigned int radixNumber);                       // Function to convert radix numbers into addresses using the division method
		unsigned int convertSequenceToRadixNumber(char* inputSequence);		        // Function to convert a given sequence into a radix number notation
		void insertInTable(unsigned int address, unsigned int radixNumber);	       // Function to insert radix number at a provided address
		bool searchTableForRadix(Node* tableAddress, unsigned int radixNumber);       // Function to search the table for a radix number
		int* characterToIntegerConversion(char* readCharacter, int* readInteger);    // Function to convert 16-char array with characters A, C, G and T to integers (0, 1, 2, 3)		
		void printNumberOfTableElements();                                          // Function to print number of elements (includes all attempts to add a fragment) added to the table
		void printNumberOfTableCollisions();                                       // Function to print the number of collisions (these collisions include fragments with different radix #s but same table address value)
		void printNumberOfTableDuplicates();	                                  // Function to print the number of dupes (includes fragments that are identical with the same radix # and therefore the same address too)		
		void printNumberOfStoredValuesInTable();                                 // Function to print the number of stored values (total elements minus duplicates)
		int countNumberOfLinesGenome();                                         // Function that counts the total number of rows in the genome file
		void saveGenomeFileInHashTable();                                      // Function to save the genome file into the hash table
		unsigned int genomeSearch(char** genome);                             // Function to search the genome for matches against seqeunces from the hash table (iterating through the genome's 16-mers)
        	unsigned int genomeSearchWithOnePercentBaseErrorRate(char** genome); // Function to search the genome for matches against seqeunces in the hash table with 1% chance of every character to be changed to another one!	
		unsigned int randomSequenceComparison();                            // Function to make completely random 16-mer sequences and then compare those against the hash table stored genome sequences
		unsigned int randomGenomeSequenceComparison(char** genome);        // Function to take random 16-mer sequences from the genome and then compare those against the hash table stored sequences
		char** saveGenomeFile(std::string genomeFilePath);                // Function to save the genome file
		int combineSplitLines(int numberOfLines, char* sequence);        // Function that combines split lines of the genome (this, and the above are helper functions for saving the genome file)				
};
