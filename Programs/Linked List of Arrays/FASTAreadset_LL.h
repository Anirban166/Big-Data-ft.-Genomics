/*-------------------------
  Author: Anirban166/Ani
  Email:  ac4743@nau.edu
--------------------------*/
#pragma once
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <cstddef>
#include <fstream>
#include <iostream>
#include "FASTAreadset_LL.h"
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

// Required node structure, with 50-character long data and the usual pointer to the next node:
struct Node 
{
    char* data = new char[50];
    Node* next;
};

// Required class:
class FASTAreadset_LL
{
    private:
        std::ifstream readSetFile, genomeFile;                        // Variables to store the readset and genome files
        std::string readSetFilePath, genomeFilePath;                 // Variables for holding the readset and genome file locations
        Node* head = NULL; Node* tail = NULL; Node* current = NULL; // Pointers for my linked list (initialized to NULL)
        int readLength, readCount;                                 // Integer variables to store the length of the read fragments and the number of reads in the dataset
        int countNumberOfLinesGenome();                           // Function that counts the total number of reads in the genome file
        int combineSplitLines(int lineCount, char* read);        // Function that combines split lines of the genome (this, and the above are helper functions for saving the genome file)

    public:
        int genomeRowCount;                                                // Variable to account for the number of rows in the genome array
        FASTAreadset_LL();                                                // Default constructor
        FASTAreadset_LL(std::string file);                               // Custom constructor
        FASTAreadset_LL(FASTAreadset_LL& object);                       // Copy constructor
        ~FASTAreadset_LL();                                            // Destructor
        void saveReadsLL();                                           // Function to save the readset as a linked list
        int countNumberOfReads();                                    // Function to count the total number of reads in the input data set file
        Node* searchFragment(std::string query);                    // Function to search for a string in the readset linked list
        char** saveGenomeFile(std::string filePath);               // Function to save the genome file
        int genomeSearchMatches(int readCount, char** genome);    // Function to search for genome vs readset matches
        void printSearchResults(Node* match, std::string query); // Function that prints the results for the query searches
};