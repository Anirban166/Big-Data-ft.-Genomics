/*-------------------------
  Author: Anirban
  Email:  ac4743@nau.edu
--------------------------*/
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
using namespace std::chrono;
// Turn off synchronization with standard streams before I/O:
#define fastio std::ios::sync_with_stdio(0); std::cin.tie(0); std::cout.tie(0)
// Stringify based macro for traditional variable debugging (prints):
#define varDebug(v) std::cout << #v << " = " << v // << "\n"
// Template to print multiple variables (comma-separated arguments) of any type each:
template<typename ...T>
void print(T&&... args) 
{
    ((std::cout << args), ...);
}

// Required class:
class FASTA_readset 
{
	private:
	std::ifstream dataFile;                                                      // Variable to store the input file or data set
    	char* filePath = new char[100];                                             // Char array to hold the file path (should not exceed 100 characters)
    	int* uniqueSequenceArray = new int[14];                                    // Arrays (one for each data set) to hold the count of unique sequences
    	int* totalSequenceArray = new int[14];                                    // Arrays (one for each data set) to hold the count of sequences (total)
    	char** readArray; int** countArray;                                      // Array of arrays to hold the reads and counts (for the entire data set)
    	char** dataSetOneReadArray; char** dataSetTwoReadArray;                 // Char array of arrays to hold the contents of data sets one and two
        int row, countCol, readCol;                                            // Variables for holding the row count and the read/count-based col counts
    	int countA, countC, countG, countT, countN;                           // Variables for holding the letter frequency counts (extra)
    	int countNumberOfReads();                                            // Function that returns the total number of reads (lines/2) in the data set
    	int getHeaderLength(char* header);                                  // Function that returns the length of the header string
     	void swap(int k, int j);                                           // Function to swap two rows based on indices (to be used in my sort function)
    	void saveReadArray (char* read, int index);                       // Function that saves the read array
    	void saveCountArray (char* header, int headerLength, int index); // Function that saves the count array

	public:
	FASTA_readset();                                     // Default constructor
	FASTA_readset(int rowCount, char* filePath);        // Custom constructor with a custom row count and filepath as input:
        ~FASTA_readset();	                           // Destructor
        void saveData();                                  // Function to save the array data after parsing (to be called first!)
        void sortReadArray();                            // Function to sort the sequence fragment rows for all the reads (quadratic)
        void sortSecondDataSetReadArray();              // Function to sort the sequence fragments or reads for the second dataset (quadratic)
        int searchFunction(char** a, char* b);         // Function to search for a particular read among reads in a data set (linear)
        int binarySearchFunction(char** a, char* b);  // Function to search for a particular sequence fragment (row) among reads in a data set (log-linear)
        void computeUniqueSequenceCount();           // Function to compute the count of unique sequences
        void computeTotalSequenceCount();           // Function to compute the count of sequences (total)
        void countLetterFrequency();               // Function (extra) to compute the frequencies/occurences of the letters A, C, G, T and N
        void computeStatistics(char problemID);   // Function to compute the statistics for a particular problemID (A, B, C), which is 9th character in argv[2]
        void printCountArray();                  // Function to print the count array
        void printReadArray();                  // Function to print the read array        
        void printUniqueSequenceCount();       // Function to print the count of unique sequences for each data set
        void printTotalSequenceCount();       // Function to print the count of sequences (total) for each data set
        void printLetterCounts();            // Function (extra) to print the frequencies of A, C, G, T and N
        void printDataSetOne();             // Function to print dataset one
        void printDataSetTwo();            // Function to print dataset two
        void createDataSetOneAndTwo();    // Function to create datasets one and two, extracting dataset-specific rows from the read array
};