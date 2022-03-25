/*-------------------------
  Author: Anirban
  Email:  ac4743@nau.edu
--------------------------*/
#include "FASTAreadset_HT.h"

// Constructor:
FASTAreadset_HT::FASTAreadset_HT(std::string filePath, unsigned int hashTableSize) 
{
    m = hashTableSize;        // Save input hash table array size
    base = 4;                // 4 character options (A, C, G, T)
    sequenceLength = 16;    // Length of fragment is 16 (16-mers) 
    genomeFilePath = filePath;    
    genomeRowCount = 0; collisions = 0; elementCount = 0; duplicates = 0;
    // Number of possible radix numbers (size) is given by base^sequenceLength (i.e. 4^16):
    size = (unsigned long long)(pow(base, sequenceLength));
    print("Size of the hash table is: ", m , "\n");
    // Allocate memory for the array of linked lists:
    T = new Node * [m];
    for(unsigned int i = 0; i < m; i++) 
    {
        T[i] = new Node;
    }
    // Setting default values for the first element: (data set to 1 so as to avoid confusion with the data value 0)
    T[0]->data = 1;
    T[0]->next = NULL;
    for(unsigned int i = 1; i < m; i++) 
    {
        T[i]->data = 0;
        T[i]->next = NULL;
    }
}

// Destructor:
FASTAreadset_HT::~FASTAreadset_HT() 
{
    Node* tempHead;
    Node* current;
    // Loop through the array and delete each pointer to my linked list:
    for(unsigned int i = 0; i < m; i++) 
    {
        tempHead = T[i];
        while(tempHead->next != NULL) 
        {
            current = tempHead->next;  // Save next element in the second pointer (wrt head)
            delete tempHead;
            tempHead = current;      // Set head to the next element
        }
        if(tempHead->next == NULL) 
        {
            delete tempHead;
        }
    }
    delete[] T;
}

// Function to convert 16-char array with characters A, C, G and T to 16-int array with integers 0, 1, 2 and 3:
int* FASTAreadset_HT::characterToIntegerConversion(char* readCharacter, int* readInteger) 
{
    for(int i = 0; i < sequenceLength; i++) 
    {
        switch(readCharacter[i])
        {
            case 'A': readInteger[i] = 0; break;
            case 'C': readInteger[i] = 1; break;     
            case 'G': readInteger[i] = 2; break;
            case 'T': readInteger[i] = 3; break;                                
             default: print("Invalid character, must be one of A,C,G,T. Exiting!\n");
                      exit(-1);
        }
    }
    return readInteger;
}

// Function to convert the integer sequence array into a radix number:
unsigned int FASTAreadset_HT::radixBaseConversion(int* readInteger) 
{
    int exponent;
    unsigned int radixNumber = 0;
    // Radix number = (value[0])(base^sequenceLength - 1) + (value[1])(base^sequenceLength - 2) + ...
    for(int i = 0; i < sequenceLength; i++) 
    {
        exponent = sequenceLength - i - 1;
        radixNumber += readInteger[i] * (unsigned int)(pow(base, exponent));
    }
    return radixNumber;
}

// Function to convert sequence into radix number:
unsigned int FASTAreadset_HT::convertSequenceToRadixNumber(char* inputSequence) 
{
    char* readCharacter = new char[sequenceLength];
    int* readInteger = new int[sequenceLength];    
    unsigned int radixNumber;
    // Saving input fragments:
    for(int j = 0; j < sequenceLength; j++) 
    {
        readCharacter[j] = inputSequence[j];
    }
    readInteger = characterToIntegerConversion(readCharacter, readInteger);
    radixNumber = radixBaseConversion(readInteger);
    delete[] readCharacter;
    delete[] readInteger;
    return radixNumber;
}

// Function to convert radix numbers into addresses using the division method:
unsigned int FASTAreadset_HT::divisionMethodHash(unsigned int radixNumber) 
{
    unsigned int tableAddress;
    tableAddress = radixNumber % m;
    return tableAddress;
}

// Function to save the genome file into the hash table:
void FASTAreadset_HT::saveGenomeFileInHashTable() 
{   
    auto startTime = high_resolution_clock::now();
    // Declaring (and allocating memory for) arrays to store the temporary input, radix number and address of all reads:
    unsigned int* radixArray = new unsigned int[genomeRowCount];
    unsigned int* addressArray = new unsigned int[genomeRowCount];
    char* tempFragment = new char[100];
    genomeFile.open(genomeFilePath);
    // Converting to radix, applying division method and inserting in the hash table for all reads:
    for(unsigned long long i = 0; i < genomeRowCount; i++) 
    {
        genomeFile >> tempFragment; // Skip Header
        genomeFile >> tempFragment;
        radixArray[i] = convertSequenceToRadixNumber(tempFragment);
        addressArray[i] = divisionMethodHash(radixArray[i]);
        insertInTable(addressArray[i], radixArray[i]);
    }
    genomeFile.close();
    delete[] tempFragment;
    delete[] radixArray;
    delete[] addressArray;
    auto stopTime = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stopTime - startTime);
    print("Time for saving the genome file to hash table is: ", duration.count(), " seconds\n");
    return;
}

// Function to save the genome file:
// (Takes as input the string which stores the filepath or the location where the genome file is stored)
// (Post saving, it breaks down into k-character fragments by shifting start location by 1, with k = 16 or 16-mers)
char** FASTAreadset_HT::saveGenomeFile(std::string filePath) 
{
    char** genome;
    genomeFilePath = filePath;
    char* read = new char[6000000];                      // Temporary placeholder space for the full 1D array genome sequence
    int kmer = 16;                                      // Size of sequence
    int j = 0;                                         // Offset for 1D array genome sequence
    int lineCount = countNumberOfLinesGenome();       // Number of lines in genome file
    int length = combineSplitLines(lineCount, read); // Length of the genome sequence
    // Setting the number of k-character long fragments to the number of rows in the genome array: (standard formula)    
    genomeRowCount = length - kmer + 1;
    print("Number of ", kmer, "-char long fragments: ", genomeRowCount, "\n");
    // Initializing the genome array, now that I know the size of it:
    genome = new char* [genomeRowCount];
    for(int i = 0; i < genomeRowCount; i++) 
    {
        genome[i] = new char[kmer];
    }
    // Split into k-mer fragments and save as a 2D character array:
    for(int r = 0; r < genomeRowCount;r++) 
    {
        for(int i = 0; i < kmer; i++) 
        {
            genome[r][i] = read[i + j];
        }
        j++;
    }
    delete[] read;
    return genome;
}

// Function to insert radix number at address:
void FASTAreadset_HT::insertInTable(unsigned int address, unsigned int radixNumber)
{
    bool isThere;                 // Boolean indicating whether the radix number is in table or not
    Node* current = T[address];  // Temporary pointer to be made to point to the node of interest
    // Checking if the radix number is already in the table:
    isThere = searchTableForRadix(T[address], radixNumber);
    if(isThere) // If it is there, then its a duplicate!
    {
        duplicates++;
        elementCount++;
        return;
    }
    else 
    {
        if((current->data) % m != radixNumber % m)
        {
            current->data = radixNumber; // Unique value
            elementCount++;
            return;
        }
        else if((current->data) % m == radixNumber % m) // Collision
        {
            T[address] = new Node;                  
            T[address]->data = radixNumber;  
            T[address]->next = current;            // Adding to the front of the list
            elementCount++;
            collisions++;
            return;
        }
        else 
        {
            print("Error in table insertion, exiting the program!\n");
            exit(-1);
        }
    }
}

// Function to search the table for radix number:
bool FASTAreadset_HT::searchTableForRadix(Node* tableAddress, unsigned int radixNumber) 
{
    Node* current;
    current = tableAddress;
    if(current->data == radixNumber) // First item (radix number found) case
    {     
        return true;                  
    }
    else 
    {
        if(current->next == NULL) // No items (radix number not found) case
        {
            return false; 
        }
        else if(current->next != NULL) 
        {
            while(current->next != NULL) 
            { 
                current = current->next;
                if(current->data == radixNumber) 
                {
                    return true;
                }
            }
            return false;
        }
        else 
        {
            print("Error in searching for a radix number, exiting!\n");
            exit(-1);
        }
    }
}

// Function to search table for given sequence:
bool FASTAreadset_HT::searchTableForSequence(char* sequence)
{
    unsigned int radixNumber, address;
    bool presence;
    radixNumber = convertSequenceToRadixNumber(sequence);
    address = divisionMethodHash(radixNumber);
    presence = searchTableForRadix(T[address], radixNumber);
    return presence;
}

// Function to insert given sequence into the table:
void FASTAreadset_HT::insertSequenceInTable(char* sequence) 
{
    unsigned int radixNumber, address;
    radixNumber = convertSequenceToRadixNumber(sequence);
    address = divisionMethodHash(radixNumber);
    insertInTable(address, radixNumber);
    print("Inserted the sequence ");
    for(int i = 0; i < sequenceLength; i++) 
    {
        print(sequence[i]);
    }
    print("!\n");
    return;
}

// Function to print number of elements (includes all attempts to add a fragment) added to the table:
void FASTAreadset_HT::printNumberOfTableElements() 
{
    print("Number of elements added: ", elementCount, "\n");
    return;
}

// Function to print the number of collisions: (these collisions include fragments with different radix #s but same table address value)
void FASTAreadset_HT::printNumberOfTableCollisions() 
{
    print("Number of collisions: ", collisions, "\n");
    return;
}

// Function to print the number of dupes: (includes fragments that are identical with the same radix # and therefore the same address too)
void FASTAreadset_HT::printNumberOfTableDuplicates() 
{
    print("Number of duplicates: ", duplicates, "\n");
    return;
}

// Function to print the number of stored elements:
// (Stored elements include collisions and unique values, but not duplicates since dupes are not saved in the list again)
void FASTAreadset_HT::printNumberOfStoredValuesInTable() 
{
    unsigned long long storedValueCount = elementCount - duplicates;
    print("Number of stored values:  ", storedValueCount, "\n");
    return;
}

// Function that counts the total number of rows in the genome file:
// (Returns the number of lines excluding the header, i.e. the first 8 lines)
int FASTAreadset_HT::countNumberOfLinesGenome()
{
    int count;
    char* temp = new char[5000000]; 
    count = -1;                    
    genomeFile.open(genomeFilePath);
    for(int i = 0; i < 8; i++)   // Skip the header (first 8 lines)
    {   
        genomeFile >> temp;
    }
    while(!genomeFile.eof()) // Count lines until the end of line
    {
        genomeFile >> temp;
        count++;
    }
    genomeFile.close();
    delete[] temp;
    print("\nTotal number of rows in the genome file: ", count, "\n");
    return count;
} 

// Function that combines split lines of the genome and returns the resultant length of it:
// (Takes as input the number of lines in the genome file and the location to save the split lines as an array)
// (Returns the number of elements in the combined 1D array, which is obtained from the separate split lines)
int FASTAreadset_HT::combineSplitLines(int lineCount, char* sequence) 
{
    char* tempSequence = new char[5000000];   
    char* header = new char[1000000];   
    int length = 0, n;
    genomeFile.open(genomeFilePath);
    for(int i = 0; i < 8; i++) 
    {  
        genomeFile >> header;
    }
    for(int j = 0; j < lineCount; j++) 
    {
        genomeFile >> tempSequence;                  // Read line from the genome file
        n = 0;
        while(tempSequence[n] != '\0') 
        {   
            sequence[length] = tempSequence[n];  // Save characters into a single read
            length++;
            n++;
        }
    }
    genomeFile.close();
    delete[] tempSequence;
    delete[] header;
    return length;
}

// Function to search the genome for matches against seqeunces in the hash table: 
unsigned int FASTAreadset_HT::genomeSearch(char** genome) 
{
    bool presence;
    unsigned int matchCount = 0;
    for(int i = 0; i < genomeRowCount; i++) 
    { 
        presence = searchTableForSequence(genome[i]);
        if(presence) matchCount++;
    }   
    return matchCount;
}

// Function to search the genome for matches against seqeunces in the hash table with 1% chance of every character to be changed to another one: 
unsigned int FASTAreadset_HT::genomeSearchWithOnePercentBaseErrorRate(char** genome) 
{   // Modify characters based on 1% per-base error rate (99% chance of no modification) before searching:
    for(int i = 0; i < genomeRowCount; i++) // Rows
    {
        for(int j = 0; j < 16; j++)       // Characters
        {        
            // 1% chance of getting 0 among 0-99: (i.e. 1 among 100 possible values)
            int lucky = rand() % 100; // Offset is 0
            // std::uniform_int_distribution<int> distribution(0, 99);
            if(!lucky) // if 0, then change the current character!
            {
                int randomLetter = rand() % 3; 
                if(genome[i][j] == 'A' && randomLetter == 0)      genome[i][j] = 'C';
                else if(genome[i][j] == 'A' && randomLetter == 1) genome[i][j] = 'G';
                else if(genome[i][j] == 'A' && randomLetter == 2) genome[i][j] = 'T';
                else if(genome[i][j] == 'C' && randomLetter == 0) genome[i][j] = 'A';
                else if(genome[i][j] == 'C' && randomLetter == 1) genome[i][j] = 'G';
                else if(genome[i][j] == 'C' && randomLetter == 2) genome[i][j] = 'T';
                else if(genome[i][j] == 'G' && randomLetter == 0) genome[i][j] = 'A';
                else if(genome[i][j] == 'G' && randomLetter == 1) genome[i][j] = 'C';
                else if(genome[i][j] == 'G' && randomLetter == 2) genome[i][j] = 'T';                                
                else if(genome[i][j] == 'T' && randomLetter == 0) genome[i][j] = 'A';
                else if(genome[i][j] == 'T' && randomLetter == 1) genome[i][j] = 'C';
                else if(genome[i][j] == 'T' && randomLetter == 2) genome[i][j] = 'G'; 
            }
        }
    }
    bool presence;
    unsigned int matchCount = 0;
    for(int i = 0; i < genomeRowCount; i++) 
    { 
        presence = searchTableForSequence(genome[i]);
        if(presence) matchCount++;
    }   
    return matchCount;
}

// Function to take random 16-mer sequences from the genome and then compare those against the hash table stored sequences:
unsigned int FASTAreadset_HT::randomGenomeSequenceComparison(char** genome)
{
    char** tempRandomGenome = nullptr;
    // Generating all the characters randomly for five million 16-mer sequences:    
    for(int i = 0; i < 5000000; i++)
    {
        for(int j = 0; j < 16; j++)
        {        
            int startRow = rand() % 5000000; // Randomly pick any one among one million rows of the original genome to start with
            tempRandomGenome[i][j] = genome[startRow][j];           
        }
    }
    bool presence;
    unsigned int matchCount = 0;
    for(int i = 0; i < 5000000; i++) 
    { 
        presence = searchTableForSequence(tempRandomGenome[i]);
        if(presence) matchCount++;
    }   
    return matchCount; 
}  

// Function to make completely random 16-mer sequences and then compare elements against the hash table stored genome sequences: 
unsigned int FASTAreadset_HT::randomSequenceComparison() 
{   
    char** tempRandomGenome = nullptr;
    // Generating all the characters randomly for five million 16-mer sequences:    
    for(int i = 0; i < 5000000; i++)
    {
        for(int j = 0; j < 16; j++)
        {        
            int randomLetter = rand() % 4;
            if(!randomLetter)          tempRandomGenome[i][j] = 'A';
            else if(randomLetter == 1) tempRandomGenome[i][j] = 'C';
            else if(randomLetter == 2) tempRandomGenome[i][j] = 'G';
            else if(randomLetter == 3) tempRandomGenome[i][j] = 'T';            
        }
    }
    bool presence;
    unsigned int matchCount = 0;
    for(int i = 0; i < 5000000; i++) 
    { 
        presence = searchTableForSequence(tempRandomGenome[i]);
        if(presence) matchCount++;
    }   
    return matchCount;    
}      