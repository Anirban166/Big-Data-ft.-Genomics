/*-------------------------
  Author: Anirban166/Ani
  Email:  ac4743@nau.edu
--------------------------*/
#include "FASTAreadset_LL.h"

// Default constructor:
FASTAreadset_LL::FASTAreadset_LL() 
{
    readCount = 1000;   // Setting an arbitrary default row count for the reads
    readLength = 50;    // Since fragments are always 50 characters long
    // Initialize head and tail pointers to the new node:    
    head = new Node;
    head->next = NULL;
    tail = head;
    for(int i = 0; i < readLength; i++) 
    {
        head->data[i] = '\0';  // Set default data to eof/terminator
    }
}

// Custom constructor: 
// (Takes as input a string which denotes the filepath for the reads dataset)
FASTAreadset_LL::FASTAreadset_LL(std::string file) 
{
    readSetFilePath = file;
    readCount = countNumberOfReads(); // Get number of reads in total file
    readLength = 50;                  // Again, fragments are 50 character long
    head = new Node;
    head->next = NULL;
    tail = head;
    for(int i = 0; i < readLength; i++) 
    {
        head->data[i] = '\0';
    }
}

// Copy constructor: 
// (Takes a FASTAreadset_LL object as input)
FASTAreadset_LL::FASTAreadset_LL(FASTAreadset_LL& object) 
{
    readLength = 50;          // Fragments are always 50 char long
    Node* forward = NULL;    // Pointer to point at the next data (for the original object)
    head = new Node;     
    int count;          

    for(int i = 0; i < readLength; i++) 
    {   // Save the original head into the copied object's head:
        head->data[i] = object.head->data[i];
    }

    current = head;
    forward = object.head->next;
    count = 1; // Offset count since I already saved one node above
    while(forward != NULL) 
    {
        current->next = new Node;
        current = current->next;
        for(int i = 0; i < readLength; i++) 
        {   // Saving original data into the new node:
            current->data[i] = forward->data[i];
        }
        forward = forward->next;
        count++;
    }
    current->next = NULL;
    readCount = count;
    print("\nNumber of reads in the copied object: ", count, "\n");
}

// Destructor:
FASTAreadset_LL::~FASTAreadset_LL() 
{
    while(head->next != NULL) 
    {
        current = head->next;
        delete[] head->data;
        delete head;
        head = current; // Set head to the next element
    }
    delete head;      // Delete head, remaining as the last node at this point
}

// Function to count the total number of reads (or rows) in the input data set file: (or the read set)
// (Returns the number of reads/rows, excluding the headers)
int FASTAreadset_LL::countNumberOfReads()
{
    int count;
    char* temp = new char[1000];          // Allocate space for temporary use (will delete when done)
    count = -1;                          // Offset count due to final "++"
    readSetFile.open(readSetFilePath);
    while(!readSetFile.eof()) 
    { 
        readSetFile >> temp;         // First in order is a header (not counting as a read/row, so no count++ for this)
        readSetFile >> temp;        // Followed by a read fragment/sequence (to be counted)
        count++;                   // Increment count (for the fragment)
    }
    readSetFile.close();
    // print("\nNumber of reads in the data file (or read set): ", count, "\n");
    delete[] temp;
    return count;
}

// Function to search for a string in the readset linked list:
// (Takes as input the 50-character long string query)
// (Returns the Node* pointing to the linked list node that matches the query)
Node* FASTAreadset_LL::searchFragment(std::string query) 
{
    Node* searchQuery = NULL;
    current = head;
    int i;
    while(current != NULL)
    {   
        // I keep counting until either query and node data matches or until I reach the end of the string length:
        i = 0;
        while(current->data[i] == query[i] && i < readLength) 
        {
            i++;
        }
        // If all 50 characters match, then its a success, so I return the pointer to node (as required):
        if(i == 50)
        {
            searchQuery = current;
            return searchQuery;
        }
        current = current->next; // Go to next node in the linked list
    }
    return searchQuery;        // If execution reaches this point, it means no match has been found, and thus I return a NULL pointer.
}

// Function to save the read set as a linked list:
void FASTAreadset_LL::saveReadsLL() 
{
    // Allocate memory (temporarily) for the header and the read:
    char* header = new char[100];
    char* read = new char[100];
    readSetFile.open(readSetFilePath);
    readSetFile >> header;                 // Skip the header
    readSetFile >> read;                  // Save the read
    for(int i = 0; i < readLength; i++)
    {
        head->data[i] = read[i];       // Saving the first read into the head node's data
    }
    for(int i = 1; i < readCount; i++) 
    {
        current = new Node;                
        readSetFile >> header;                 
        readSetFile >> read;                    
        for(int j = 0; j < readLength; j++)
        {
            current->data[j] = read[j];    // Saving the read into the new node's data
        }
        tail->next = current;            // The former last node should now point to the new node
        current->next = NULL;           // And the new node should point to NULL
        tail = current;                // Finally, I update my tail pointer to point at the new node
    }
    delete[] header;
    delete[] read;
    print("\nStatus update: The entire read set file has been saved successfully.\n");  
}

// Function that counts the total number of reads/rows in the genome file:
// (Returns the number of lines excluding the header, i.e. the first 8 lines)
int FASTAreadset_LL::countNumberOfLinesGenome()
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
int FASTAreadset_LL::combineSplitLines(int lineCount, char* read) 
{
    char* tempRead = new char[5000000];   
    char* header = new char[1000000];   
    int length = 0, n;
    genomeFile.open(genomeFilePath);
    for(int i = 0; i < 8; i++) 
    {  
        genomeFile >> header;
    }
    for(int j = 0; j < lineCount; j++) 
    {
        genomeFile >> tempRead;              // Read line from the genome file
        n = 0;
        while(tempRead[n] != '\0') 
        {   
            read[length] = tempRead[n];  // Save characters into a single read
            length++;
            n++;
        }
    }
    genomeFile.close();
    delete[] tempRead;
    delete[] header;
    return length;
}

// Function to save the genome file:
// (Takes as input the string which stores the filepath or the location where the genome file is stored)
// (Post saving, it breaks down into 50-character fragments by shifting start location by 1, with k = 50 or 50-mers)
char** FASTAreadset_LL::saveGenomeFile(std::string filePath) 
{
    char** genome;
    genomeFilePath = filePath;
    char* read = new char[6000000];                      // Temporary placeholder space for the full 1D array genome sequence
    int kmer = 50;                                      // Size of fragments
    int j = 0;                                         // Offset for 1D array genome sequence
    int lineCount = countNumberOfLinesGenome();       // Number of lines in genome file
    int length = combineSplitLines(lineCount, read); // Length of the genome sequence

    // Setting the number of 50-character long fragments to the number of rows in the genome array: (standard formula)
    genomeRowCount = length - kmer + 1;
    print("\nNumber of 50-character long fragments: ", genomeRowCount, "\n");
    
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

// Function to search genome vs readset matches:
// (Takes as input the number of reads and returns the match count)
int FASTAreadset_LL::genomeSearchMatches(int readCount, char** genome)
{     
    int k, matchCount;
    current = head;
    matchCount = 0;
    auto startTime = high_resolution_clock::now();   
    for(int j = 0; j < readCount; j++)             // For the subsection of the readset linked list
    {            
        for(int i = 0; i < genomeRowCount; i++)  // For the entire genome array
        {  
            k = 0;
            while((genome[i][k] == current->data[k]) && (k < 50)) // Keep count of equal characters
            {
                k++;
            }
            if(k == 50) // All characters match - search successful, match found
            {   
                matchCount++;
            }
        }
        current = current->next; // Proceed to the next node (for the next read)
    }
    auto stopTime = high_resolution_clock::now();              
    auto duration = duration_cast<seconds>(stopTime - startTime);
    print("\nSearch time for ", readCount, " reads: ", duration.count(), " seconds\n");
    return matchCount;
}

// Function that prints the results for the query search(es):
// (Takes as input the pointer to match against in the readset linked list and the string query to search for)
void FASTAreadset_LL::printSearchResults(Node* match, std::string query) 
{
    if(match == NULL) 
        print("\nQuery sequence ", query, " was not found in the linked list.\n");
    else 
    {
        print("\nQuery sequence ");
        for(int i = 0; i < 50; i++) // Print the entire 50-character long sequence
        {
            print(query[i]);
        } 
        print(" was found in the linked list, with the returned pointer data being ");
        for(int i = 0; i < 50; i++) 
        {
            print(match->data[i]);
        }
        print(".\n");
    }
}