/*-------------------------
  Author: Anirban166/Ani
  Email:  ac4743@nau.edu
--------------------------*/
#include "suffixTree.h"

// Default constructor:
suffixTree::suffixTree() 
{
    treeNodeCount = 0;
    genomeLength = 50;
    genome = new char[genomeLength];
    saveGenome();
    root = new Node;
    root->A = NULL; root->C = NULL; root->G = NULL; root->T = NULL; root->$ = NULL;
}

// Custom constructor:
suffixTree::suffixTree(std::string inputGenome) 
{
    treeNodeCount = 0;
    genomeFilePath = inputGenome;
    genomeLength = getLengthOfGenome();
    genome = new char[genomeLength];
    saveGenome();
    root = new Node;
    root->A = NULL; root->C = NULL; root->G = NULL; root->T = NULL; root->$ = NULL;
}

// Destructor:
suffixTree::~suffixTree() 
{
    Node* currentNode;
    std::stack <Node*> deconStack;

    // Pushing all the non-NULL children of the root to my stack:
    currentNode = root;
    if(currentNode->A != NULL) deconStack.push(currentNode->A);
    if(currentNode->C != NULL) deconStack.push(currentNode->C);
    if(currentNode->G != NULL) deconStack.push(currentNode->G);
    if(currentNode->T != NULL) deconStack.push(currentNode->T);
    if(currentNode->$ != NULL) deconStack.push(currentNode->$);

    // Popping, deleting, and then pushing all non-NULL children: (repeated until the stack is empty)
    while(!deconStack.empty()) 
    {
        root = (deconStack.top());
        deconStack.pop();
        delete currentNode;
        currentNode = root;
        if(currentNode->A != NULL) deconStack.push(currentNode->A);
        if(currentNode->C != NULL) deconStack.push(currentNode->C);
        if(currentNode->G != NULL) deconStack.push(currentNode->G);
        if(currentNode->T != NULL) deconStack.push(currentNode->T);
        if(currentNode->$ != NULL) deconStack.push(currentNode->$);
    }
    delete[] genome;
    return;
}

// Function to return the number of characters in the genome:
int suffixTree::getLengthOfGenome() 
{
    char* tempRead = new char[5000000];
    char* header = new char[100];
    int genomeLength = 0;
    
    genomeFile.open(genomeFilePath);
    genomeFile >> header;             
    genomeFile >> tempRead;

    while(tempRead[genomeLength] != '\0') 
    {
        genomeLength++;
    }

    genomeFile.close();
    delete[] tempRead;
    delete[] header;
    return genomeLength;
}

// Function to save the genome file as a character array:
void suffixTree::saveGenome() 
{
    char* tempRead = new char[2000000];
    char* header = new char[100];

    genomeFile.open(genomeFilePath);
    genomeFile >> header;
    genomeFile >> tempRead;

    for(int i = 0; i < genomeLength - 1; i++) 
    {
        genome[i] = tempRead[i];
    }
    genome[genomeLength - 1] = '$'; // Adding the terminating character

    genomeFile.close();
    delete[] tempRead;
    delete[] header;
    return;
}

// Function to save the genome as a suffix tree:
void suffixTree::createGenomeSuffixTree() 
{
    char currentCharacter;
    Node* currentNode;
    for(int i = 0; i < genomeLength; i++) 
    {   // For each suffix in the genome (number of suffixes == number of characters in genome)
        currentNode = root;
        for(int j = 0; j < genomeLength - i; j++) 
        {   // For each character in current suffix:
            currentCharacter = genome[i + j]; // Current character to insert is genome[i + j] == first character in the current suffix
            switch(currentCharacter) 
            {
                #define nullAllNodes() \
                    currentNode->A = NULL; \
                    currentNode->C = NULL; \
                    currentNode->G = NULL; \
                    currentNode->T = NULL; \
                    currentNode->$ = NULL

                case 'A':
                if(currentNode->A == NULL) 
                {
                    treeNodeCount++;
                    currentNode->A = new Node;
                    currentNode = currentNode->A;
                    nullAllNodes();
                }
                else currentNode = currentNode->A;
                break;

                case 'C':                    
                if(currentNode->C == NULL) 
                {  
                    treeNodeCount++;
                    currentNode->C = new Node;
                    currentNode = currentNode->C;
                    nullAllNodes();
                }
                else currentNode = currentNode->C;
                break;

                case 'G':                  
                if(currentNode->G == NULL) 
                {
                    treeNodeCount++;
                    currentNode->G = new Node;
                    currentNode = currentNode->G;
                    nullAllNodes();
                }
                else currentNode = currentNode->G;
                break;

                case 'T':
                if(currentNode->T == NULL) 
                { 
                    treeNodeCount++;
                    currentNode->T = new Node;
                    currentNode = currentNode->T;
                    nullAllNodes();
                }
                else currentNode = currentNode->T;
                break;

                case '$':              
                if(currentNode->$ == NULL) 
                {
                    treeNodeCount++;
                    currentNode->$ = new Node;
                    currentNode = currentNode->$;
                    nullAllNodes();
                }
                else currentNode = currentNode->$;
                break;

                default:
                print("Error! Valid characters include A, C, G, T (also $) only.\n");
                exit(-1);
            }
        }
    }
    print("Number of nodes: ", treeNodeCount, "\n");
    return;
}

// Function to save traversal information in a struct:
// (Saving the node address, the number of mismatches, the remaining query string and the length of it in a struct whose purpose is to get pushed to the stack during the trie search)
TraversalInformation* suffixTree::saveTraversalInformation(TraversalInformation* currentTraversalInformation, Node* currentNode, int currentRemainingLength, char* currentQuery, char whichPtr) 
{
    // Determine which child address to save:
    if(whichPtr == 'A')      currentTraversalInformation->address = currentNode->A;
    else if(whichPtr == 'C') currentTraversalInformation->address = currentNode->C;
    else if(whichPtr == 'G') currentTraversalInformation->address = currentNode->G;
    else if(whichPtr == 'T') currentTraversalInformation->address = currentNode->T;
    else 
    {
        print("Error! Valid characters include A, C, G and T only.\n");
        exit(-1);
    }
    // Save the remaining information:
    currentTraversalInformation->remainingLength = currentRemainingLength;
    currentTraversalInformation->remainingQuery = new char[currentRemainingLength];
    for(int i = 0; i < currentRemainingLength; i++) 
    {
        currentTraversalInformation->remainingQuery[i] = currentQuery[i + 1];
    }
    return currentTraversalInformation;
}

// Function to search the suffix tree: (for a given character array of a given length; returns boolean true if found, else false)
bool suffixTree::searchSuffixTree(char* inputQuery, int inputQueryLength) 
{
    Node* currentNode = root;                            // Node pointer to track where I am on the tree
    std::stack <TraversalInformation*> myStack;         // Stack ftw :)    
    TraversalInformation* currentTraversalInformation; // Pointer to struct that holds the information pushed onto the stack
    char* currentQuery = inputQuery;                  // Stores the current query character array
    char currentQueryCharacter = inputQuery[0];      // Current query character to look for (starting with the first character of the input query)
    int currentRemainingLength = 0;                 // The remaining query characters to search
    int characterMatchCount = 0;                   // Number of characters in the input query found in suffix tree

    // Obtaining the remaining query length:
    for(int i = 1; i < inputQueryLength; i++) 
    {
        currentRemainingLength++;
    }

    // Pushing the non-NULL children of the root with traversal information onto my stack:
    switch(currentQueryCharacter) 
    {
        case 'A': // If character matches with the base A, I save and push the traversal information if it ain't NULL
        if(currentNode->A != NULL) 
        {
            currentTraversalInformation = new TraversalInformation;
            currentTraversalInformation = saveTraversalInformation(currentTraversalInformation, currentNode, currentRemainingLength, currentQuery, 'A');
            myStack.push(currentTraversalInformation);
        }
        else return false;
        break;

        case 'C':
        if(currentNode->C != NULL) 
        {
            characterMatchCount++;
            currentTraversalInformation = new TraversalInformation;
            currentTraversalInformation = saveTraversalInformation(currentTraversalInformation, currentNode, currentRemainingLength, currentQuery, 'C');
            myStack.push(currentTraversalInformation);
        }
        else return false;
        break;

        case 'G':
        if(currentNode->G != NULL) 
        {
            characterMatchCount++;
            currentTraversalInformation = new TraversalInformation;
            currentTraversalInformation = saveTraversalInformation(currentTraversalInformation, currentNode, currentRemainingLength, currentQuery, 'G');
            myStack.push(currentTraversalInformation);
        }
        else return false;
        break;

        case 'T':
        if(currentNode->T != NULL) 
        {
            characterMatchCount++;
            currentTraversalInformation = new TraversalInformation;
            currentTraversalInformation = saveTraversalInformation(currentTraversalInformation, currentNode, currentRemainingLength, currentQuery, 'T');
            myStack.push(currentTraversalInformation);
        }
        else return false;
        break;

        default:
        print("Error! Valid characters include A, C, G and T only.\n");
        exit(-1);
    }

    // Now that something is in the stack, I can repeat the above steps for the rest of the tree:
    while(!myStack.empty() && characterMatchCount!= inputQueryLength) 
    {
        // Popping from the stack and extracting relevent information:
        currentNode = myStack.top()->address;
        currentRemainingLength = (myStack.top()->remainingLength) - 1;
        currentQueryCharacter = myStack.top()->remainingQuery[0];
        currentQuery = myStack.top()->remainingQuery;
        delete myStack.top();
        myStack.pop();

        switch(currentQueryCharacter) 
        {
            case 'A':
            if(currentNode->A != NULL)
            {
                characterMatchCount++;
                currentTraversalInformation = new TraversalInformation;
                currentTraversalInformation = saveTraversalInformation(currentTraversalInformation, currentNode, currentRemainingLength, currentQuery, 'A');
                myStack.push(currentTraversalInformation);
            }
            else return false;
            break;

            case 'C':
            if(currentNode->C != NULL) 
            {
                characterMatchCount++;
                currentTraversalInformation = new TraversalInformation;
                currentTraversalInformation = saveTraversalInformation(currentTraversalInformation, currentNode, currentRemainingLength, currentQuery, 'C');
                myStack.push(currentTraversalInformation);
            }
            else return false;
            break;

            case 'G':
            if(currentNode->G != NULL)
            {
                characterMatchCount++;
                currentTraversalInformation = new TraversalInformation;
                currentTraversalInformation = saveTraversalInformation(currentTraversalInformation, currentNode, currentRemainingLength, currentQuery, 'G');
                myStack.push(currentTraversalInformation);
            }
            else return false;
            break;

            case 'T':
            if(currentNode->T != NULL) 
            {
                characterMatchCount++;
                currentTraversalInformation = new TraversalInformation;
                currentTraversalInformation = saveTraversalInformation(currentTraversalInformation, currentNode, currentRemainingLength, currentQuery, 'T');
                myStack.push(currentTraversalInformation);
            }
            else return false;
            break;

            default:
            print("Error! Valid characters include A, C, G and T only.\n");
            exit(-1);
        }
    }
    return true;
}

// Function to search a sequence in the genome from my suffix tree:
void suffixTree::searchSuffixTreeForRead(char* sequence, int length) 
{
    int currentReadLength = length;
    // Searching my suffix tree for the current read:
    print("Searching for ");
    for(int i = 0; i < currentReadLength; i++) 
    {
        print(sequence);
    } 
    print("\n");
    bool isThere = searchSuffixTree(sequence, currentReadLength);
    if(isThere) print("Read found - Perfect match!\n");
    else        print("Read not found!\n");

    return;
}

// Function to create the readset by selecting readCount number of random reads from genome of length readLength:
char** suffixTree::select_N_RandomFragmentsFromGenome(int readCount, int readLength) 
{
    int maximumValue, randomNumber;

    // Allocating space for my reads:
    char** readset = new char*[readCount];
    for(int i = 0; i < readCount; i++) 
    {
        readset[i] = new char[readLength];
    }

    // In determining the maximum value for the genome index, I cannot select a genome index higher than (length - readLength), and
    // I add one because of the offset (starting at index 0), then subtract one because of the terminating character in the genome:
    maximumValue = genomeLength - readLength + 1 - 1; // Otherwise yeah, its meaningless to write that.

    // Creating random reads:
    for(int i = 0; i < readCount; i++) 
    {
        randomNumber = (rand() % maximumValue);
        for(int j = 0; j < readLength; j++) 
        {
            readset[i][j] = genome[randomNumber + j];
        }
    }
    return readset;
}

// Function to search my suffix tree for readCount amount of random reads made from the genome: (returns the number of reads found)
int suffixTree::searchSuffixTrieFor_N_RandomReadsFromGenome(int readCount, int readLength) 
{
    bool isThere;
    int matchCount = 0;
    // Creating readCount number of random reads from the genome:
    char** readset = select_N_RandomFragmentsFromGenome(readCount, readLength);

    // Searching the suffix tree for reads:
    auto start = high_resolution_clock::now();
    for(int currentReadIndex = 0; currentReadIndex < readCount; currentReadIndex++) 
    { 
        isThere = searchSuffixTree(readset[currentReadIndex], readLength);
        // Incrementing counter by one if the read is found:
        matchCount += (isThere) ? 1 : 0;
    }
    auto stop = high_resolution_clock::now();
    // Casting to milliseconds since I observed timings less than a second for low read count:
    auto duration = duration_cast<milliseconds>(stop - start);
    print("Search time: \n");
    print(duration.count(), " milliseconds\n");

    for(int i = 0; i < readCount; i++) 
    {
        delete[] readset[i];
    }

    return matchCount;
}