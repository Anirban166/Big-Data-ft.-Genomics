/*-------------------------
  Author: Anirban166/Ani
  Email:  ac4743@nau.edu
--------------------------*/
#include "prefixTree.h"

// Default constructor for prefixTree:
prefixTree::prefixTree() 
{
    trieNodeCount = 0;
    genomeLength = 50;
    // Save genome and reads:
    genome = new char[genomeLength];
    saveGenome();
    rowCountForReads = 50;
    readLength = 36 + 1;  // 36-character as required + one for the terminating character ($)
    readset = new char*[rowCountForReads];
    for(int i = 0; i < rowCountForReads; i++) 
    {
        readset[i] = new char[readLength];
    }
    saveNRandomFragmentsFromGenomeAsReadset();

    root = new Node;
    root->A = NULL; root->C = NULL; root->G = NULL; root->T = NULL; root->$ = NULL;
}

// Custom constructor:
prefixTree::prefixTree(std::string genomefilepath, int readCount, int readlength) 
{
    trieNodeCount = 0;
    // Save genome, reads and the genome kmer array:
    genomeFilePath = genomefilepath;
    genomeLength = getLengthOfGenome();
    genome = new char[genomeLength];
    saveGenome();
    rowCountForReads = readCount;
    readLength = readlength;
    readset = new char* [rowCountForReads];
    for(int i = 0; i < rowCountForReads; i++) 
    {
        readset[i] = new char[readLength];
    }
    saveNRandomFragmentsFromGenomeAsReadset();
    splitGenomeIntoKmer(readLength);

    root = new Node;
    root->A = NULL; root->C = NULL; root->G = NULL; root->T = NULL; root->$ = NULL;
}

// Custom constructor with per-base error taken as an additional input:
prefixTree::prefixTree(std::string genomefilepath, int readCount, int readlength, int percentError) 
{
    trieNodeCount = 0;
    genomeFilePath = genomefilepath;
    genomeLength = getLengthOfGenome();
    genome = new char[genomeLength];
    saveGenome();
    rowCountForReads = readCount;
    readLength = readlength;
    readset = new char* [rowCountForReads];
    for(int i = 0; i < rowCountForReads; i++) 
    {
        readset[i] = new char[readLength];
    }
    saveNRandomFragmentsFromGenomeAsReadset();
    splitGenomeIntoKmer(readLength);

    // Add percentError% per-base error to reads:
    introducePercentError(percentError);

    root = new Node;
    root->A = NULL; root->C = NULL; root->G = NULL; root->T = NULL; root->$ = NULL;
}

// Copy constructor:
prefixTree::prefixTree(prefixTree& object)
{
    trieNodeCount = object.trieNodeCount;
    genomeFilePath = object.genomeFilePath;   
    genomeLength = object.genomeLength;
    genome = object.genome;
    rowCountForReads = object.rowCountForReads;
    readLength = object.readLength;
    readset = object.readset;
    root = object.root;
}

// Destructor:
prefixTree::~prefixTree() 
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
    return;
}

// Function to return the number of characters in the genome:
int prefixTree::getLengthOfGenome() 
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
void prefixTree::saveGenome() 
{
    char* tempRead = new char[2000000];
    char* header = new char[100];

    genomeFile.open(genomeFilePath);
    genomeFile >> header;
    genomeFile >> tempRead;

    for(int i = 0; i < genomeLength; i++) 
    {
        genome[i] = tempRead[i];
    }

    genomeFile.close();
    delete[] tempRead;
    delete[] header;
    return;
}

// Function to save the readset as a prefix tree:
void prefixTree::createTrieReadset() 
{
    char currentCharacter;
    Node* currentNode;
    for(int n = 0; n < rowCountForReads; n++) 
    {
        currentNode = root; // Starting at root
        for(int i = 0; i < readLength; i++) 
        {
            currentCharacter = readset[n][i]; // Setting current character to insert
            switch(currentCharacter) 
            {
                #define nullAllNodes() \
                    currentNode->A = NULL; \
                    currentNode->C = NULL; \
                    currentNode->G = NULL; \
                    currentNode->T = NULL; \
                    currentNode->$ = NULL

                #define buildNodeCase(character, whichNode) \
                case character:                             \
                if(currentNode->whichNode == NULL)          \
                {                                           \
                    trieNodeCount++;                        \
                    currentNode->whichNode = new Node;      \
                    currentNode = currentNode->whichNode;   \
                    nullAllNodes();                         \
                }                                           \
                else currentNode = currentNode->whichNode;  \
                break

                buildNodeCase('A', A);
                buildNodeCase('C', C);
                buildNodeCase('G', G);
                buildNodeCase('T', T);
                buildNodeCase('$', $);
                
                default:
                print("Error! Valid characters include A, C, G, T (also $) only.\n");
                exit(-1);
            }
        }
    }
    print("Number of nodes in the trie: ", trieNodeCount, "\n");
    return;
}

// Function to save traversal information in a struct:
// (Saving the node address, the number of mismatches, the remaining query string and the length of it in a struct whose purpose is to get pushed onto the stack during the trie search)
TraversalInformation* prefixTree::saveTraversalInformation(TraversalInformation* currentTraversalData, Node* currentNode, int currentMismatchCount, int currentRemainingLength, char* currentQuery, char whichPtr) 
{
    // Determining which child address to save:
    if(whichPtr == 'A')      currentTraversalData->address = currentNode->A;
    else if(whichPtr == 'C') currentTraversalData->address = currentNode->C;
    else if(whichPtr == 'G') currentTraversalData->address = currentNode->G;
    else if(whichPtr == 'T') currentTraversalData->address = currentNode->T;
    else 
    {
        print("Error! Valid characters include A, C, G and T only.\n");
        exit(-1);
    }
    currentTraversalData->mismatchCount = currentMismatchCount;
    currentTraversalData->remainingLength = currentRemainingLength;
    currentTraversalData->remainingQuery = new char[currentRemainingLength];
    for(int i = 0; i < currentRemainingLength; i++) 
    {
        currentTraversalData->remainingQuery[i] = currentQuery[i + 1];
    }
    return currentTraversalData;
}

// Function to search the trie: (using the input genome)
int prefixTree::fuzzySearchTrie(char* inputGenomeQuery) 
{
    Node* currentNode = root;
    TraversalInformation* currentTraversalData;
    std::stack <TraversalInformation*> myStack;
    char* currentQuery = inputGenomeQuery;
    char currentQueryCharacter;
    int currentMismatchCount = 0, currentRemainingLength = 0, matchCount = 0, tolerateErrorThreshold = 1; // Tolerating up to one mismatch

    currentQueryCharacter = inputGenomeQuery[0];
    for(int i = 1; i < readLength; i++) 
    {
        currentRemainingLength++;
    }

    // Pushing the non-NULL children of the root with traversal information onto my stack:
    switch(currentQueryCharacter) 
    {
        case 'A': // If character matches with the base A, I save and push the traversal information if it ain't NULL
        if(currentNode->A != NULL) 
        {
            currentTraversalData = new TraversalInformation;
            currentTraversalData = saveTraversalInformation(currentTraversalData, currentNode, currentMismatchCount, currentRemainingLength, currentQuery, 'A');
            myStack.push(currentTraversalData);
        }
        // Otherwise, I increment the mismatch counter and save + push all the other non-NULL children to my stack:
        currentMismatchCount++;
        if(currentNode->C != NULL)
        {
            currentTraversalData = new TraversalInformation;
            currentTraversalData = saveTraversalInformation(currentTraversalData, currentNode, currentMismatchCount, currentRemainingLength, currentQuery, 'C');
            myStack.push(currentTraversalData);
        }
        if(currentNode->G != NULL)
        {
            currentTraversalData = new TraversalInformation;
            currentTraversalData = saveTraversalInformation(currentTraversalData, currentNode, currentMismatchCount, currentRemainingLength, currentQuery, 'G');
            myStack.push(currentTraversalData);
        }
        if(currentNode->T != NULL) 
        {
            currentTraversalData = new TraversalInformation;
            currentTraversalData = saveTraversalInformation(currentTraversalData, currentNode, currentMismatchCount, currentRemainingLength, currentQuery, 'T');
            myStack.push(currentTraversalData);
        }
        break;

        case 'C':
        if(currentNode->C != NULL)
        {
            currentTraversalData = new TraversalInformation;
            currentTraversalData = saveTraversalInformation(currentTraversalData, currentNode, currentMismatchCount, currentRemainingLength, currentQuery, 'C');
            myStack.push(currentTraversalData);
        }
        currentMismatchCount++;
        if(currentNode->A != NULL) 
        {
            currentTraversalData = new TraversalInformation;
            currentTraversalData = saveTraversalInformation(currentTraversalData, currentNode, currentMismatchCount, currentRemainingLength, currentQuery, 'A');
            myStack.push(currentTraversalData);
        }
        if(currentNode->G != NULL) 
        {
            currentTraversalData = new TraversalInformation;
            currentTraversalData = saveTraversalInformation(currentTraversalData, currentNode, currentMismatchCount, currentRemainingLength, currentQuery, 'G');
            myStack.push(currentTraversalData);
        }
        if(currentNode->T != NULL) 
        {
            currentTraversalData = new TraversalInformation;
            currentTraversalData = saveTraversalInformation(currentTraversalData, currentNode, currentMismatchCount,  currentRemainingLength, currentQuery, 'T');
            myStack.push(currentTraversalData);
        }
        break;

        case 'G':
        if(currentNode->G != NULL) 
        {
            currentTraversalData = new TraversalInformation;
            currentTraversalData = saveTraversalInformation(currentTraversalData, currentNode, currentMismatchCount, currentRemainingLength, currentQuery, 'G');
            myStack.push(currentTraversalData);
        }
        currentMismatchCount++;
        if(currentNode->A != NULL)
        {
            currentTraversalData = new TraversalInformation;
            currentTraversalData = saveTraversalInformation(currentTraversalData, currentNode, currentMismatchCount, currentRemainingLength, currentQuery, 'A');
            myStack.push(currentTraversalData);
        }
        if(currentNode->C != NULL) 
        {
            currentTraversalData = new TraversalInformation;
            currentTraversalData = saveTraversalInformation(currentTraversalData, currentNode, currentMismatchCount, currentRemainingLength, currentQuery, 'C');
            myStack.push(currentTraversalData);
        }
        if(currentNode->T != NULL) 
        {
            currentTraversalData = new TraversalInformation;
            currentTraversalData = saveTraversalInformation(currentTraversalData, currentNode, currentMismatchCount, currentRemainingLength, currentQuery, 'T');
            myStack.push(currentTraversalData);
        }
        break;

        case 'T':
        if(currentNode->T != NULL) 
        {
            currentTraversalData = new TraversalInformation;
            currentTraversalData = saveTraversalInformation(currentTraversalData, currentNode, currentMismatchCount, currentRemainingLength, currentQuery, 'T');
            myStack.push(currentTraversalData);
        }
        currentMismatchCount++;
        if(currentNode->A != NULL) 
        {
            currentTraversalData = new TraversalInformation;
            currentTraversalData = saveTraversalInformation(currentTraversalData, currentNode, currentMismatchCount, currentRemainingLength, currentQuery, 'A');
            myStack.push(currentTraversalData);
        }
        if(currentNode->C != NULL) 
        {
            currentTraversalData = new TraversalInformation;
            currentTraversalData = saveTraversalInformation(currentTraversalData, currentNode, currentMismatchCount, currentRemainingLength, currentQuery, 'C');
            myStack.push(currentTraversalData);
        }
        if(currentNode->G != NULL) 
        {
            currentTraversalData = new TraversalInformation;
            currentTraversalData = saveTraversalInformation(currentTraversalData, currentNode, currentMismatchCount, currentRemainingLength, currentQuery, 'G');
            myStack.push(currentTraversalData);
        }
        break;

        case '$':  // Leaf node case (match found)
        matchCount++;
        break;

    default:
        print("Error! Valid characters include A, C, G and T only.\n");
        exit(-1);
    }

    // Now that something is in the stack, I can repeat the above steps for the rest of the trie:
    while(!myStack.empty()) 
    {
        // Popping from the stack and extracting relevent information:
        currentNode = myStack.top()->address;
        currentMismatchCount = myStack.top()->mismatchCount;
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
                currentTraversalData = new TraversalInformation;
                currentTraversalData = saveTraversalInformation(currentTraversalData, currentNode, currentMismatchCount, currentRemainingLength, currentQuery, 'A');
                myStack.push(currentTraversalData);
            }
            currentMismatchCount++;
            // Check if the number of mismatches exceed the threshold:
            if(currentMismatchCount > tolerateErrorThreshold) 
            {
                break;
            }
            // If not, then continue with the same deal, i.e. save and push all the other non-NULL children to my stack:
            else 
            {
                if(currentNode->C != NULL) 
                {
                    currentTraversalData = new TraversalInformation;
                    currentTraversalData = saveTraversalInformation(currentTraversalData, currentNode, currentMismatchCount, currentRemainingLength, currentQuery, 'C');
                    myStack.push(currentTraversalData);
                }
                if(currentNode->G != NULL) 
                {
                    currentTraversalData = new TraversalInformation;
                    currentTraversalData = saveTraversalInformation(currentTraversalData, currentNode, currentMismatchCount, currentRemainingLength, currentQuery, 'G');
                    myStack.push(currentTraversalData);
                }
                if(currentNode->T != NULL) 
                {
                    currentTraversalData = new TraversalInformation;
                    currentTraversalData = saveTraversalInformation(currentTraversalData, currentNode, currentMismatchCount, currentRemainingLength, currentQuery, 'T');
                    myStack.push(currentTraversalData);
                }
                break;
            }
            case 'C':
            if(currentNode->C != NULL) 
            {
                currentTraversalData = new TraversalInformation;
                currentTraversalData = saveTraversalInformation(currentTraversalData, currentNode, currentMismatchCount, currentRemainingLength, currentQuery, 'C');
                myStack.push(currentTraversalData);
            }
            currentMismatchCount++;
            if(currentMismatchCount > tolerateErrorThreshold) 
            {
                break;
            }
            else 
            {
                if(currentNode->A != NULL) 
                {
                    currentTraversalData = new TraversalInformation;
                    currentTraversalData = saveTraversalInformation(currentTraversalData, currentNode, currentMismatchCount, currentRemainingLength, currentQuery, 'A');
                    myStack.push(currentTraversalData);
                }
                if(currentNode->G != NULL) 
                {
                    currentTraversalData = new TraversalInformation;
                    currentTraversalData = saveTraversalInformation(currentTraversalData, currentNode, currentMismatchCount, currentRemainingLength, currentQuery, 'G');
                    myStack.push(currentTraversalData);
                }
                if(currentNode->T != NULL) 
                {
                    currentTraversalData = new TraversalInformation;
                    currentTraversalData = saveTraversalInformation(currentTraversalData, currentNode, currentMismatchCount, currentRemainingLength, currentQuery, 'T');
                    myStack.push(currentTraversalData);
                }
                break;
            }
            case 'G':
            if(currentNode->G != NULL) 
            {
                currentTraversalData = new TraversalInformation;
                currentTraversalData = saveTraversalInformation(currentTraversalData, currentNode, currentMismatchCount, currentRemainingLength, currentQuery, 'G');
                myStack.push(currentTraversalData);
            }
            currentMismatchCount++;
            if(currentMismatchCount > tolerateErrorThreshold) 
            {
                break;
            }
            else 
            {
                if(currentNode->A != NULL) 
                {
                    currentTraversalData = new TraversalInformation;
                    currentTraversalData = saveTraversalInformation(currentTraversalData, currentNode, currentMismatchCount, currentRemainingLength, currentQuery, 'A');
                    myStack.push(currentTraversalData);
                }
                if(currentNode->C != NULL) 
                {
                    currentTraversalData = new TraversalInformation;
                    currentTraversalData = saveTraversalInformation(currentTraversalData, currentNode, currentMismatchCount, currentRemainingLength, currentQuery, 'C');
                    myStack.push(currentTraversalData);
                }
                if(currentNode->T != NULL) 
                {
                    currentTraversalData = new TraversalInformation;
                    currentTraversalData = saveTraversalInformation(currentTraversalData, currentNode, currentMismatchCount, currentRemainingLength, currentQuery, 'T');
                    myStack.push(currentTraversalData);
                }
                break;
            }
            case 'T':
            if(currentNode->T != NULL) 
            {
                currentTraversalData = new TraversalInformation;
                currentTraversalData = saveTraversalInformation(currentTraversalData, currentNode, currentMismatchCount, currentRemainingLength, currentQuery, 'T');
                myStack.push(currentTraversalData);
            }
            currentMismatchCount++;
            if(currentMismatchCount > tolerateErrorThreshold) 
            {
                break;
            }
            else 
            {
                if(currentNode->A != NULL) 
                {
                    currentTraversalData = new TraversalInformation;
                    currentTraversalData = saveTraversalInformation(currentTraversalData, currentNode, currentMismatchCount, currentRemainingLength, currentQuery, 'A');
                    myStack.push(currentTraversalData);
                }
                if(currentNode->C != NULL) 
                {
                    currentTraversalData = new TraversalInformation;
                    currentTraversalData = saveTraversalInformation(currentTraversalData, currentNode, currentMismatchCount, currentRemainingLength, currentQuery, 'C');
                    myStack.push(currentTraversalData);
                }
                if(currentNode->G != NULL) 
                {
                    currentTraversalData = new TraversalInformation;
                    currentTraversalData = saveTraversalInformation(currentTraversalData, currentNode, currentMismatchCount, currentRemainingLength, currentQuery, 'G');
                    myStack.push(currentTraversalData);
                }
                break;
            }

            case '$':
            matchCount++;
            break;

            default:
            print("Error! Valid characters include A, C, G and T only.\n");
            exit(-1);
        }
    }
    return matchCount;
}

// Function to create the readset by selecting readCount amount of random fragments from the genome (of length readLength):
void prefixTree::saveNRandomFragmentsFromGenomeAsReadset() 
{
    int randomNumber, maximumBound = genomeLength - readLength + 1 + 1;  // Cannot select a genome index higher than (length - readLength)
    srand(time(0));
    for(int i = 0; i < rowCountForReads; i++) 
    {
        randomNumber = (rand() % maximumBound);
        for(int j = 0; j < readLength - 1; j++) 
        {
            readset[i][j] = genome[randomNumber + j];
        }
        readset[i][readLength - 1] = '$';
    }
    return;
}

// Function to update the readset by introducing x% per-base error rate:
void prefixTree::introducePercentError(int percentError) 
{
    int randomNumberOne, randomNumberTwo;
    char original;
    srand(time(0));
    for(int i = 0; i < rowCountForReads; i++) 
    {
        for(int j = 0; j < readLength - 1; j++) 
        {
            // Checking if the number obtained from 0 to 99 lies within the safe zone for x% error, i.e. within 0 to (99 - x): (if it doesn't, then it hits the x% error rate)
            randomNumberOne = rand() % 100 + 1;
            if(randomNumberOne > 100 - percentError) 
            {
                original = readset[i][j]; // Saving the original character
                while(readset[i][j] == original) 
                {   
                    randomNumberTwo = (rand() % 4);
                    if(randomNumberTwo == 0)      readset[i][j] = 'A';
                    else if(randomNumberTwo == 1) readset[i][j] = 'C';
                    else if(randomNumberTwo == 2) readset[i][j] = 'G';
                    else if(randomNumberTwo == 3) readset[i][j] = 'T';
                    else 
                    {
                        print("Invalid input.\nExiting!\n");
                        exit(-1);
                    }
                }
            }
        }
    }
    return;
}

// Function to split a fragment into kmers, save it as an array and then return it:
void prefixTree::splitGenomeIntoKmer(int kmer) 
{
    int offset = 0;
    rowCountForGenomeKmerArray = genomeLength - kmer + 1 + 1; 

    // Initializing the kmer array, now that I know the size of it:
    genomeKmerArray = new char*[rowCountForGenomeKmerArray];
    for(int i = 0; i < rowCountForGenomeKmerArray; i++) 
    {
        genomeKmerArray[i] = new char[kmer];
    }

    // Splitting into kmer fragments and save as a 2D array: 
    for(int r = 0; r < rowCountForGenomeKmerArray; r++) 
    {
        for(int i = 0; i < kmer - 1; i++) 
        {
            genomeKmerArray[r][i] = genome[i + offset];
        }
        genomeKmerArray[r][kmer - 1] = '$';
        offset++;
    }
    return;
}

// Function to search the trie for all genome kmers:
int prefixTree::searchTrieForGenomeKmers() 
{
    int totalMatchCount = 0, currentKmerMatchCount, genomeKmerMatchCount = 0;
    for(int i = 0; i < rowCountForGenomeKmerArray; i++) 
    {
        currentKmerMatchCount = fuzzySearchTrie(genomeKmerArray[i]);
        totalMatchCount += currentKmerMatchCount;
        genomeKmerMatchCount += (currentKmerMatchCount > 0) ? 1 : 0;
    }
    print("Number of genome kmers with at least one match: ", genomeKmerMatchCount, "\n");
    return totalMatchCount;
}