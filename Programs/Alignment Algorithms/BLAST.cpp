/*-------------------------
  Author: Anirban166/Ani
  Email:  ac4743@nau.edu
--------------------------*/
#include "BLAST.h"

// BLAST class constructor:
BLAST::BLAST(int scoreForMatch, int penalityForMismatch, int penaltyForGap) 
{
    readLength = 50;
    matchScore = scoreForMatch;
    gapScore = penaltyForGap, 
    mismatchScore = penalityForMismatch;
}

// textAlignment class constructor:
textAlignment::textAlignment(int len) 
{
	length = len;
    code = new char[length];
	sequenceOneText = new char[length];
	sequenceTwoText = new char[length];
}

// Function to save the alignment:
// (Note that during traceback, the alignment was saved backwards, so this function reverses the order!)
void textAlignment::saveAlignment(int alignmentLength, char* sequenceOne, char* sequenceTwo, char* alignmentCode) 
{
    // Save the first sequence:
    for(int i = 0; i < alignmentLength; i++) 
    {
        sequenceOneText[i] = sequenceOne[alignmentLength - i - 1];
    }
    // Then the second one:
    for(int i = 0; i < alignmentLength; i++) 
    {
        sequenceTwoText[i] = sequenceTwo[alignmentLength - i - 1];
    }
    // Save the code:
    for(int i = 0 ;i < alignmentLength; i++) 
    {
        if(alignmentCode[alignmentLength - i - 1] == 'L' || alignmentCode[alignmentLength - i - 1] == 'U') 
        {
            code[i] = ' ';
        }
        else 
        {
            code[i] = (sequenceOneText[i] != sequenceTwoText[i]) ? 'x' : '|';
        }
    }
    return;
}

// Function to print the text alignment:
void textAlignment::printAlignment(int alignmentLength) 
{
    for(int p = 0; p < alignmentLength; p++) 
    {
        print(sequenceOneText[p]);
    }
    print("\n");
    for(int p = 0; p < alignmentLength; p++) 
    {
        print(code[p]);
    }
    print("\n");
    for(int p = 0; p < alignmentLength; p++) 
    {
        print(sequenceTwoText[p]);
    }
    print("\n");
    return;
}

// Function to perform the traceback computation in accordance with the final phase of the BLAST algorithm:
int BLAST::performTraceback(int** D, char** traceback, int* indexOne, int* indexTwo, char* tempSequenceOne, char* tempSequenceTwo, char* tempCode, char* sequenceOne, char* sequenceTwo) 
{
    int alignmentIndex;
    alignmentIndex = 0; 
    while(D[*indexTwo][*indexOne] != 0) 
    {
        if(traceback[*indexTwo][*indexOne] == 'D') 
        {   // Save diagonal (D), and nucleotide sequences one and two:
            tempCode[alignmentIndex] = 'D'; // Save diagonal (D)
            tempSequenceOne[alignmentIndex] = sequenceOne[*indexOne - 1];
            tempSequenceTwo[alignmentIndex] = sequenceTwo[*indexTwo - 1];
            *indexOne = *indexOne - 1; // diagonal move = j - 1 (column)
            *indexTwo -= 1;           // diagonal move = i - 1 (row)
        }
        else if(traceback[*indexTwo][*indexOne] == 'L') 
        {   // L (left), save nucelotide sequence one and gap for sequence two
            tempCode[alignmentIndex] = 'L';
            tempSequenceOne[alignmentIndex] = sequenceOne[*indexOne - 1];
            tempSequenceTwo[alignmentIndex] = '_';
            *indexOne -= 1;     // left move = j-1 (col)
        }
        else 
        {   // U (up), reverse of L (i.e., save gap for sequence one and nucelotide sequence two)
            tempCode[alignmentIndex] = 'U'; 
            tempSequenceOne[alignmentIndex] = '_';
            tempSequenceTwo[alignmentIndex] = sequenceTwo[*indexTwo - 1];
            *indexTwo -= 1;
        }
        alignmentIndex++;
    }
    return alignmentIndex;
}

// Function to get the maximum score in the diagonal array:
int BLAST::findMaxInArray(int** D, int* indexOne, int* indexTwo, int sequenceOneLength, int sequenceTwoLength)
{
    int max = D[*indexTwo][*indexOne]; // Storing the first array element as the maximum initially.
    // Iterating through the array whilst skipping the first row and column:
    for(int i = 1; i < sequenceTwoLength + 1; i++)  
    {
        for(int j = 1; j < sequenceOneLength + 1; j++) 
        {
            if(D[i][j] > D[*indexTwo][*indexOne]) 
            {
                *indexTwo = i; 
                *indexOne = j;
                max = D[*indexTwo][*indexOne]; // Update maximum score value.
            }
        }
    }
    return max;
}

// Function to perform the array computation for the alignment algorithm: (creates both the traceback and diagonal arrays)
void BLAST::computeArray(int** D, char** traceback, int matchScore, int gapScore, int mismatchScore, char* sequenceOne, char* sequenceTwo, int sequenceOneLength, int sequenceTwoLength) 
{
    int scoreOne, scoreTwo, scoreThree;
    // Iterating through the array whilst skipping the first row and column:
    for(int j = 1; j < sequenceOneLength + 1; j++) 
    {
        for(int i = 1; i < sequenceTwoLength + 1; i++) 
        {
            // Computing the scores:
            if(sequenceOne[j - 1] == sequenceTwo[i - 1]) // Case I: nucelotides match
            {
                scoreOne = D[i - 1][j - 1] + matchScore;
            }
            else // Mismatch!
            {                                         
                scoreOne = D[i - 1][j - 1] + mismatchScore;
            }
            scoreTwo = D[i - 1][j] + gapScore;
            scoreThree = D[i][j - 1] + gapScore;

            // Computing the diagonal and traceback arrays:
            if(scoreOne >= scoreTwo) 
            {
                if(scoreOne >= scoreThree) 
                { 
                    D[i][j] = scoreOne;     
                    traceback[i][j] = 'D';
                }
                else 
                {         
                    D[i][j] = scoreThree;     
                    traceback[i][j] = 'L';
                }
            }
            else 
            {
                if(scoreTwo >= scoreThree) 
                { 
                    D[i][j] = scoreTwo;  
                    traceback[i][j] = 'U';
                }
                else 
                {    
                    D[i][j] = scoreThree;    
                    traceback[i][j] = 'L';
                }
            } 
        }
    }
    return;
}

// Function to perform all of the steps in the BLAST algorithm:
textAlignment* BLAST::BLASTAlignment(char* sequenceOne, char* sequenceTwo, int sequenceOneLength, int sequenceTwoLength,  int matchScore, int gapScore, int mismatchScore, int* max)
{
    int** D;
    char** traceback;
    int* indexOne;      // j (col)
    int* indexTwo;     // i (row)
    indexOne = new int;
    indexTwo = new int;
    *indexOne = 1;
    *indexTwo = 1;
    int alignmentLength;
    // Temporary arrays for saving the text alignment:
    char* tempArrayOne = new char[100];
    char* tempArrayTwo = new char[100];
    char* tempArrayThree = new char[100];
    // Allocate memory for, and initialize the diagonal and traceback arrays:
    D = new int* [sequenceTwoLength + 1];
    traceback = new char* [sequenceTwoLength + 1];
    for(int i = 0; i < sequenceTwoLength + 1; i++) 
    {
        D[i] = new int[sequenceOneLength + 1];
        traceback[i] = new char[sequenceOneLength + 1];
    }
    for(int i = 0; i < sequenceTwoLength + 1; i++) 
    {
        for(int j = 0; j < sequenceOneLength + 1; j++) 
        {
            D[i][j] = 0;
            traceback[i][j] = '0';
        }
    }

    // BLAST time:
    computeArray(D, traceback, matchScore, gapScore, mismatchScore, sequenceOne, sequenceTwo, sequenceOneLength, sequenceTwoLength);    
    *max = findMaxInArray(D, indexOne, indexTwo, sequenceOneLength, sequenceTwoLength);
    alignmentLength = performTraceback(D, traceback, indexOne, indexTwo, tempArrayOne, tempArrayTwo, tempArrayThree, sequenceOne, sequenceTwo);

    // Save the text alignment:
    textAlignment* output = new textAlignment(alignmentLength);
    output->saveAlignment(alignmentLength, tempArrayOne, tempArrayTwo, tempArrayThree);

    // Deallocate/free memory for used variables:
    for(int i = 0; i < sequenceTwoLength + 1; i++) 
    {
        delete[] D[i];
        delete[] traceback[i];
    }
    delete[] D;
    delete[] traceback;
    delete indexOne; delete indexTwo;
    delete[] tempArrayOne; delete[] tempArrayTwo; delete[] tempArrayThree;

    return output;
}

// Function to return the number of characters in the genome file:
int BLAST::getLengthOfGenome(std::string genomeFilePath) 
{
    std::ifstream genomeFile;
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

// Function to save the genome file (.fa one) as a character array:
void BLAST::saveGenome(std::string genomeFilePath, char* genome, int genomeLength) 
{
    std::ifstream genomeFile;
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

// Function to return the number of rows (or reads) in the readset file:
int BLAST::countNumberOfReads(std::string readFilePath) 
{
    std::ifstream readFile;
    int count = 0;
    char* temp = new char[1000];
    readFile.open(readFilePath);
    while(!readFile.eof()) 
    {
        readFile >> temp;
        readFile >> temp;
        count++;
    }
    readFile.close();
    // varDebug(count);
    delete[] temp;
    return count;
}

// Function to save the reads into a character array:
void BLAST::saveReads(int readCount, int readLength, char** readset, std::string readFilePath) 
{
    std::ifstream readFile;
    char* temp = new char[1000];
    readFile.open(readFilePath);
    
    for(int i = 0; i < readCount; i++) 
    {
        readFile >> temp;
        readFile >> temp;
        for(int j = 0; j < readLength; j++) 
        {
            readset[i][j] = temp[j];
        }
    }
    
    readFile.close();
    delete[] temp;
    return;
}

// Function to create specified number (based on input readCount) of random reads (with the four bases A/C/G/T) of specified (input) length:
void BLAST::createRandomReads(char** readset, int readCount, int readLength) 
{
    srand(time(0));
    int randomNumber;
    for(int i = 0; i < readCount; i++) 
    {
        for(int j = 0; j < readLength; j++) 
        {
            randomNumber = (rand() % 4);
            if(randomNumber == 0)      readset[i][j] = 'A';
            else if(randomNumber == 1) readset[i][j] = 'C';
            else if(randomNumber == 2) readset[i][j] = 'G';
            else if(randomNumber == 3) readset[i][j] = 'T';
            else 
            {
                print("Error, invalid input.\nExiting!\n");
                exit(-1);
            }
        }
    }
    return;
}

// Function to split fragment into kmers, and to return that saved array:
char** BLAST::splitFragmentIntoKmer(char* input, int inputLength, int kmer, int* rowCount) 
{
    char** output;
    int j = 0;
    *rowCount = inputLength - kmer + 1;

    // Initializing the kmer array, now that I know the size of it:
    output = new char* [*rowCount];
    for(int i = 0; i < *rowCount; i++) 
    {
        output[i] = new char[kmer];
    }

    // Splitting into kmer fragments and save as a 2D array: 
    for(int r = 0; r < *rowCount; r++) 
    {
        for(int i = 0; i < kmer; i++) 
        {
            output[r][i] = input[i + j];
        }
        j++;
    }
    return output;
}

// Function to create and initialize the hash table: (with collision-based chaining implemented with linked lists)
Node** BLAST::createHashTable(int hashTableSize) 
{
    Node** T;
    // Allocating memory for the array of linked lists:
    T = new Node * [hashTableSize];
    for(int i = 0; i < hashTableSize; i++) 
    {
        T[i] = new Node;
    }

    // Initializing em':
    T[0]->data = 1;
    T[0]->next = NULL;
    for(unsigned int i = 1; i < hashTableSize; i++) 
    {
        T[i]->data = 0;
        T[i]->next = NULL;
    }
    return T;
}

// Function to insert a character array into my hash table:
void BLAST::insertFragmentInTable(char* fragment, int kmer, int hashTableSize, Node** T, int* genomeIndex) 
{
    int radixNumber, address;
    radixNumber = convertFragmentToRadixNumber(fragment, kmer);
    address = divisionMethodHash(radixNumber, hashTableSize);
    insertInTable(address, radixNumber, T, hashTableSize, genomeIndex);
    return;
}

// Function to convert the sequence fragment (character array) into a radix number:
int BLAST::convertFragmentToRadixNumber(char* inputRead, int kmer) 
{
    int* readInteger = new int[kmer];
    char* readCharacter = new char[kmer];    
    int radixNumber; // Radix number for a single fragment

    for(int j = 0; j < kmer; j++) 
    {
        readCharacter[j] = inputRead[j];
    }
    readInteger = characterToIntegerConversion(readCharacter, readInteger, kmer);
    radixNumber = radixBaseConversion(readInteger, kmer);

    delete[] readInteger;
    delete[] readCharacter;
    return radixNumber;
}

// Function to convert 16-char array with char (A, C, G, T) to 16-int array (0, 1, 2, 3):
int* BLAST::characterToIntegerConversion(char* readCharacter, int* readInteger, int kmer) 
{
    for(int i = 0; i < kmer; i++) 
    {
        if(readCharacter[i] == 'A')      readInteger[i] = 0;
        else if(readCharacter[i] == 'C') readInteger[i] = 1;
        else if(readCharacter[i] == 'G') readInteger[i] = 2;
        else if(readCharacter[i] == 'T') readInteger[i] = 3;
        else 
        {
            print("Invalid character, must be one of A,C,G,T. Exiting!\n");
            exit(-1);
        }
    }
    return readInteger;
}

// Function to convert the integer sequence array into a radix number:
int BLAST::radixBaseConversion(int* readInteger, int kmer) 
{
    int radixNumber = 0, base = 4, exponent;
    // Radix number = (value[0])(base^sequenceLength - 1) + (value[1])(base^sequenceLength - 2) + ...
    for(int i = 0; i < kmer; i++) 
    {
        exponent = kmer - i - 1;
        radixNumber += readInteger[i] * (unsigned int)(pow(base, exponent));
    }
    return radixNumber;
}

// Function to convert radix numbers into addresses using the division method:
int BLAST::divisionMethodHash(int radixNumber, int hashTableSize) 
{
    int tableAddress;
    tableAddress = radixNumber % hashTableSize;
    return tableAddress;
}

// Function to insert radix number at specified (input) address:
void BLAST::insertInTable(int address, int radixNumber, Node** T, int hashTableSize, int* genomeIndex) 
{
    bool isThere;                 // Boolean indicating whether the radix number is in table or not
    Node* current = T[address];  // Temporary pointer to be made to point to the node of interest
    // Checking if the radix number is already in the table:
    isThere = searchTableForRadix(T[address], radixNumber, genomeIndex);
    if(isThere) return;       // Duplicates
    else 
    {
        if ((current->data) % hashTableSize != radixNumber % hashTableSize) 
        {
            current->data = radixNumber;
            current->position = *genomeIndex;
            return;
        }
        else if ((current->data) % hashTableSize == radixNumber % hashTableSize) 
        {
            T[address] = new Node;
            T[address]->data = radixNumber;
            T[address]->position = *genomeIndex;
            T[address]->next = current;
            return;
        }
        else 
        {
            print("Error in table insertion, exiting the program!\n");
            exit(-1);
        }
    }
}

// Function to search my hash table for a radix number:
bool BLAST::searchTableForRadix(Node* tableAddress, int radixNumber, int* genomeIndex) 
{
    Node* current;
    current = tableAddress;
    if(current->data == radixNumber)
    {
        *genomeIndex = current->position;
        return true;
    }
    else 
    {
        if(current->next == NULL) 
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
                    *genomeIndex = current->position;
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

// Function to search the hash table for a given query within the BLAST class:
bool BLAST::searchTableForFragment(char* fragment, int kmer, int hashTableSize, Node** T, int* genomeIndex) 
{
    int radixNumber = convertFragmentToRadixNumber(fragment, kmer);
    int address = divisionMethodHash(radixNumber, hashTableSize);
    bool isThere = searchTableForRadix(T[address], radixNumber, genomeIndex);
    return isThere;
}

// Function to return the length of the sub genome: (sets lower and upper bounds prior to that)
// (Any expansion below the minimum genome index (0) above the maximum (genome length) is capped off at zero and max respectively)
int BLAST::getSubGenomeLength(int* lowerBound, int* upperBound, int* genomeIndex, int currentReadKmer, int kmer, int currentReverseIndex, int genomeLength) 
{
    int subGenomeLength;
    *lowerBound = *genomeIndex - currentReadKmer;
    if(*lowerBound < 0) 
    {
        *lowerBound = 0;
    }
    *upperBound = *genomeIndex + kmer + currentReverseIndex;
    if(*upperBound > genomeLength) 
    {
        *upperBound = genomeLength;
    }
    subGenomeLength = *upperBound - *lowerBound;
    return subGenomeLength;
}

// Function to expand the seed and save the sub-genome:
void BLAST::expandSeed(int* lowerBound, int* upperBound, char* genome, char* subGenome) 
{
    int index = 0;
    for(int i = *lowerBound; i < *upperBound; i++) 
    {
        subGenome[index] = genome[i];
        index++;
    }
    return;
}

// Function to search and expand a seed:
char* BLAST::searchAndExpandSeed(int readNumber, char** readset, int readLength, int kmer, int hashTableSize, int* genomeIndex, Node** T, int genomeLength, char* genome, int* subGenomeLength, bool* seedFound, char* subGenome, bool printFlag) 
{
    int* kmerReadRowCount = new int;
    int* lowerBound = new int;
    int* upperBound = new int;
    bool isThere;
    char** readKmer = splitFragmentIntoKmer(readset[readNumber], readLength, kmer, kmerReadRowCount);
    int notFoundCount = 0, currentReadKmer = 0, currentReverseIndex = *kmerReadRowCount - 1;
    *subGenomeLength = 1;
    while(currentReadKmer < *kmerReadRowCount) 
    {
        isThere = searchTableForFragment(readKmer[currentReadKmer], kmer, hashTableSize, T, genomeIndex);
        if(isThere) 
        {
            if(printFlag) 
            {
                print("First read found:\n");
            }
            *subGenomeLength = getSubGenomeLength(lowerBound, upperBound, genomeIndex, currentReadKmer, kmer, currentReverseIndex, genomeLength);
            subGenome = new char[*subGenomeLength];
            expandSeed(lowerBound, upperBound, genome, subGenome);
            *seedFound = true;
            break; // Stopping at the first seed that has been found.
        }
        else 
        {
            notFoundCount++;
            currentReadKmer++;
            currentReverseIndex--;
        }
        if(notFoundCount == *kmerReadRowCount) 
        {
            if(printFlag) 
            {
                print("No matches found, sorry!\n");
            }
            subGenome = new char;
            *subGenome = '0'; // Default
            *seedFound = false;
        }
    }
    for(int i = 0; i < *kmerReadRowCount; i++) 
    {
        delete[] readKmer[i];
    }
    delete[] readKmer;
    delete lowerBound;
    delete upperBound;
    delete kmerReadRowCount;
    return subGenome;
}

// Function to perform seed-based BLAST for specified (readCount) number of random reads:
void BLAST::performBLASTForNRandomReads(int readCount, int readLength, int kmer, int hashTableSize, Node** T, char* genome, int genomeLength, int* genomeIndex, int matchScore, int gapScore, int mismatchScore, int* max, bool printFlag) 
{
    int* subGenomeLength = new int;
    char* subGenome = nullptr;
    char** readset;
    bool* seedFound = new bool;
    print("For ", readCount, " reads:\n");
    readset = new char* [readCount];
    for(int i = 0; i < readCount; i++) 
    {
        readset[i] = new char[readLength];
    }
    createRandomReads(readset, readCount, readLength);

    // For each read, I split it into kmers, search for each kmer seed, and perform my alignment algorithm if applicable: (also extracting the time taken to complete the process)
    auto start = high_resolution_clock::now();
    for(int readNum = 0; readNum < readCount; readNum++) 
    {
        subGenome = searchAndExpandSeed(readNum, readset, readLength, kmer, hashTableSize, genomeIndex, T, genomeLength, genome, subGenomeLength, seedFound, subGenome, printFlag);
        if(*seedFound) 
        {
            textAlignment* obj = BLASTAlignment(subGenome, readset[readNum], *subGenomeLength, readLength, matchScore, gapScore, mismatchScore, max);
        }
        delete[] subGenome;
    }
    auto stop = high_resolution_clock::now();                     
    auto duration = duration_cast<seconds>(stop - start);
    print("Alignment time: \n");
    print(duration.count(), " seconds\n");

    for(int i = 0; i < readCount; i++) 
    {
        delete[] readset[i];
    }
    delete[] readset;
    delete seedFound;    
    delete subGenomeLength;
    return;
}

// Function to create the readset by selecting specified number (based on the input readCount) of random fragments from the genome: (randomly selects the specified amount of reads from genome of length readLength)
char** BLAST::selectNRandomFragmentsFromGenome(int readCount, int readLength, char* genome, int genomeLength) 
{
    char** readset;
    int maxBound;
    int randNum;
    readset = new char* [readCount];
    for(int i = 0; i < readCount; i++) 
    {
        readset[i] = new char[readLength];
    }
    maxBound = genomeLength - readLength + 1;  // Since I can't select a genome index higher than (length - readLength)
    srand(time(0));
    for(int i = 0; i < readCount; i++) 
    {
        randNum = (rand() % maxBound);
        for(int j = 0; j < readLength; j++) 
        {
            readset[i][j] = genome[randNum + j];
        }
    }
    return readset;
}

// Function to update the readset by introducing 5% error rate of base change:
void BLAST::introducePercentError(int readCount, int readLength, char** readset, int percentError) 
{
    int randomNumberOne, randomNumberTwo;
    char original;
    srand(time(0));
    for(int i = 0; i < readCount; i++) 
    {
        for(int j = 0; j < readLength; j++) 
        {    // x% error rate introduction can be done in multiple ways as I can think of. The one I was used for the previous problem (check hash table constructs) was with the range, i.e. it would be (95 + 1 - 0) + 0 here
            // For this one, I will check if the number obtained from 0 to 99 lies within the safe zone for (100 - 95) = 5% error, i.e. within 0 to 94. If it doesn't, then it hits the 5% error rate.
            randomNumberOne = rand() % 100;
            if(randomNumberOne > 100 - percentError) 
            {
                original = readset[i][j]; // Saving original character
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

// Function to delete the hash table:
void BLAST::deleteHashTable(int hashTableSize, Node** T) 
{
    Node* tempHead;
    Node* current;
    // Looping through the array and deleting each pointer to the linked list:
    for(int i = 0; i < hashTableSize; i++) 
    {
        tempHead = T[i];
        while(tempHead->next != NULL) 
        {
            current = tempHead->next;
            delete tempHead;
            tempHead = current;
        }
        if(tempHead->next == NULL) 
        {
            delete tempHead;
        }
    }
    return;
}