/*-------------------------
  Author: Anirban
  Email:  ac4743@nau.edu
--------------------------*/
#include "matrix.h"
// Compile and run locally using: g++ -std=c++17 <filename>.cpp -o <executablename> && ./<executablename> <problemname> <filepath>
// Example: g++ -std=c++17 arrayVersion.cpp -o arrayVersion && ./arrayVersion <problemID> sampleDataset.fa
// Alternatively, use the make file: make && ./arrayVersion <problemID> sampleDataset.fa
// Slurm-based compute cluster (Monsoon for instance): make && sbatch <runnerscript>.sh
// Use -O3 flag during compilation for better compiler optimization, and/or multiple threads for obtaining results faster.

// Count the total number of rows present in our entire data set: (outputs the total row count in that file)
int FASTA_readset::countNumberOfReads()
{
    int count = -1;
    char* tempData = new char[1000];
    dataFile.open(filePath);
    while (!dataFile.eof()) 
	{   // One for the header, one for the read: (i.e. the read count would be half the number of lines in the file)
        dataFile >> tempData;
        dataFile >> tempData;
        count++;
    }
    dataFile.close();
    return count;
}

// Function to save a row into the countArray: (inputs: pointer to header, length of header, row index of countArray)
void FASTA_readset::saveCountArray(char* header, int headerLength, int index)
{
    int a, readNum, start, end, m = 1, count = 0;
    char* buffer = new char[100];                
    while (m < headerLength)   
    {                          
        start = m + 1; // Offset header index by 1 (to skip '_'/'R') 
        end = start;
        while(header[end] != '_' && end < headerLength)
        {
            end++;
        }
        for (a = start; a < end; a++) // From start of the subsection up till the end, save into buffer
        {
            buffer[a - start] = header[a];
        }
        buffer[a] = '\0';              
        readNum = atoi(buffer);       
        countArray[index][count]= readNum; 
        for (int p = 0; p < (a - start) + 1; p++)
        {
            buffer[p] = '\0';
        }
        m = end;
        count++;
    }
}
   
// Function to save a row into the readArray: (inputs: pointer to temp read, row index of readArray)
void FASTA_readset::saveReadArray(char* read, int index)
{
    int readLength = 50; // Since each read is exactly 50 characters (or rather nucleotides) long
    for (int j = 0; j < readLength; j++) 
    {
        readArray[index][j] = read[j];
    }
}

// Function to save data from the file (saves header as countArray and reads as readArray)
void FASTA_readset::saveData()
{
    char* header = new char[1000]; 
    char* read = new char[1000];
    int index;                                            
    int headerLength;                                      
    dataFile.open(filePath);                             
    for (int i = 0; i < row; i++)
    {                                                      
        dataFile >> header;                            
        dataFile >> read;                            
        headerLength = getHeaderLength(header); 
        index = i;                                    
        saveCountArray(header, headerLength, index); 
        saveReadArray(read, index);                  
    }
    dataFile.close();
    print("\nStatus update: The data has been saved successfully.\n");
}

// Function to create arrays for holding datasets one and two, and to extract the specific fragments into them from the entire read array:
void FASTA_readset::createDataSetOneAndTwo()
{
    // readArray dimensions for dataset one/two = [row < 5000][50] = [uniqueSequenceArray[0]/[1]][50]
    int dataSetOneLength = uniqueSequenceArray[0], dataSetTwoLength = uniqueSequenceArray[1]; 
    dataSetOneReadArray = new char*[dataSetOneLength];
    dataSetTwoReadArray = new char*[dataSetTwoLength];
    // Initialize the char array of arrays that hold data sets one and two:
    for(int i = 0; i < dataSetOneLength; i++)
    {
        dataSetOneReadArray[i] = new char[readCol]; // [dataset one size]x[50] char array
    }
    for(int i = 0; i < dataSetTwoLength; i++)
    {
        dataSetTwoReadArray[i] = new char[readCol]; // [dataset two size]x[50] char array
    }   
    // I loop through all the rows with it being the first index, then select ones with second column of countArray to be less than
    // or not equal to 0, which would pertain to data set one's fragments. Then, I extract that sequence from the readArray and increment
    // the value of dataset one's index (has to be separate - initial mistake) within the condition's scope to get my desired array:
    for(int i = 0, j = 0; i < row, j < dataSetOneLength; i++)
    {
        if(countArray[i][1] > 0)
        {
            dataSetOneReadArray[j] = readArray[i];
            j++;
        }
    }
    // Likewise, I apply the same logic to get dataset two:
    for(int i = 0, j = 0; i < row, j < dataSetTwoLength; i++)
    {
        if(countArray[i][2] != 0)
        {
            dataSetTwoReadArray[j] = readArray[i];
            j++;
        }
    }
}

// Function to print dataset one: (must be called after createDataSetOneAndTwo())
void FASTA_readset::printDataSetOne()
{   
    print("\nReads in dataset 1:\n");
    for(int i = 0; i < uniqueSequenceArray[0]; i++)
    {
        print("Read number ", i + 1, " in dataset 1: ", dataSetOneReadArray[i], "\n");
    }
}

// Function to print dataset one and two: (must be called after createDataSetOneAndTwo())
void FASTA_readset::printDataSetTwo()
{   
    print("\nReads in dataset 2:\n");
    for(int i = 0; i < uniqueSequenceArray[1]; i++)
    {
        print("Read number ", i + 1, " in dataset 2: ", dataSetTwoReadArray[i], "\n");
    }
}

// Function to count the length of the header: (input: pointer to header string, output: length of header)
int FASTA_readset::getHeaderLength(char* header)
{
    int headerLength;
    headerLength = 0;
    // Increment header length until terminator (end of string) is reached:
    while(header[headerLength] != '\0') 
    {
        headerLength++;
    }
    return headerLength;
}

// Function to count the unique number of sequences per data set: (for each of the 14)
void FASTA_readset::computeUniqueSequenceCount()
{
    int count = 0;                  
    int unique;                       
    for (int i = 1; i < countCol; i++)
    {  
        unique = 0;                    
        for (int j = 0; j < row; j++) 
        {
            if (countArray[j][i] != 0) 
            {
                unique++;
            }
        }
        uniqueSequenceArray[count] = unique;
        count++;                          
    }
}

// Function to count the total number of sequences per data set: (for each of the 14)
void FASTA_readset::computeTotalSequenceCount()
{
    int count = 0;                 
    int total;                    
    for (int i = 1; i < countCol; i++)  
    {                                   
        total = 0;                   
        for (int j = 0; j < row; j++)
        {
            total = total + (countArray[j][i]); 
        }                                       
        totalSequenceArray[count] = total; 
        count++; 
    }
}

// Function to get the frequencies/counts of each letter (A, C, G, T, N) in the read array (entire data file):
void FASTA_readset::countLetterFrequency()
{
    int frequencyA(0), frequencyC(0), frequencyG(0), frequencyT(0), frequencyN(0);
    // int frequencyA = frequencyC = frequencyG = frequencyT = frequencyN = 0;
    for (int i = 0; i < readCol; i++) // Iterate through all the columns
    { 
        for (int j = 0; j < row; j++) // Iterate through all the rows
        {
            switch(readArray[j][i])
            {
                case 'A': frequencyA++; break;
                case 'C': frequencyC++; break;
                case 'G': frequencyG++; break;
                case 'T': frequencyT++; break;
                case 'N': frequencyN++; break;
                default: print("Invalid Sequence: Acceptable characters include A, C, G, T or N only.\n");
                         print("Exiting the program!\n"); exit(-1);                               
            }
        }
    }
    countA, countC, countG, countT, countN = frequencyA, frequencyC, frequencyG, frequencyT, frequencyN; 
}

// Function to compute statistics (calls other functions, based on problemID):
void FASTA_readset::computeStatistics(char problemID)
{
    if(problemID == 'A')
    {
        computeUniqueSequenceCount();
        computeTotalSequenceCount();
        print("\nNumber of unique sequence fragments in each data set:\n");
        printUniqueSequenceCount();
        print("\nNumber of total sequence fragments in each data set:\n");
        printTotalSequenceCount();
        // countLetterFrequency();
    }
    else if(problemID == 'B')
    {
        computeUniqueSequenceCount();
        createDataSetOneAndTwo();
        auto start = high_resolution_clock::now();    
        int f, count = 0;
        printf("\n");        
        // Search for an element in dataset one against the elements of dataset two, then loop through dataset one elements to search for all of one in two:
        // #pragma omp parallel for private(i, count) shared(dataSetTwoReadArray, dataSetOneReadArray) schedule(static)       
        for(int i = 0; i < uniqueSequenceArray[0]; i++)
        {
            f = searchFunction(dataSetTwoReadArray, dataSetOneReadArray[i]);
            if(f != -1)
            {
                count++;
                print("Read Sequence number ", i, " of dataset one found to match with that in dataset two at ", f, " position.\n");
            }    
        }
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        print("\nTime taken for search: ", duration.count(), " microseconds.\n");  
        print("\nNumber of sequence fragments in dataset one that match with those (unsorted 50 nucleotides-long data) in dataset two: ", count, "\n");       
    }
    else if(problemID == 'C')
    {
        computeUniqueSequenceCount();
        createDataSetOneAndTwo();
        auto sortStartTime = high_resolution_clock::now();        
        sortSecondDataSetReadArray(); // printDataSetTwo();
        auto sortStopTime = high_resolution_clock::now();
        auto sortDuration = duration_cast<microseconds>(sortStopTime - sortStartTime);
        print("\nTime taken for sorting dataset two: ", sortDuration.count(), " microseconds.\n");        

        auto searchStartTime = high_resolution_clock::now();
        int f, count = 0;
        printf("\n");        
        // Perform a binary search for an element in dataset one against the elements of dataset two, then loop to search for all of one in two:
        // #pragma omp parallel for private(i, count) shared(dataSetTwoReadArray, dataSetOneReadArray) schedule(static)
        for(int i = 0; i < uniqueSequenceArray[0]; i++)
        {   
            f = binarySearchFunction(dataSetTwoReadArray, dataSetOneReadArray[i]);
            if(f != -1)
            {
                count++;
                print("Read Sequence number ", i, " of dataset one found to match with that in dataset two at ", f ," position.\n"); 
            }   
        }               
        auto searchStopTime = high_resolution_clock::now();
        auto searchDuration = duration_cast<microseconds>(searchStopTime - searchStartTime);
        print("\nTime taken for search (sorted fragments): ", searchDuration.count(), " microseconds.\n");
        print("\nNumber of sequence fragments in dataset one that match with those (sorted 50 nucleotides-long data) in dataset two: ", count, "\n");
    }
}

// Function that implements a bruteforce search to find a row among multiple rows in a dataset, or in between datasets:
// Complexity: O(50*N) for one fragment, O(50*N*M) for entire search. (N, M = fragments in dataset two, one)
int FASTA_readset::searchFunction(char** DataSetTwoArrayOfArrays, char* DataSetOneSingleArray)
{   
    /* std::strcmp() version:
    for(int i = 0; i < uniqueSequenceArray[1]; i++)
    {
        // print("Comparing: ", DataSetTwoArrayOfArrays[i], " and ", DataSetOneSingleArray, " for this round.\n");
        if((std::strcmp(DataSetTwoArrayOfArrays[i], DataSetOneSingleArray) == 0)) 
            return 1;
    }*/           
    for (int i = 0; i < readCol; i++)
    { 
        for (int j = 0; j < uniqueSequenceArray[1]; j++)
        {
            int counter = 0;
            for(int j = 0; j < 50; j++)
            {   
                if(DataSetTwoArrayOfArrays[j][i] == DataSetOneSingleArray[j]) // Comparing character by character
                    counter++;
            }
            if(counter == 50) // If all 50 characters/nucleotides are the same, then its a match! (element found)
                return 1;                             
        }
    }   
    return -1;
}

// Function that implements a binary search to search for one read among multiple rows in a dataset, or in between datasets:
// Use elements of dataset one to search against elements of dataset two (or search for dataset one elements in readArray)
int FASTA_readset::binarySearchFunction(char* DataSetTwoArrayOfArrays[], char* DataSetOneSingleArray)
{  
    int left = 0, right = uniqueSequenceArray[1] - 1;
    while (left <= right) 
    {
        int middle = left + (right - left) / 2;
        if (DataSetTwoArrayOfArrays[middle] == DataSetOneSingleArray)
            return middle;
        if (DataSetTwoArrayOfArrays[middle] < DataSetOneSingleArray)
            left = middle + 1;
        else
            right = middle - 1;
    }
    return -1;
}

// Function to print the count array: (must be called after saveData())
void FASTA_readset::printCountArray()
{
    for (int i = 0; i < row; i++) 
    {
        for (int j = 0; j < countCol; j++) 
        {
            print(countArray[i][j], " ");
        }
        print("\n");
    }
}

// Function to print the read array: (must be called after saveData())
void FASTA_readset::printReadArray()
{
    for (int i = 0; i < row; i++) 
    {
        print(i, " ");
        for (int j = 0; j < readCol; j++) 
        {
            print(readArray[i][j]);
        }
        print("\n");
    }
}

// Function to print the number of unique sequences in each of the 14 data sets: (must be called after computeUniqueSequenceCount())
void FASTA_readset::printUniqueSequenceCount()
{
    for (int i = 0; i < 14; i++) 
    {
        print("There are ", uniqueSequenceArray[i], " unique sequence fragments in dataset number ", i + 1, ".\n");
    }
}

// Function to print the total number of sequences in each of the 14 data sets: (must be called after computeTotalSequenceCount())
void FASTA_readset::printTotalSequenceCount()
{
    int length = 14;  // length of totalSequenceArray
    for (int i = 0; i < length; i++) 
    {
        print("There are ", totalSequenceArray[i], " total sequence fragments in dataset number ", i + 1, ".\n");
    }
}

// Function to print the total count of each letter in entire data file: (must be called after countLetterFrequency()):
void FASTA_readset::printLetterCounts()
{
    print("There are a total of ", countA, " 'A' characters in the data.\n");
    print("There are a total of ", countC, " 'C' characters in the data.\n");
    print("There are a total of ", countG, " 'G' characters in the data.\n");
    print("There are a total of ", countT, " 'T' characters in the data.\n");
    print("There are a total of ", countN, " 'N' characters in the data.\n");
}

// Function to sort the sequence fragment array for dataset two (entire rows, not characters within one): (simple bubble sort)
void FASTA_readset::sortSecondDataSetReadArray()
{   
    char tempArray[readCol];
    int size = uniqueSequenceArray[1];
    for(int i = 1; i < size; i++)
    {
        for(int j = 0; j < size - i; j++)
        {
            if(strcmp(dataSetTwoReadArray[j], dataSetTwoReadArray[j + 1]) > 0)
            {
                strcpy(tempArray, dataSetTwoReadArray[j]);
                strcpy(dataSetTwoReadArray[j], dataSetTwoReadArray[j + 1]);
                strcpy(dataSetTwoReadArray[j + 1], tempArray);
            }
        }
    }
    print("\nStatus update: The sort of sequence fragments (reads) for dataset two has been successful.\n");    
}                

// Function to swap two rows in the read array: (input: indices of the rows)
// Alternatives: std::swap(k, j), XOR, x = x + y - (y = x) with ascii shifts
void FASTA_readset::swap(int j, int k)
{
    char tempArray[readCol];
    for (int i = 0; i < readCol; i++)
    {
        tempArray[i] = readArray[j][i];
        readArray[j][i] = readArray[k][i];
        readArray[k][i] = tempArray[i];
    }
}

// Function to sort the sequence fragment arrays for the entire read array (again, for the entire 50-character long rows): (selection)
void FASTA_readset::sortReadArray()
{ 
    int p;
    for (int k = 0; k < row - 1; k++) 
    {   
        for (int j = k + 1; j < row; j++) 
        {
            p = 0; 
            while (p < readCol)
            {
                if (readArray[j][p] > readArray[k][p]) 
                {
                    break;
                } 
                else if (readArray[j][p] < readArray[k][p]) 
                {
                    swap(k, j);
                    break;
                } 
                else if (readArray[j][p] == readArray[k][p]) 
                {
                    p++;
                } 
                else 
                {
                    print("\nError!\n");
                    exit(-1);
                }
            }
        }
    }
    print("\nStatus update: The sort of sequence fragments for the read array has been successful.\n");
}

// Default constructor:
FASTA_readset::FASTA_readset()
{    
    // Initialize everything to 0/NULL:
    countCol = readCol = row = 0;
    // Could use the default values for variables countCol (15) and readCol (50), then a random value for row above
    // following with the initialization of the array elements to zero below, but it seems like a waste of memory.
    countArray = new int*[row];
    readArray = new char*[row];
    for(int i = 0; i < row; ++i) 
    {
        countArray[i] = new int[countCol];
        readArray[i] = new char[readCol];
    }
}

// Custom constructor:
FASTA_readset::FASTA_readset(int rowCount, char* filePathArray)
{
    filePath = filePathArray;  // Save filepath from argv[2]
    countCol = 15;            // 1 identifier (R prefixed) + 14 copy numbers = 15 characters
    readCol = 50;            // 50 characters (nucleotides) long
    row = rowCount;         // Set rows to consider/read based on user input
    int allRows = countNumberOfReads();
    // varDebug(allRows);
    // If custom constructor called at main has rowCount < 0, then consider all rows to read, else consider the user provided value:
    row = (rowCount < 0) ? allRows : rowCount;
    // Now that row count is taken care of, I allocate memory and initialize the required arrays:
    countArray = new int*[row];
    readArray = new char*[row];
    // Made a simple parallelized version (just to test faster) using OpenMP: (for instance, each thread here gets rows given by row/threadCount)
    // Add stuff to:
    // Homework1.cpp: #include <omp.h> 
    // Makefile: -fopenmp -lpthread $threadCount
    // Runner scripts: #SBATCH --cpus-per-task threadCount
    // tried for threadCount := 2, 4 and 8; (number of threads/cores while testing on Monsoon, i.e. 1 thread per CPU core therein)
    // #pragma omp parallel for private(i) shared(countArray, readArray) schedule(static) // collapse(threadCount)
    for(int i = 0; i < row; i++) 
    {
        countArray[i] = new int[countCol];
        readArray[i] = new char[readCol];
    }
}

// Destructor:
FASTA_readset::~FASTA_readset()
{
    // Uncomment the commented lines (with code) below to know the time taken by the destructor:
	// auto startDestruction = high_resolution_clock::now();
    // Delete the arrays used, thereby deallocating the memory reserved for them:
    delete[] countArray, readArray, uniqueSequenceArray, totalSequenceArray, dataSetOneReadArray, dataSetTwoReadArray;
    // Delete rest of the variables for the same purpose:
    delete[] filePath, countA, countC, countG, countT, countN, readCol, row;                                             
    // auto wrapUp = high_resolution_clock::now();
    // auto duration1 = duration_cast<microseconds>(startDestruction - wrapUp);
    // print("\nTime taken by the Destructor: ", duration1.count(), " microseconds.\n");
}

int main(int argc, char *argv[])
{
    if (argc != 3) 
	{
		print("Error: Two input parameters are expected.\nProper usage:\n./Sequence <problemFlag> <filePath>\nExiting the program!\n");  
		exit(-1); // equivalent, but return 1/0; might be more graceful
	}
    else print("The number of arguments passed is: ", argc, "\nThe first argument is: ", argv[0], 
	          "\nThe second argument is: ", argv[1], "\nThe third argument is: ", argv[2], "\n"); 
              
    std::string problemNumber = argv[1], filePath = argv[2];
    char* filePathArray = new char[100];
    char* problemNumberArray = new char[100];

    // Saving the length of problemNumber and filePath, since I'm dealing with an std::string object instead of a const char*:
    int problemNumberlength = 0, filePathLength = 0;
    while (problemNumber[problemNumberlength] != '\0')
    {
        problemNumberlength++;
    }
    while (filePath[filePathLength] != '\0')
    {
        filePathLength++;
    }
    // Saving problemNumber and filePath as arrays:
    for(int i = 0; i < problemNumberlength; i++)
    {
        problemNumberArray[i] = problemNumber[i];
    }
    for(int i = 0; i < filePathLength; i++)
    {
        filePathArray[i] = filePath[i];
    }

    // Saving the problem ID (A, B or C) character to compare against, and then passing it as argument to my custom constructor:
    char problemID = problemNumberArray[8]; // 9th character, right after 'Problem1'

    switch(problemID)
    {
        case 'A': print("\nProblem1A selected: Computing the number of unique and total sequence fragments!\n"); break;
        case 'B': print("\nProblem1B selected: Comparing the contents of data sets 1 and 2! (unsorted data)\n"); break;
        case 'C': print("\nProblem1C selected: Comparing the contents of data sets 1 and 2! (sorted data)\n");   break;
         default: print("\nInvalid second argument. Expecting one of these: Problem1A or Problem1B or Problem1C\n");
                  print("Exiting the program!\n"); exit(-1);          
    }    

    // Common for all problem cases, with the specific problemID (A/B/C) specifying which statistics to compute:
    // Change to the row count you want, if you do not want to read the entire datafile (path to which you're passing as the third argument):    
    int rowCount = -1; // Going with the default, i.e. considering all the rows of the combined dataset
    // fastio;
    FASTA_readset test(rowCount, filePathArray);
    test.saveData();
    // Perform all the desired actions for a particular problemID: (inclusive of computing and printing statistics)
    test.computeStatistics(problemID);    

    return 0;
}
