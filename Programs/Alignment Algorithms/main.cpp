/*-------------------------
  Author: Anirban166/Ani
  Email:  ac4743@nau.edu
--------------------------*/
// Compile and run locally using: 
// make && ./<executablename> -o <problemname> <readsetfilepath> <genomefilepath>
// Example: make && ./alignmentConstructs Problem1A sampleDataset.fa sampleGenome.fasta
// Usage on a Slurm-based compute cluster: make && sbatch <runnerscript>.sh
#include "BLAST.h"

int main(int argc, char* argv[]) 
{
    if(argc != 4) 
    {
        print("\nError: Three input parameters (post the executable name) are expected.\n\nProper usage looks like:\n", "\n./alignmentConstructs <problemID> <filepath> <filepath>\n\n", 
              "For example, try this:\n./alignmentConstructs Problem1A /scratch/ac4743/sampleDataset.fa /scratch/ac4743/sampleGenome.fa\n\nExiting the program!\n"); 
        exit(-1);
    }    
    else print("\nThe number of arguments passed is: ", argc, "\nThe first argument is: ", argv[0], "\nThe second argument is: ", argv[1], "\nThe third argument is: ", argv[2], "\nThe fourth argument is: ", argv[3], "\n");

    // fastio;
    std::string problem = argv[1], readFilePath = argv[2], genomeFilePath = argv[3];
    // varDebug(problem[8]);
    int* max = new int;
    char* genome;
    int genomeLength, readCount;

    if(!problem.compare("Problem1A")) 
    {
        // Instantiating an object of my class with a match score of 2, a mismatch penalty of -1 and a gap penalty of -1 as required:
        BLAST testObject(2, -1, -1);
        // Using 11-mer seeds and a hash table of size (m) 10,000,000 as required:
        int kmer = 11, hashTableSize = 10000000, matchCount = 0;
        char** readset = nullptr;              // Array for readset
        int* genomeIndex = new int;           // Pointer for indexing genome
        char** genomeKmer;                   // Array for genome splits into kmer
        int* genomeRowCount = new int;      // Number of rows in the genome kmer
        char* subGenome = nullptr;         // To expand portion of genome for seed search
        int* subGenomeLength = new int;   // Length of the sub-genome
        bool* seedFound = new bool;      // Boolean indicating whether the seed was found or not
        bool printFlag = false;         // Boolean indicating whether to print text alignment or not
        Node** T;                      // Pointer to hash table
        int percentError;             // Percent of error (to be set 5% to work on the second part of 1A)

        // Saving the genome from its file:
        genomeLength = testObject.getLengthOfGenome(genomeFilePath);
        genome = new char[genomeLength];
        testObject.saveGenome(genomeFilePath, genome, genomeLength);
        // Splitting it into a kmer array:
        genomeKmer = testObject.splitFragmentIntoKmer(genome, genomeLength, kmer, genomeRowCount);
        print("Number of rows in genome: ", *genomeRowCount, "\n");
        // Storing that genome kmer array in the hash table:
        T = testObject.createHashTable(hashTableSize);
        for(int i = 0; i < *genomeRowCount; i++) 
        {
            *genomeIndex = i;
            testObject.insertFragmentInTable(genomeKmer[i], kmer, hashTableSize, T, genomeIndex);
        }      

        // Randomly selecting five million 50-mer fragments from the genome (Anthrax):
        readCount = 5000000; 
        matchCount = 0;
        readset = testObject.selectNRandomFragmentsFromGenome(readCount, testObject.readLength, genome, genomeLength);
        // Performing seed based alignment for each fragment from the BLAST class object:
        auto startTime = high_resolution_clock::now();         
        for (int fragment = 0; fragment < readCount; fragment++) 
        {
            subGenome = testObject.searchAndExpandSeed(fragment, readset, testObject.readLength, kmer, hashTableSize, genomeIndex, T, genomeLength, genome, subGenomeLength, seedFound, subGenome, printFlag);
            if (*seedFound) 
            {
                matchCount++;
                textAlignment* obj = testObject.BLASTAlignment(subGenome, readset[fragment], *subGenomeLength, testObject.readLength, testObject.matchScore, testObject.gapScore, testObject.mismatchScore, max);
            }
        }        
        print("Number of 50-mer fragments considered: ", readCount);
        print("Number of matches: ", matchCount, ".\n");
        auto stopTime = high_resolution_clock::now();
        auto duration = duration_cast<seconds>(stopTime - startTime);
        print("\nTime taken for querying the genome: ", duration.count(), " seconds.\n");         

        // Adding 5% per-base random errors and re-selecting for the second part of Problem1A:
        percentError = 5;

        // Randomly selecting one million 50-mer fragments from the genome (Anthrax) with 5% per-based error rate:
        readCount = 1000000; 
        matchCount = 0; // Resetting the counter (using the same variable to represent the number of matches)
        testObject.introducePercentError(readCount, testObject.readLength, readset, percentError);
        // Performing seed based alignment for each fragment from the BLAST class obejct:
        for(int fragment = 0; fragment < readCount; fragment++) 
        {
            subGenome = testObject.searchAndExpandSeed(fragment, readset, testObject.readLength, kmer, hashTableSize, genomeIndex, T, genomeLength, genome, subGenomeLength, seedFound, subGenome, printFlag);
            if(*seedFound) 
            {
                matchCount++;
                textAlignment* obj = testObject.BLASTAlignment(subGenome, readset[fragment], *subGenomeLength, testObject.readLength, testObject.matchScore, testObject.gapScore, testObject.mismatchScore, max);
            }
        }
        print("Number of 50-mer fragments considered: ", readCount);        
        print("Number of matches with 5% error rate: ", matchCount, ".\n");

        // Deallocate memory for used pointers/objects:
        testObject.deleteHashTable(hashTableSize, T);
        delete[] T;
        delete max;
        delete[] genome;        
        for(int i = 0; i < readCount; i++) 
        {
            delete[] readset[i];
        }
        delete[] readset;        
        for(int i = 0; i < *genomeRowCount; i++) 
        {
            delete[] genomeKmer[i];
        }
        delete[] genomeKmer;
        delete seedFound;        
        delete genomeIndex;
        delete genomeRowCount;
        delete[] subGenome;
        delete subGenomeLength;
    }
    else if(!problem.compare("Problem1B"))
    {
        char** readset;
        int alignmentLength;        
        BLAST testObject(2, -1, -1);

        // Saving the genome and reads:
        genomeLength = testObject.getLengthOfGenome(genomeFilePath);
        genome = new char[genomeLength];
        testObject.saveGenome(genomeFilePath, genome, genomeLength);
        readCount = testObject.countNumberOfReads(readFilePath);
        readset = new char* [readCount];
        for(int i = 0; i < readCount; i++) 
        {
            readset[i] = new char[testObject.readLength];
        }
        testObject.saveReads(readCount, testObject.readLength, readset, readFilePath);
        print("\n");

        // Variable to store the number of perfect hits: (i.e., when all 50 characters are a match, giving a score of 100 for a match score of +2)
        int perfectHitsCounter = 0;
        auto startTime = high_resolution_clock::now(); 
        // Perform my alignment algorithm for each read:
        for(int i = 0; i < readCount; i++) 
        {
            textAlignment* obj = testObject.BLASTAlignment(genome, readset[i], genomeLength, testObject.readLength, testObject.matchScore, testObject.gapScore, testObject.mismatchScore, max);
            alignmentLength = obj->length;
            print("Printing alignment for read #", (i + 1), ":\n");
            obj->printAlignment(alignmentLength);
            if((*max) == 100) perfectHitsCounter++;            
            print("Maximum score: ", *max, "\n");
        }
        print("\nPerfect hits: ", perfectHitsCounter, ".\n");
        auto stopTime = high_resolution_clock::now();
        auto duration = duration_cast<seconds>(stopTime - startTime);
        print("\nTime taken for querying the genome stored in a BLAST object against the entire reads dataset: ", duration.count(), " seconds.\n");          

        delete max;
        delete[] genome;        
        for(int i = 0; i < readCount; i++) 
        {
            delete[] readset[i];
        }
        delete[] readset;
    }
    else 
    {
        print("Invalid problem ID ('Problem1A' or 'Problem1B' are the valid ones).\nExiting!\n");
        exit(-1);
    }
}