/*-------------------------
  Author: Anirban
  Email:  ac4743@nau.edu
--------------------------*/
// Compile and run locally using: 
// make && ./<executablename> -o <problemname> <genomefilepath>
// Examples:
// make && ./hashTableConstructs Problem1A sampleGenome.fasta
// Usage on a Slurm-based compute cluster: make && sbatch <runnerscript>.sh
#include "FASTAreadset_HT.h"

int main(int argc, char* argv[]) 
{   // Argument count must be three: (executable name followed by problem ID and the genome file path)
    if(argc != 3) 
    {
        print("\nError: Two or three input parameters are expected.\n\nProper usage looks like:\n", "./hashTableConstructs <problemID> <filepath>\n", 
		"For example, try this:\n./hashTableConstructs Problem1A /scratch/ac4743/HW2/test_genome.fasta\n", "\n\nExiting the program!\n"); 
        exit(-1);
    }
    else 
    {
        print("\nThe number of arguments passed is: ", argc, "\nThe first argument is: ", argv[0], 
	          "\nThe second argument is: ", argv[1], "\nThe third argument is: ", argv[2], "\n"); 
    }

    // fastio;
    std::string problem = argv[1], genomeFilePath = argv[2];
    // varDebug(problem[8]);

    if(!problem.compare("Problem1A")) 
    {   // Creating an object of my hash table class to store the genome sequences (using a table size of five million):
        FASTAreadset_HT requiredObject(genomeFilePath, 5000000);        
        char** genome = requiredObject.saveGenomeFile(genomeFilePath);
        requiredObject.saveGenomeFileInHashTable();  

        // Generating five million randomized genome-based 16-mers and searching for them in my hash table populated with the Anthrax genome:      
        unsigned int matchCount = requiredObject.randomGenomeSequenceComparison(genome);
        print("Number of matches in between the one million 16-mers randomly picked from the genome versus and the sequences from the hash table: ", matchCount, " \n");

        // Generating five million completely random (not from the genome) 16-mers and searching for them in my hash table populated with the Anthrax genome:
        matchCount = requiredObject.randomSequenceComparison();
        print("Number of matches in between the one million completely random (not from the genome) 16-mers versus the sequences from the hash table: ", matchCount, " \n");        
        for(int i = 0; i < requiredObject.genomeRowCount; i++) delete[] genome[i];
        delete genome;
    }
    else if(!problem.compare("Problem1B")) 
    {   // Setting the size of hash table to ten million elements, then creating an object using that value of m passed to my constructor:
        unsigned int m = 10000000, matchCount;
        FASTAreadset_HT requiredObject(genomeFilePath, m);        
        char** genome = requiredObject.saveGenomeFile(genomeFilePath);
        requiredObject.saveGenomeFileInHashTable();  

        // Iterating through my genome and searching for the 16-mer sequences in my hash table:    
        matchCount = requiredObject.genomeSearch(genome);
        print("\nNumber of matches for the genome's 16-mers and the entries in the hash table: ", matchCount, " \n");
        
        // Iterating through my genome and searching for the chance-based modified 16-mer sequences in my hash table:
        matchCount = requiredObject.genomeSearchWithOnePercentBaseErrorRate(genome);
        print("\nNumber of matches for the genome's 16-mers (after undergoing 1% per-base error inclusion for all the characters) and the entries in the hash table: ", matchCount, " \n");
        for(int i = 0; i < requiredObject.genomeRowCount; i++) delete[] genome[i];
        delete genome;
    }        
    else 
    {
        print("\nInvalid second argument. Expecting one of: 'Problem1A' or 'Problem1B' or 'ProblemC'\n");
        exit(-1);
    }
}
