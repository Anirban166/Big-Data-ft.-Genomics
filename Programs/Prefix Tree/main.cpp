/*-------------------------
  Author: Anirban166/Ani
  Email:  ac4743@nau.edu
--------------------------*/
// Compile and run locally using: 
// make && ./<executablename> -o <problemname> <genomefilepath>
// Example: make && ./prefixTreeConstructs Problem1A SARSCoV2.txt
// Usage on a Slurm-based compute cluster: make && sbatch <runnerscript>.sh
#include "prefixTree.h"

int main(int argc, char* argv[]) 
{
    if(argc != 3) 
    {
        print("\nError: Three input parameters (post the executable name) are expected.\n\nProper usage looks like:\n", "\n./prefixTreeConstructs <problemID> <filepath>\n\n", 
              "For example, try this:\n./prefixTreeConstructs Problem1A /scratch/ac4743/SARSCoV2.txt\n\nExiting the program!\n"); 
        exit(-1);
    }    
    else print("\nThe number of arguments passed is: ", argc, "\nThe first argument is: ", argv[0], "\nThe second argument is: ", argv[1], "\nThe third argument is: ", argv[2], "\nThe fourth argument is: ", argv[3], "\n");

    // fastio;
    std::string problem = argv[1], genomeFilePath = argv[2];
    // varDebug(problem[8]);

    if(!problem.compare("Problem1A"))
    { 
        // Considering ten million random 36-mers:
        int readCount = 10000000, totalMatchCount, readlength = 36 + 1; // 36-mer + one for the terminating character ($)
        print("\nSearching the prefix tree for ", readCount, " random reads: \n");
        prefixTree testObject(genomeFilePath, readCount, readlength);
        testObject.createTrieReadset();
        totalMatchCount = testObject.searchTrieForGenomeKmers();
        print("Total number of matches (including multiples): ", totalMatchCount, "\n");
    }
    else if(!problem.compare("Problem1B"))
    {
        // Considering ten million random 36-mers each with a 5% per-base error rate:
        int readCount = 10000000, totalMatchCount, readlength = 37, percentError = 5;
        print("\nSearching the prefix tree for ", readCount, " random reads: (each with a ", percentError, "% per-base error rate)\n");
        prefixTree testObject(genomeFilePath, readCount, readlength, percentError);
        testObject.createTrieReadset();
        totalMatchCount = testObject.searchTrieForGenomeKmers();
        print("Total number of matches (including multiples): ", totalMatchCount, "\n");        
    }
    else 
    {
        print("Invalid problem ID ('Problem1A' or 'Problem1B' are the valid ones).\nExiting!\n");
        exit(-1);
    }
}