/*-------------------------
  Author: Anirban166/Ani
  Email:  ac4743@nau.edu
--------------------------*/
// Compile and run locally using: 
// make && ./<executablename> -o <problemname> <genomefilepath>
// Example: make && ./suffixTreeConstructs Problem1A SARSCoV2.txt
// Usage on a Slurm-based compute cluster: make && sbatch <runnerscript>.sh
#include "suffixTree.h"

int main(int argc, char* argv[]) 
{
    if(argc != 3) 
    {
        print("\nError: Three input parameters (post the executable name) are expected.\n\nProper usage looks like:\n", "\n./suffixTreeConstructs <problemID> <filepath>\n\n", 
              "For example, try this:\n./suffixTreeConstructs Problem1A /scratch/ac4743/SARSCoV2.txt\n\nExiting the program!\n"); 
        exit(-1);
    }    
    else print("\nThe number of arguments passed is: ", argc, "\nThe first argument is: ", argv[0], "\nThe second argument is: ", argv[1], "\nThe third argument is: ", argv[2], "\nThe fourth argument is: ", argv[3], "\n");

    // fastio;
    std::string problem = argv[1], genomeFilePath = argv[2];
    // varDebug(problem[8]);

    if(!problem.compare("Problem1A"))
    {
        // Considering ten million random 36-mers:        
        int matchCount, readCount = 10000000, readLength = 36;
        suffixTree testObject(genomeFilePath);
        testObject.createGenomeSuffixTree();
        print("\nSearching the suffix tree for ", readCount, " reads: \n");
        matchCount = testObject.searchSuffixTrieFor_N_RandomReadsFromGenome(readCount, readLength);
        print("Total number of matchCount (reads found): ", matchCount, "\n");
    }
    else if(!problem.compare("Problem1B"))
    {
        char* sequence = "ACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        suffixTree test(genomeFilePath);
        test.createGenomeSuffixTree();
        test.searchSuffixTreeForRead(sequence, 36);       
    }    
    else 
    {
        print("Invalid problem ID ('Problem1A' or 'Problem1B' are the valid ones).\nExiting!\n");
        exit(-1);
    }
}
