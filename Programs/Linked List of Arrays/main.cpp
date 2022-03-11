/*-------------------------
  Author: Anirban
  Email:  ac4743@nau.edu
--------------------------*/
// Compile and run locally using: 
// make && ./<executablename> -o <problemname> <readsetfilepath> <genomefilepath>
// Examples:
// make && ./linkedListConstructs Problem1A sampleDataset.fa
// make && ./linkedListConstructs Problem1B sampleDataset.fa sampleGenome.fasta
// Monsoon (or the real deal): make && sbatch <runnerscript>.sh
#include "FASTAreadset_LL.h"

int main(int argc, char* argv[]) 
{   // Argument count must be three or four: (for problem part A or B respectively)
    if((argc < 3) || (argc > 4)) 
    {
        print("\nError: Two or three input parameters are expected.\n\nProper usage looks like:\n",
        "./Homework2 <problemID> <filepath>\n", "or\n./Homework2 <problemID> <filepath> <filepath>\n\n", 
		"For example, try this:\n./Homework2 Problem1A /scratch/ac4743/HW2/sample_dataset.fa\n", 
        "or\n./Homework2 Problem1B /scratch/ac4743/HW2/sample_dataset.fa /scratch/ac4743/HW2/sample_genome.fa", 
        "\n\nExiting the program!\n"); 
        exit(-1);
    }
    else 
    {
        print("\nThe number of arguments passed is: ", argc, "\nThe first argument is: ", argv[0], 
	          "\nThe second argument is: ", argv[1], "\nThe third argument is: ", argv[2], "\n"); 
        if(argc == (1 << 2)) // Expressing my predilection for bitwise shifts and powers of two :)
            print("The fourth argument is: ", argv[3], "\n");
    }

    // fastio;
    std::string problem = argv[1], readSetFilePath = argv[2];
    // varDebug(problem[8]);

    if(!problem.compare("Problem1A"))
    {
        print("\nProblem1A selected: Saving the entire ~36 million read set in an object of FASTAreadset_LL!\n");
        FASTAreadset_LL test(readSetFilePath);
        test.saveReadsLL();      
        // Use copy constructor if required:
        // print("\nCopying the entire read set:\n");
        // FASTAreadset_LL copiedObject(test);

        Node* match = NULL;
        std::string query;

        print("\nBeginning the search for the required five fragments!\n");
        // Search for sequence (i):
        query = "CTAGGTACATCCACACACAGCAGCGCATTATGTATTTATTGGATTTATTT";
        match = test.searchFragment(query);
        test.printSearchResults(match, query);
        // Search for sequence (ii):
        query = "GCGCGATCAGCTTCGCGCGCACCGCGAGCGCCGATTGCACGAAATGGCGC";
        match = test.searchFragment(query);
        test.printSearchResults(match, query);
        // Search for sequence (iii):
        query = "CGATGATCAGGGGCGTTGCGTAATAGAAACTGCGAAGCCGCTCTATCGCC";
        match = test.searchFragment(query);
        test.printSearchResults(match, query);
        // Search for sequence (iv):
        query = "CGTTGGGAGTGCTTGGTTTAGCGCAAATGAGTTTTCGAGGCTATCAAAAA";
        match = test.searchFragment(query);
        test.printSearchResults(match, query);
        // Search for sequence (v):
        query = "ACTGTAGAAGAAAAAAGTGAGGCTGCTCTTTTACAAGAAAAAGTNNNNNN";
        match = test.searchFragment(query);
        test.printSearchResults(match, query);
    }
    else if(!problem.compare("Problem1B")) 
    {
        char** genome;
        int matchCount;        
        std::string genomeFilePath = argv[3];

        print("\nProblem1B selected: Saving the read set and the genome file!\n");
        FASTAreadset_LL test(readSetFilePath);
        test.saveReadsLL();
        print("\nComputing the count of 50-character fragments that can be obtained from the Bacillus anthracis genome!\n");       
        genome = test.saveGenomeFile(genomeFilePath);
        print("\nStatus update: The read set and the genome file have been saved successfully (as a linked list and as an array respectively).\n");         
        
        // Estimating the total time for 1,000, 10,000 and 100,000 queries:        
        print("\nIterating through the genome to find matching sequences within the read set!\n");
        for(int queryCount = 1000; queryCount < 1000000; queryCount *= 10)
        {
            print("\nSearching for matches in the entire genome versus ", queryCount, " read fragments:\n");
            matchCount = test.genomeSearchMatches(queryCount, genome);
            print("\nNumber of matches for ", queryCount, " reads: ", matchCount, "\n");
        }      
        // Kick genome out of existence:
        for(int i = 0; i < test.genomeRowCount; i++) 
        {
            delete[] genome[i];
        }
        delete genome;
    }
    else 
    {
        print("\nInvalid second argument. Expecting either 'Problem1A' or 'Problem1B'.\n");
        exit(-1);
    }
    return 0;
}