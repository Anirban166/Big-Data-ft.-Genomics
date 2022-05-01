General
---
Venturing into big data via various data structures and algorithms for analysis, organization and querying of large genomic datasets using C++, Slurm and a compute cluster.

Contents
---
- [Array Constructs](https://github.com/Anirban166/Big-Data-ft.-LSDS/tree/main/Programs/Array%20of%20Arrays): Using an array of arrays to store, sort and search through 50-character sequence fragments from a dataset in FASTA format with over 36 million reads.
- [Linked List Constructs](https://github.com/Anirban166/Big-Data-ft.-LSDS/tree/main/Programs/Linked%20List%20of%20Arrays): Using an unidirectional singly linked list of character arrays to store and search through sequence fragments, with these in turn being used as queries against 50-mer sequences of the Bacillus anthracis genome, with over 36 million and 5.2 million entries in the read set and genome file respectively.
- [Hash Table Constructs](https://github.com/Anirban166/Big-Data-ft.-LSDS/tree/main/Programs/Hash%20Table%20with%20Collision%20Chaining): Using a hash table (with collision-based chaining) to store 16-mer sequences of the Anthrax genome, with these in turn being compared against millions of randomized genome-based 16-mers, completely randomized equivalents (i.e., ones not from the genome) and modified 16-mers with changes in the base characters given by a 1% per-base error rate. 
- [Alignment Algorithms](https://github.com/Anirban166/Big-Data-ft.-LSDS/tree/main/Programs/Alignment%20Algorithms): Using a chaining hash table with an implementation of the BLAST algorithm for aligning k-mer fragments of a genome (such as SARS-CoV-2 and Anthrax) and for comparing against sequences from a readset (to find number of matches and perfect hits), with the former being stored as entries in the table (comparison is also done for the genome having a 5% change of modification for each base). Includes seed-based local alignment and a contrasting human-friendly output (computed and organized via tracebacks) for the result with respect to the subject sequence.
- [Prefix Tree Constructs](https://github.com/Anirban166/Big-Data-ft.-LSDS/tree/main/Programs/Prefix%20Tree): Using a trie to store 36-mers from the SARS-CoV-2 genome and using both regular and error induced random sequences (ten million of each) to search for matches against the stored entries.
- [Suffix Tree Constructs](https://github.com/Anirban166/Big-Data-ft.-LSDS/tree/main/Programs/Suffix%20Tree): Using a suffix tree to store 36-mer sequences from the SARS-CoV-2 genome and using randomized versions of these (again, in millions) as search queries against the entries stored in the tree, apart from providing the ability to search the stored entries for a perfect match with a particular sequence of user-specified length.
  
Data
--- 
Sample datasets to test my code against can be found [here](https://github.com/Anirban166/Big-Data-ft.-LSDS/tree/main/Sample%20Data). Files that exceed a hundred megabytes have been brought under the bar with a small subset of the original data being extracted, as the full versions for those go upto several gigabytes of data, and would need paid use of GitHub's large file system.
