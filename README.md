MPI_Deterministic_Sample_Sort
=============================

##Assignment 4 COMP 4009


Write a parallel MPI program for p processors to sort n integers via the  Deterministic Sample Sort method discussed in class. For the local sequential array sort, use in-place Heapsort (array implementation of the heap). For communication, use only MPI_AllToAll or MPI_AllToAllv; no single message passing such as MPI_send.

Each processor should read an input file of n/p integers in the following format. The input file for processor i should be named "input-i.txt" (input-00.txt, input-01, input-02, input-03, ...). Following the sort, each processor will output a file "output-i.txt" with its n/p integers of the result (one integer per line). 
