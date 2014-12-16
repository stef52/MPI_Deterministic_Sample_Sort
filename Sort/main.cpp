#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "mpi.h"

#define InputFilePrefix "./inputFiles/"      //input file location 
#define OutputFilePrefix "./outputFiles/"    //output file location
#define MPI_Type  MPI_INT
#define IsSmaller(x, y) ((x) < (y)) 		//is smaller macros

static void fix(int data[], register int m, register int n); //readjust heap
void heapsort(int data[], int n); // heap sort
void sort(int data[], int localSize); //mpi sort

int main(int argc,char *argv[])
{
	int myid, prosCount;
	FILE *infile, *outfile;
	char infilepostfix[7], outfilepostfix[8];
	char infilename[100] = InputFilePrefix;
	char outfilename[100] = OutputFilePrefix;
	int *data;	// input data 
	int n;     	// input size 
	int p;     	// input proc 

	int buf, tmp, i, startTotal, stopTotal;

	MPI_Init(&argc, &argv);

	startTotal = MPI_Wtime();
	MPI_Comm_size(MPI_COMM_WORLD, &prosCount);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	infilepostfix[0] = 'i';
	infilepostfix[1] = 'n';
	infilepostfix[2] = 'p';
	infilepostfix[3] = 'u';
	infilepostfix[4] = 't';
	infilepostfix[5] = '-';
	infilepostfix[6] = ((myid % 100) / 10) + 48;	//found this on google
	infilepostfix[7] = (myid % 10) + 48;
	infilepostfix[8] = '.';
	infilepostfix[9] = 't';
	infilepostfix[10] = 'x';
	infilepostfix[11] = 't';
	infilepostfix[12] = '\0';
	  
	strcat(infilename, infilepostfix);
	infile = fopen(infilename, "r");
  
	if (infile == NULL) 
	{
		printf("can't open infile(s) \n");
		exit(0);
	}
	  
	fscanf(infile, "%d\n", &n);
	fscanf(infile, "%d\n", &p);
	n = n/p; //  n/p
	data = 	((int *) malloc(n*sizeof(int)+1));
	 
	for (i=0; i<n; i++) 
	{
		fscanf(infile, "%d", &buf); 
		data[i]=buf;
	}
  
	tmp = fclose(infile);
	if (tmp != 0) 
	{
		printf("can't close %s \n", infilename);
		exit(0);
	}

	sort(data, n);

	outfilepostfix[0] = 'o';
	outfilepostfix[1] = 'u';
	outfilepostfix[2] = 't';
	outfilepostfix[3] = 'p';
	outfilepostfix[4] = 'u';
	outfilepostfix[5] = 't';
	outfilepostfix[6] = ((myid % 100) / 10) + 48;
	outfilepostfix[7] = (myid % 10) + 48;
	outfilepostfix[8] = '.';
	outfilepostfix[9] = 't';
	outfilepostfix[10] = 'x';
	outfilepostfix[11] = 't';
	outfilepostfix[12] = '\0';
  
	strcat(outfilename, outfilepostfix);
	outfile = fopen(outfilename, "w");
	fprintf(outfile,"proc. %d: \n n = %d \n p = %d \n", myid, n, p);
	 
	for (i=0; i<n; i++) 
	{
		buf = data[i];
		fprintf(outfile,"%d \n", buf);
	}

	tmp = fclose(outfile);  
	if (tmp != 0) 
	{
		printf("can't close outfile %s \n", outfilename);
		exit(0);
	}
	  
	stopTotal = MPI_Wtime();
	printf("Proc. %d: Total Time = %d sec.\n", myid, (stopTotal - startTotal));
	MPI_Finalize();
}

void sort(int data[], int n)
{
	int myid, prosCount, resultSize,i, j, tmp, l, r, Left, Right;
	int *scounts, *displs, *sdispl, *recvtype, *bloc, *bsize, *counts, 
	*send, *receive, *result, *sendBuf, *recvcounts;

	MPI_Comm_size(MPI_COMM_WORLD, &prosCount);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	result      = (int *) malloc((2*n)*sizeof(int)+1);
	sendBuf    	= (int *) malloc(prosCount*sizeof(int)+1);
	recvcounts 	= (int *) malloc((prosCount * prosCount)*sizeof(int)+1);
	bloc        = (int *) malloc(prosCount*sizeof(int)+1); 
	bsize       = (int *) malloc(prosCount*sizeof(int)+1); 
	counts      = (int *) malloc(prosCount*sizeof(int)+1);
	displs      = (int *) malloc(prosCount*sizeof(int)+1);
	receive     = (int *) malloc(prosCount*sizeof(int)+1);
	recvtype    = (int *) malloc(prosCount*sizeof(int)+1);
	scounts     = (int *) malloc(prosCount*sizeof(int)+1);
	sdispl      = (int *) malloc(prosCount*sizeof(int)+1);
	send        = (int *) malloc(prosCount*sizeof(int)+1);

	heapsort(data, n);

	for (i=0; i<prosCount; i++) {
		sendBuf[i] = data[i * n / prosCount];
		recvtype[i] = i * prosCount;
		displs[i] = prosCount;
	}

	MPI_Gatherv(sendBuf, prosCount, MPI_Type, recvcounts, displs,
	      recvtype, MPI_INT, 0, MPI_COMM_WORLD);

	if (myid == 0)
    {
		heapsort(recvcounts, prosCount * prosCount);
		for (i=0; i<prosCount; i++) 
		{
			sendBuf[i] = recvcounts[i * prosCount];
		}
    }

	MPI_Bcast(sendBuf, prosCount, MPI_Type, 0, MPI_COMM_WORLD);

	j=0;
	bloc[0] = 0;
	for (i=1; i<prosCount; i++)
    {
		while (IsSmaller(data[j],sendBuf[i])) {
			j++;
		}
		bloc[i]=j;
    }
  
	for (i=0; i<prosCount-1; i++)
    {
      bsize[i] = bloc[i+1]-bloc[i];
    }
	
	bsize[prosCount-1] = n - bloc[prosCount-1];

	for (i=0; i<prosCount; i++)
    {
      scounts[i] 	= 1;
      sdispl[i] 	= i;
      displs[i] 	= 1;
      recvtype[i] 	= i;
    }
    
	MPI_Alltoallv(bsize, scounts, sdispl, MPI_INT, counts, displs,
		recvtype,  MPI_Type, MPI_COMM_WORLD);

	resultSize = 1;
	for (i=0; i<prosCount; i++)
    {
      scounts[i] 	= bsize[i];
      sdispl[i] 	=  bloc[i];
      displs[i] 	= counts[i];
      recvtype[i] 	= resultSize-1; 
      resultSize 	+= counts[i];
    }
	resultSize--;

	MPI_Alltoallv(data, scounts, sdispl, MPI_Type, result, displs,
		recvtype,  MPI_Type, MPI_COMM_WORLD);

	heapsort(result, resultSize);

	for (i=0; i<prosCount; i++)
	{
		displs[i] = 1;
		recvtype[i]  = i;
	}

	MPI_Allgatherv(&resultSize, 1, MPI_INT, counts, displs, recvtype,
		 MPI_INT, MPI_COMM_WORLD);

	Left=0;
	for (i=0; i<myid; i++){
		Left += counts[i];
	}
	Right = Left + counts[myid] - 1;

	for (i=0; i<prosCount; i++)
	{
		l = i * n;
		r = ((i+1) * n) - 1;
		send[i] = 0;
      if ((Left <= l) && (l <= Right) && (Right <= r)) send[i] = Right-l+1;
      if ((l <= Left) && (Left <= r) && (r <= Right)) send[i] = r-Left+1;
      if ((l <= Left) && (Right <= r)) send[i] = Right-Left+1;
      if ((Left <= l) && (r <= Right)) send[i] = r-l+1;
	}

	for (i=0; i<prosCount; i++)
	{
		scounts[i] 	= 1;
		sdispl[i] 	= i;
		displs[i] 	= 1;
		recvtype[i] 	= i;
	}
	MPI_Alltoallv(send, scounts, sdispl, MPI_INT, receive, displs,
		recvtype,  MPI_Type, MPI_COMM_WORLD);

	l = 0;
	r = 0;
	for (i=0; i<prosCount; i++)
    {
		scounts[i] = send[i];
		sdispl[i] =  l; l += scounts[i];
		displs[i] = receive[i];
		recvtype[i] = r; r += displs[i];
    }

	MPI_Alltoallv(result, scounts, sdispl, MPI_Type, data, displs,
		recvtype,  MPI_Type, MPI_COMM_WORLD);
}

// sort into increasing order
// heapsort logic from Stack Overflow
// https://stackoverflow.com/questions/7520133/heap-sort-in-c
void heapsort(int data[], int n)
{
	int *b;
	b = data - 1;
	// put array into heap form
	for (int j = n/2; j > 0; j--)
		fix(data, j, n);
	for (int j = n-1; j > 0; j--)
	{
	    int temp = b[1];
	    b[1]=b[j+1];
	    b[j+1]=temp;
	    fix(data, 1, j);
	}
}

static void fix(int data[], int m, int n)
{
	int j, k;
	int *b;
	//origin is 1
	b = data - 1;  
	j = m;
	k = m * 2;
	while (k <= n)
	{
		if (k < n && IsSmaller(b[k],b[k+1])) ++k;
		if (IsSmaller(b[j],b[k]))
		{
			int temp = b[j];
			b[j]=b[k];
			b[k]=temp;
		}
		j = k;
		k *= 2;
	}
}
