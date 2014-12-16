#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include "mpi.h"

#define InputFilePrefix "./inputFiles/"      //input file location 
#define OutputFilePrefix "./outputFiles/"    //output file location
#define MPI_Type  MPI_INT
#define Smaller(x, y) ((x) < (y))

int *MakeArray(int n)
{
	return((int *) malloc(n*sizeof(int)+1));
}

static void fix(int v[], register int m, register int n);
void heapsort(int v[], int n);
void sort(int data[], int localSize);

int main(int argc,char *argv[])
{
	int myid, prosCount;
	FILE *infile, *outfile;
	char infilepostfix[7], outfilepostfix[8];
	char infilename[100] = InputFilePrefix;
	char outfilename[100] = OutputFilePrefix;
	int *data;    // input data 
	int n;     		  	// input size 
	int p;     			// input proc 

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
	infilepostfix[6] = ((myid % 100) / 10) + 48;
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
	data = MakeArray(n);
	 
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
	printf("Proc. %d: Total Time = %d sec.\n", myid, (startTotal - stopTotal));
	MPI_Finalize();
}

// sort into increasing order
// heapsort logic from Stack Overflow
// https://stackoverflow.com/questions/7520133/heap-sort-in-c
void heapsort(int v[], int n)
{
	int j, temp;
	int *b;
	b = v - 1;
	// put array into heap form
	for (j = n/2; j > 0; j--)
		fix(v, j, n);
	for (j = n-1; j > 0; j--)
	{
	  temp = b[1];
	  b[1]=b[j+1];
	  b[j+1]=temp;
	  fix(v, 1, j);
	}
}

static void fix(int v[], int m, int n)
{
	int j, k, temp;
	int *b;
	//origin is 1
	b = v - 1;  
	j = m;
	k = m * 2;
	while (k <= n)
	{
		if (k < n && Smaller(b[k],b[k+1])) ++k;
		if (Smaller(b[j],b[k]))
		{
			temp = b[j];
			b[j]=b[k];
			b[k]=temp;
		}
		j = k;
		k *= 2;
	}
}

void sort(int data[], int n)
{
	int myid, prosCount, resultSize,i, j, buf, tmp, l, r, L, R;
	int *scounts, *rcounts, *sdispl, *rdispl, *bloc, *bsize, *counts, 
	*send, *receive, *result, *splits, *allSplits;

	MPI_Comm_size(MPI_COMM_WORLD, &prosCount);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	result      = MakeArray((2 * n));
	splits    	= MakeArray(prosCount);
	allSplits 	= MakeArray((prosCount * prosCount));
	scounts     = (int *) malloc(prosCount*sizeof(int)+1);
	rcounts     = (int *) malloc(prosCount*sizeof(int)+1);
	sdispl      = (int *) malloc(prosCount*sizeof(int)+1);
	rdispl      = (int *) malloc(prosCount*sizeof(int)+1);
	counts      = (int *) malloc(prosCount*sizeof(int)+1);
	send        = (int *) malloc(prosCount*sizeof(int)+1);
	receive     = (int *) malloc(prosCount*sizeof(int)+1);
	bloc        = (int *) malloc(prosCount*sizeof(int)+1); 
	bsize       = (int *) malloc(prosCount*sizeof(int)+1); 

	heapsort(data, n);

	for (i=0; i<prosCount; i++) {
		splits[i] = data[i * n / prosCount];
		rdispl[i] = i * prosCount;
		rcounts[i] = prosCount;
	}

	MPI_Gatherv(splits, prosCount, MPI_Type, allSplits, rcounts,
	      rdispl, MPI_INT, 0, MPI_COMM_WORLD);

	if (myid == 0)
    {
		heapsort(allSplits, prosCount * prosCount);
		for (i=0; i<prosCount; i++) 
		{
			splits[i] = allSplits[i * prosCount];
		}
    }

	MPI_Bcast(splits, prosCount, MPI_Type, 0, MPI_COMM_WORLD);

	j=0;
	bloc[0] = 0;
	for (i=1; i<prosCount; i++)
    {
		while (Smaller(data[j],splits[i])) {
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
      rcounts[i] 	= 1;
      rdispl[i] 	= i;
    }
    
	MPI_Alltoallv(bsize, scounts, sdispl, MPI_INT, counts, rcounts,
		rdispl,  MPI_Type, MPI_COMM_WORLD);

	resultSize = 1;
	for (i=0; i<prosCount; i++)
    {
      scounts[i] 	= bsize[i];
      sdispl[i] 	=  bloc[i];
      rcounts[i] 	= counts[i];
      rdispl[i] 	= resultSize-1; 
      resultSize 	+= counts[i];
    }
	resultSize--;

	MPI_Alltoallv(data, scounts, sdispl, MPI_Type, result, rcounts,
		rdispl,  MPI_Type, MPI_COMM_WORLD);

	heapsort(result, resultSize);

	for (i=0; i<prosCount; i++)
	{
		rcounts[i] = 1;
		rdispl[i]  = i;
	}

	MPI_Allgatherv(&resultSize, 1, MPI_INT, counts, rcounts, rdispl,
		 MPI_INT, MPI_COMM_WORLD);

	L=0;
	for (i=0; i<myid; i++){
		L += counts[i];
	}
	R = L + counts[myid] - 1;

	for (i=0; i<prosCount; i++)
	{
		l = i * n;
		r = ((i+1) * n) - 1;
		send[i] = 0;
      if ((L <= l) && (l <= R) && (R <= r)) send[i] = R-l+1;
      if ((l <= L) && (L <= r) && (r <= R)) send[i] = r-L+1;
      if ((l <= L) && (R <= r)) send[i] = R-L+1;
      if ((L <= l) && (r <= R)) send[i] = r-l+1;
	}

	for (i=0; i<prosCount; i++)
	{
		scounts[i] 	= 1;
		sdispl[i] 	= i;
		rcounts[i] 	= 1;
		rdispl[i] 	= i;
	}
	MPI_Alltoallv(send, scounts, sdispl, MPI_INT, receive, rcounts,
		rdispl,  MPI_Type, MPI_COMM_WORLD);

	l = 0;
	r = 0;
	for (i=0; i<prosCount; i++)
    {
		scounts[i] = send[i];
		sdispl[i] =  l; l += scounts[i];
		rcounts[i] = receive[i];
		rdispl[i] = r; r += rcounts[i];
    }

	MPI_Alltoallv(result, scounts, sdispl, MPI_Type, data, rcounts,
		rdispl,  MPI_Type, MPI_COMM_WORLD);

//free all the tihngs
	free(result);
	free(splits);
	free(allSplits);
	free(scounts);
	free(rcounts);
	free(sdispl);
	free(rdispl);
	free(counts);
	free(send);
	free(receive);
	free(bloc);
	free(bsize);
}
