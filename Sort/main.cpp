#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include "mpi.h"

#define InputFilePrefix "./inputFiles/"      //input file location 
#define OutputFilePrefix "./outputFiles/"    //output file location
#define ArrayData int
#define MPI_Type  MPI_INT
#define Smaller(x, y) ((x) < (y))
#define SWAP(x, y) (temp = (x), (x) = (y), (y) = temp)

ArrayData *MakeArray(int n)
{
	return((ArrayData *) malloc(n*sizeof(ArrayData)+1));
};

static void fix(ArrayData v[], register int m, register int n);
void heapsort(ArrayData v[], int n);
void sort(ArrayData data[], int localSize);

int main(int argc,char *argv[])
{
	int myid, prosCount;
	FILE *infile, *outfile;
	char infilepostfix[7], outfilepostfix[8];
	char infilename[100] = InputFilePrefix;
	char outfilename[100] = OutputFilePrefix;
	ArrayData *data;    // input data 
	int n;     		  	// input size 
	int p;     		// input proc 

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
	fprintf(outfile,"proc. %d: \n n = %d \n p = %d \n ", myid, n, p);
	 
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

static void fix(ArrayData v[], int m, int n)
{
	int j, k, temp;
	ArrayData *b;
	//origin is 1
	b = v - 1;  
	j = m;
	k = m * 2;
	while (k <= n)
	{
		if (k < n && Smaller(b[k],b[k+1])) ++k;
		if (Smaller(b[j],b[k]))  SWAP(b[j], b[k]);
		j = k;
		k *= 2;
	}
}

// sort into increasing order
void heapsort(ArrayData v[], int n)
{
	int j, temp;
	ArrayData *b;
	b = v - 1;
	// put array into heap form
	for (j = n/2; j > 0; j--)
		fix(v, j, n);
	for (j = n-1; j > 0; j--)
	{
	  SWAP(b[1], b[j+1]);
	  fix(v, 1, j);
	}
}

void sort(ArrayData data[], int n)
{
  int myid, prosCount;
  ArrayData *result;               
  int resultSize;               
  ArrayData *splitters, *allSplitters;
  int i, j, buf, tmp, l, r, L, R;
  int *scounts, *rcounts, *sdispl, 
      *rdispl, *bloc, *bsize, *counts, 
      *send, *receive;

  MPI_Comm_size(MPI_COMM_WORLD, &prosCount);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  result       = MakeArray((2 * n));
  splitters    = MakeArray(prosCount);
  allSplitters = MakeArray((prosCount * prosCount));
  scounts      = (int *) malloc(prosCount*sizeof(int)+1);
  rcounts      = (int *) malloc(prosCount*sizeof(int)+1);
  sdispl       = (int *) malloc(prosCount*sizeof(int)+1);
  rdispl       = (int *) malloc(prosCount*sizeof(int)+1);
  counts       = (int *) malloc(prosCount*sizeof(int)+1);
  send         = (int *) malloc(prosCount*sizeof(int)+1);
  receive      = (int *) malloc(prosCount*sizeof(int)+1);
  bloc         = (int *) malloc(prosCount*sizeof(int)+1); 
  bsize        = (int *) malloc(prosCount*sizeof(int)+1); 

  /* ----------------------- compute splitters ----------*/
  heapsort(data, n);

  for (i=0; i<prosCount; i++) {splitters[i] = data[i * n / prosCount];};
  for (i=0; i<prosCount; i++) {rdispl[i] = i * prosCount;};
  for (i=0; i<prosCount; i++) {rcounts[i] = prosCount;};

  MPI_Gatherv(splitters, prosCount, MPI_Type, allSplitters, rcounts,
	      rdispl, MPI_INT, 0, MPI_COMM_WORLD);

  if (myid == 0)
    {
      heapsort(allSplitters, prosCount * prosCount);
      for (i=0; i<prosCount; i++) {splitters[i] = allSplitters[i * prosCount];};
    };

  MPI_Bcast(splitters, prosCount, MPI_Type, 0, MPI_COMM_WORLD);

  /* ---- locate splitters in sorted local data ----------*/
  j=0;
  bloc[0] = 0;
  for (i=1; i<prosCount; i++)
    {
      while (Smaller(data[j],splitters[i])) {j++;};
      bloc[i]=j;
    }
  

  /* ---- compute #el. in each sub-bucket ----------*/
  for (i=0; i<prosCount-1; i++)
    {
      bsize[i] = bloc[i+1]-bloc[i];
    }
  bsize[prosCount-1] = n - bloc[prosCount-1];

  /* ----- communicate bucket sizes ---------------*/
  for (i=0; i<prosCount; i++)
    {
      scounts[i] = 1;
      sdispl[i] = i;
      rcounts[i] = 1;
      rdispl[i] = i;
    }
  MPI_Alltoallv(bsize, scounts, sdispl, MPI_INT, counts, rcounts,
		rdispl,  MPI_Type, MPI_COMM_WORLD);

  /* ----- send buckets to destination proc. ---------------*/
  resultSize = 1;
  for (i=0; i<prosCount; i++)
    {
      scounts[i] = bsize[i];
      sdispl[i] =  bloc[i];
      rcounts[i] = counts[i];
      rdispl[i] = resultSize-1; resultSize += counts[i];
    };
  resultSize--;

  MPI_Alltoallv(data, scounts, sdispl, MPI_Type, result, rcounts,
		rdispl,  MPI_Type, MPI_COMM_WORLD);

  /* ----- sort data at dest.proc. ---------------*/
  heapsort(result, resultSize);


  /* ----- balance data  ------------------------------------------*/

  /* compute array 'counts' with number of data items at each proc.*/
  for (i=0; i<prosCount; i++)
  {
      rcounts[i] = 1;
      rdispl[i] = i;
  };

  MPI_Allgatherv(&resultSize, 1, MPI_INT, counts, rcounts, rdispl,
		 MPI_INT, MPI_COMM_WORLD);

  /* compute array 'send' with number data items to be sent to each proc.*/
  L=0;
  for (i=0; i<myid; i++){L += counts[i];};
  R = L + counts[myid] - 1;

  for (i=0; i<prosCount; i++)
  {
      l = i * n;
      r = ((i+1) * n) - 1;
      send[i] = 0;
      /*-- L l R r ---*/
      if ((L <= l) && (l <= R) && (R <= r)) send[i] = R-l+1;
      /*-- l L r R ---*/
      if ((l <= L) && (L <= r) && (r <= R)) send[i] = r-L+1;
      /*-- l L R r ---*/
      if ((l <= L) && (R <= r)) send[i] = R-L+1;
      /*-- L l r R ---*/
      if ((L <= l) && (r <= R)) send[i] = r-l+1;
  };

  /* compute array 'receive' with number data items to be received each proc.*/
  for (i=0; i<prosCount; i++)
  {
      scounts[i] = 1;
      sdispl[i] = i;
      rcounts[i] = 1;
      rdispl[i] = i;
  }
  MPI_Alltoallv(send, scounts, sdispl, MPI_INT, receive, rcounts,
		rdispl,  MPI_Type, MPI_COMM_WORLD);

  /* move data */

  l = 0;
  r = 0;
  for (i=0; i<prosCount; i++)
    {
      scounts[i] = send[i];
      sdispl[i] =  l; l += scounts[i];
      rcounts[i] = receive[i];
      rdispl[i] = r; r += rcounts[i];
    };

  MPI_Alltoallv(result, scounts, sdispl, MPI_Type, data, rcounts,
		rdispl,  MPI_Type, MPI_COMM_WORLD);

  /* ------ free up memory space ---------------------------------*/
  free(result);
  free(splitters);
  free(allSplitters);
  free(scounts);
  free(rcounts);
  free(sdispl);
  free(rdispl);
  free(counts);
  free(send);
  free(receive);
  free(bloc);
  free(bsize);

};
