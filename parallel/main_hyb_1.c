/*TODO: Add a header  */
#define _POSIX_C_SOURCE 199309L
#include<stdio.h>
#include<stdlib.h>
#include<errno.h>
#include<string.h>
#include<sys/time.h>
#include<time.h>
#include<math.h>
#include"inc/mmio.h"
#include"mpi.h"
#include"omp.h"
extern int errno ;
/* Structs Declaration */
typedef struct cscMatrix{
	int *csc_r;
	int *csc_c;
	int nnz;
	int rows;
	int i;
	int j;
}cscMat;
typedef struct cooMatrix{
	int *coo_r;
	int *coo_c;
	int nnz;
	int rows;
	int i;
	int j;
}cooMat;
/* Function Predeclarations */
void checkArgsNum(int argc);
void coo2csc(int * const row,int * const col,int const * const row_coo,
		int const * const col_coo,int const nnz,int const n,int const isOneBased);
cscMat readFile(char *filename);
cooMat sortMat(cooMat Matrix);
void matrixMult(cscMat MatA, cscMat MatB, int **totalHits, int *totalHitsSize, int firstFlag);
int imin(int a, int b);

/* <<< --- MAIN --- >>> */
int
main(int argc, char *argv[]){
	/* Check arguments if correct (Must be 2 names of .mtx files) */
	checkArgsNum(argc);
	/* Get MatA */
	cscMat MatA;
	MatA = readFile(argv[1]);
	/* Get MatB */
	cscMat MatB;
	MatB = readFile(argv[2]);
	/* Init MatC*/
	cooMat MatC;
	MatC.rows=MatA.rows;
	MatC.nnz=0;
	MatC.coo_r=NULL;
	MatC.coo_c=NULL;
	/* Timing variables< */
	struct timespec ts_start;
	struct timespec ts_end;
	clock_gettime(CLOCK_MONOTONIC,&ts_start);


	/* MPI */
	int numtasks,rank;
	struct cooMatrix *MatrixCOOArrA, *MatrixCOOArrB, MatrixCOOArrC;
	struct cscMatrix *privMatrixArrA, *MatrixArrB;

	MPI_Init(NULL,NULL);
	MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	//printf("tasks %d rank %d\n",numtasks,rank);

	/* Block size for each axis */
	int block = numtasks;
	/* Check the size of the block if is greater that rows */
	if(MatC.rows/block<2){
		printf("Block too large Block Size %d Rows %d\n Exiting\n",block,MatC.rows);
		exit(EXIT_FAILURE);
	}
	int blockalloc=block; // if the rows%block is larger put an extra bucket
	if (MatC.rows%block!=0 ) blockalloc++;

	/* Assign the workload to every thread */
	/* Every process gets a block of columns to calculate from the matrix C */
	int workStart,workEnd;
	/* Check if the blocks allocated are equally divaded with num of processes */
	if (blockalloc%numtasks==0){
		int workload = (int) blockalloc/numtasks;
		workStart= rank * workload;
		workEnd = workStart + workload;
	}
	/* If not then balance the work of each process 
	   i.e num proc = 3, blockalloc=10 p1|4 p2|3 p3|3 total 10*/
	else {
		int workload = floor(blockalloc/numtasks);
		if (rank < blockalloc - workload*numtasks){
			workStart = (workload+1) *rank;
			workEnd = workStart+workload+1;
		}
		else {
			int mod = blockalloc - workload*numtasks;
			workStart = mod + workload*mod + (rank - mod)*workload;
			workEnd  =workStart +workload;
		}
	}
	int totalWorkload=workEnd-workStart;
	printf("rank %d start %d end %d\n",rank,workStart,workEnd);

	/* Allocate an arrays that will contain the  blocks COO and CSC  */
	MatrixCOOArrA = malloc(blockalloc*totalWorkload*sizeof(struct cooMatrix));
	MatrixCOOArrB = malloc(blockalloc*totalWorkload*sizeof(struct cooMatrix));
	if (MatrixCOOArrA==NULL || MatrixCOOArrB==NULL) {
		printf("Memory allocation at line %d Failed\n",__LINE__);
		exit(EXIT_FAILURE);
	}
	for (int i=0;i<blockalloc*totalWorkload;i++){
		MatrixCOOArrA[i].nnz=0;
		MatrixCOOArrA[i].coo_r=NULL;
		MatrixCOOArrA[i].coo_c=NULL;
		MatrixCOOArrA[i].rows=(int) floor(MatA.rows/block);
		MatrixCOOArrB[i].nnz=0;
		MatrixCOOArrB[i].coo_r=NULL;
		MatrixCOOArrB[i].coo_c=NULL;
		MatrixCOOArrB[i].rows=(int) floor(MatB.rows/block);
	}

	/* Disect each Matrix */
	// A
	for(int j=workStart*MatrixCOOArrA[0].rows;j<imin(workEnd*MatrixCOOArrA[0].rows,MatA.rows);j++){
		for(int r=MatA.csc_c[j];r<MatA.csc_c[j+1];r++){
			int i = MatA.csc_r[r];
			/* Find Arri and Arrj  */
			int Arri = (int) floor(i/MatrixCOOArrA[0].rows);
			int Arrj = (int) floor(j/MatrixCOOArrA[0].rows);
			//printf("i %d j %d Arri %d Arrj %d \n",i,j,Arri,Arrj);

			/* Check if structure exist in MatrixCOOArrA and allocate it */	  
			int ptr = Arri +blockalloc*(Arrj-workStart);
			MatrixCOOArrA[ptr].nnz ++;
			MatrixCOOArrA[ptr].coo_r = (int *)realloc(MatrixCOOArrA[ptr].coo_r,
					MatrixCOOArrA[ptr].nnz*sizeof(int));
			MatrixCOOArrA[ptr].coo_c = (int *)realloc(MatrixCOOArrA[ptr].coo_c,
					MatrixCOOArrA[ptr].nnz*sizeof(int));
			if (MatrixCOOArrA[ptr].coo_c == NULL || MatrixCOOArrA[ptr].coo_r==NULL)
			{fprintf(stderr,"Line %d: Error Alocating MatrixCOOArrA",__LINE__);
				exit(EXIT_FAILURE);}
			/* The final i and j that re pushed to the coo array are shifted to their relative
			   positions in respect to the block coordinates */
			MatrixCOOArrA[ptr].coo_r[MatrixCOOArrA[ptr].nnz-1] = i - Arri*MatrixCOOArrA[ptr].rows;
			MatrixCOOArrA[ptr].coo_c[MatrixCOOArrA[ptr].nnz-1] = j - Arrj*MatrixCOOArrA[ptr].rows;;
		}
	}
	// B
	for(int j=workStart*MatrixCOOArrA[0].rows;j<imin(workEnd*MatrixCOOArrA[0].rows,MatB.rows);j++){
		for(int r=MatB.csc_c[j];r<MatB.csc_c[j+1];r++){
			int i = MatB.csc_r[r];
			/* Find Arri and Arrj  */
			int Arri = (int) floor(i/MatrixCOOArrB[0].rows);
			int Arrj = (int) floor(j/MatrixCOOArrB[0].rows);
			//printf("Arri %d Arrj %d for i %d j %d \n",Arri,Arrj,i,j);

			/* Check if structure exist in MatrixCOOArrB and allocate it */	  
			int ptr = Arri +blockalloc*(Arrj-workStart);
			MatrixCOOArrB[ptr].nnz ++;
			MatrixCOOArrB[ptr].coo_r = (int *)realloc(MatrixCOOArrB[ptr].coo_r,MatrixCOOArrB[ptr].nnz*sizeof(int));
			MatrixCOOArrB[ptr].coo_c = (int *)realloc(MatrixCOOArrB[ptr].coo_c,MatrixCOOArrB[ptr].nnz*sizeof(int));
			if (MatrixCOOArrB[ptr].coo_c == NULL || MatrixCOOArrB[ptr].coo_r==NULL)
			{fprintf(stderr,"Line %d: Error Alocating MatrixCOOArrB",__LINE__);
				exit(EXIT_FAILURE);}
			/* The final i and j that re pushed to the coo array are shifted to their relative
			   positions in respect to the block coordinates */
			MatrixCOOArrB[ptr].coo_r[MatrixCOOArrB[ptr].nnz-1] = i - Arri*MatrixCOOArrB[ptr].rows;
			MatrixCOOArrB[ptr].coo_c[MatrixCOOArrB[ptr].nnz-1] = j - Arrj*MatrixCOOArrB[ptr].rows;
		}
	}
	/* Covert All COO to CSC */
	//TODO: Remember to free
	privMatrixArrA =  malloc(blockalloc*totalWorkload*sizeof(struct cscMatrix));
	MatrixArrB =  malloc(blockalloc*totalWorkload*sizeof(struct cscMatrix));
	if (privMatrixArrA==NULL || MatrixArrB==NULL) {
		printf("Memory allocation at line %d Failed\n",__LINE__);
		exit(EXIT_FAILURE);
	}
	#pragma omp parallel for schedule(dynamic,1) nowait
	for(int i=0;i<blockalloc;i++){
		for(int j=0;j<totalWorkload;j++){
			int ptr = i + blockalloc*j;
			/* Declaration of CSC format arrays  */ 
			int row = MatrixCOOArrA[ptr].rows;
			int nnz = MatrixCOOArrA[ptr].nnz;
			int isOneBased =0; 
			/*Allocate memory for csc  matrices */
			privMatrixArrA[ptr].nnz = nnz;
			privMatrixArrA[ptr].rows =row;
			privMatrixArrA[ptr].i=i*row;
			privMatrixArrA[ptr].j=(j+workStart)*row;
			privMatrixArrA[ptr].csc_r = (int *)malloc(nnz     * sizeof(int));
			privMatrixArrA[ptr].csc_c = (int *)malloc((row+1) * sizeof(int));
			if (privMatrixArrA[ptr].csc_r==NULL || privMatrixArrA[ptr].csc_c==NULL) {
				printf("Memory allocation at line %d Failed\n",__LINE__);
				exit(EXIT_FAILURE);
			}
			/* Conversion */	
			if(MatrixCOOArrA[ptr].nnz==0){	// If the array is empty nnz=0 then the arrays are NULL
				privMatrixArrA[ptr].csc_r = NULL;
				privMatrixArrA[ptr].csc_c = NULL;
			}
			else{
				coo2csc(privMatrixArrA[ptr].csc_r,privMatrixArrA[ptr].csc_c,
						MatrixCOOArrA[ptr].coo_r,MatrixCOOArrA[ptr].coo_c,
						nnz,row,isOneBased);
			}
			/* Free Memory */
			free(MatrixCOOArrA[ptr].coo_r);
			free(MatrixCOOArrA[ptr].coo_c);
			// B
			/* Declaration of CSC format arrays */ 
			row = MatrixCOOArrB[ptr].rows;
			nnz = MatrixCOOArrB[ptr].nnz;
			isOneBased =0; 
			/*Allocate memory for csc  matrices */
			MatrixArrB[ptr].csc_r = (int *)malloc(nnz     * sizeof(int));
			MatrixArrB[ptr].csc_c = (int *)malloc((row+1) * sizeof(int));
			MatrixArrB[ptr].nnz = nnz;
			MatrixArrB[ptr].rows =row;
			MatrixArrB[ptr].i=i*row;
			MatrixArrB[ptr].j=(j+workStart)*row;
			if (MatrixArrB[ptr].csc_r==NULL || MatrixArrB[ptr].csc_c==NULL) {
				printf("Memory allocation at line %d Failed\n",__LINE__);
				exit(EXIT_FAILURE);
			}
			/* Conversion */	
			if(MatrixCOOArrB[ptr].nnz==0){	// If the array is empty nnz=0 then the arrays are NULL
				MatrixArrB[ptr].csc_r = NULL;
				MatrixArrB[ptr].csc_c = NULL;
			}
			else{
				coo2csc(MatrixArrB[ptr].csc_r,MatrixArrB[ptr].csc_c,
						MatrixCOOArrB[ptr].coo_r,MatrixCOOArrB[ptr].coo_c,
						nnz,row,isOneBased);
			}
			/* Free Memory */
			free(MatrixCOOArrB[ptr].coo_r);
			free(MatrixCOOArrB[ptr].coo_c);

		}
	}
	/* Free Memory */
	free(MatrixCOOArrA);
	free(MatrixCOOArrB);
	/* Allocate memory fore the common MatrixArrA */
	struct cscMatrix *MatrixArrA =  malloc(blockalloc*blockalloc*sizeof(struct cscMatrix));
	/* rankMat is an array containing the workEnd value of every proccess in order to find root */
	int *rankMat = (int *)malloc(numtasks*sizeof(int));
       	MPI_Allgather(&workEnd,1,MPI_INT,rankMat,1,MPI_INT,MPI_COMM_WORLD);
	/* Gather the Results */
	for(int i=0;i<blockalloc;i++){
		for(int j=0;j<blockalloc;j++){
			int ptr = i+blockalloc*j;
			int root;
			/* Sender code */
			if(j>=workStart && j<workEnd){
				int privptr = i+blockalloc*(j-workStart);
				root=rank;
				/* WARNING: Do not free anything yet */
				/* put the values in place */
				MatrixArrA[ptr].csc_r = privMatrixArrA[privptr].csc_r; 
				MatrixArrA[ptr].csc_c = privMatrixArrA[privptr].csc_c; 
				MatrixArrA[ptr].nnz = privMatrixArrA[privptr].nnz; 
				MatrixArrA[ptr].rows = privMatrixArrA[privptr].rows;
				MatrixArrA[ptr].i = privMatrixArrA[privptr].i;
				MatrixArrA[ptr].j = privMatrixArrA[privptr].j;

			} 
			else{
				/* Check the rankMat array for the root */
				for(int r=0;r<numtasks;r++){
					if(j<rankMat[r]) {
						root =r; 
						break;
					}
				}
			
			}
			/* broadcast values */
			MPI_Bcast(&MatrixArrA[ptr].nnz,1,MPI_INT,root,MPI_COMM_WORLD);
			MPI_Bcast(&MatrixArrA[ptr].rows,1,MPI_INT,root,MPI_COMM_WORLD);
			MPI_Bcast(&MatrixArrA[ptr].i,1,MPI_INT,root,MPI_COMM_WORLD);
			MPI_Bcast(&MatrixArrA[ptr].j,1,MPI_INT,root,MPI_COMM_WORLD);
			/* allocate memory for the csc_c and csc_r if the nnz !=0*/
			if(rank!=root && MatrixArrA[ptr].nnz!=0){
				MatrixArrA[ptr].csc_r = (int *)malloc(MatrixArrA[ptr].nnz*sizeof(int));
				MatrixArrA[ptr].csc_c = (int *)malloc((MatrixArrA[ptr].rows+1)*sizeof(int));
			}
			/* Check if nnz are 0 and pass the NULL values */
			if (MatrixArrA[ptr].nnz==0){
				MatrixArrA[ptr].csc_r =NULL;	
				MatrixArrA[ptr].csc_c =NULL;	
			}
			else{
				/* If the nnz are not 0 then broadcast the arrays */
			MPI_Bcast(MatrixArrA[ptr].csc_r,MatrixArrA[ptr].nnz,MPI_INT,root,MPI_COMM_WORLD);
			MPI_Bcast(MatrixArrA[ptr].csc_c,MatrixArrA[ptr].rows+1,MPI_INT,root,MPI_COMM_WORLD);
			}
		}
	}

	/* DEBUG Print the buckets */
	/*
	   if (rank==0){
	   for(int i=0;i<blockalloc;i++){
	   for(int j=0;j<blockalloc;j++){
	   int ptr = i + blockalloc*j;
	   printf("i %d j %d arri %d arrj%d \n",i,j,MatrixArrA[ptr].i,MatrixArrA[ptr].j);
	   if(MatrixArrA[ptr].csc_c==NULL) {printf("NULL\n");continue;}
	   for (int c=0;c<MatrixArrA[ptr].rows;c++){
	   for(int ci=MatrixArrA[ptr].csc_c[c];ci<MatrixArrA[ptr].csc_c[c+1];ci++){
	   int r = MatrixArrA[ptr].csc_r[ci];
	   printf("%d %d \n",r+MatrixArrA[ptr].i,c+MatrixArrA[ptr].j);
	   }
	   }
	   puts("\n");
	   }
	   }
	   }
	   exit(EXIT_SUCCESS) ;
	   */
	MatrixCOOArrC.nnz=0;
	MatrixCOOArrC.coo_c=NULL;
	MatrixCOOArrC.coo_r=NULL;
	/* printf("Rank %d Start %d End %d\n",rank,workStart,workEnd); */
	/* Start Calculating the result */
	/* Assign columns of C to processes */
	for (int Ccol=workStart;Ccol<workEnd;Ccol++){
		/* Iterate the rows of C */
		for (int Crow=0;Crow<blockalloc;Crow++){
			/* This is the temp matrix that contain the total elements for a value*/
			int Cptr = Crow+blockalloc*Ccol;
			cooMat tempSubC;
			tempSubC.coo_r = NULL;
			tempSubC.coo_c = NULL;
			tempSubC.nnz = 0;
			tempSubC.i = Crow * MatrixArrA[Cptr].rows;
			tempSubC.j = Ccol * MatrixArrA[Cptr].rows;
			/* Allocate an array that contains the hits for each column */
			int **totalHits = (int **) malloc(MatrixArrA[Cptr].rows*sizeof(int*));
			int *totalHitsSize = (int *) malloc(MatrixArrA[Cptr].rows*sizeof(int));
			if (totalHits==NULL || totalHitsSize==NULL){
				printf("Memory allocation at line %d Failed\n",__LINE__);
				exit(EXIT_FAILURE);
			}
			int firstFlag = 1;
			/* if (rank==0) printf("Arri %d Arrj %d rows %d\n",Crow,Ccol,MatrixArrA[Cptr].rows); */
			/* Iterate the blocks in respect to row and column */
			for(int Ck=0;Ck<blockalloc;Ck++){
				int Aptr = Crow+blockalloc*Ck;
				int Bptr = Ck+blockalloc*(Ccol-workStart);
				// If block is empty skip
				if (MatrixArrA[Aptr].nnz>0 && MatrixArrB[Bptr].nnz>0){
					/* Multiply the arrays A x B = Ctemp */
					matrixMult(MatrixArrA[Aptr],MatrixArrB[Bptr],totalHits,totalHitsSize,firstFlag);
					firstFlag=0;
					// DEBUG
					/* if (rank ==0){ */
					/* 	printf("nnz %d\n",tempC.nnz); */
					/* for (int i=0;i<tempC.nnz;i++){ */
					/* 	if(rank==0) printf("Arri %d Arrj %d i %d j %d\n",Crow,Ccol,tempC.coo_r[i]+MatrixArrA[Cptr].i+1,tempC.coo_c[i]+MatrixArrA[Cptr].j+1); */
					/* } */
					/* } */
					/* Merge the results of temporary matrix to MatrixCOOArrC*/
				} 
			}
			/* The Last time create the MatrixCOOArrC */
			/* Allocate memory */
			if(firstFlag!=1){
				for(int j=0;j<MatrixArrA[Cptr].rows;j++){
					int nnz= 0;
					for (int nnz=0;nnz<totalHitsSize[j];nnz++){
						tempSubC.nnz ++;
						tempSubC.coo_r = (int *)realloc(tempSubC.coo_r,tempSubC.nnz*sizeof(int));
						tempSubC.coo_c = (int *)realloc(tempSubC.coo_c,tempSubC.nnz*sizeof(int));
						if (tempSubC.coo_c == NULL || tempSubC.coo_r==NULL)
						{fprintf(stderr,"Line %d: Error Alocating Matrix C",__LINE__);
							exit(EXIT_FAILURE);}
						tempSubC.coo_c[tempSubC.nnz-1] = j;
						tempSubC.coo_r[tempSubC.nnz-1] = totalHits[j][nnz];
					}
					free(totalHits[j]);
				}
			}
			/* Free arrays allocated for omp */
			free(totalHits);
			free(totalHitsSize);
			/* Realloc space for the new elements */
			int startIndex = MatrixCOOArrC.nnz;
			MatrixCOOArrC.nnz+=tempSubC.nnz;
			MatrixCOOArrC.coo_r = (int *)realloc(MatrixCOOArrC.coo_r,MatrixCOOArrC.nnz*sizeof(int));
			MatrixCOOArrC.coo_c = (int *)realloc(MatrixCOOArrC.coo_c,MatrixCOOArrC.nnz*sizeof(int));
			/* Merge the current tempSubC*/
			for (int i=0 ;i<tempSubC.nnz;i++){
				int ptr = startIndex +i;
				MatrixCOOArrC.coo_c[ptr] = tempSubC.coo_c[i] + tempSubC.j;
				MatrixCOOArrC.coo_r[ptr] = tempSubC.coo_r[i] + tempSubC.i;
			}
			/* tempSubC */
			free(tempSubC.coo_r);
			free(tempSubC.coo_c);
		}
	}
	/* if(rank==2){ */
	/* 	puts("\nMATRIX C\n"); */
	/* 	printf("nnz %d\n",MatrixCOOArrC.nnz); */
	/* 	for (int i =0 ;i<MatrixCOOArrC.nnz;i++) printf("%d %d\n",MatrixCOOArrC.coo_r[i]+1,MatrixCOOArrC.coo_c[i]+1); */
	/*  }  */
	/* sortMat(MatrixCOOArrC); */
	/* Get the displacement to fit the arrays */
	int recvcount[numtasks], displs[numtasks];
	MPI_Gather(&MatrixCOOArrC.nnz,1,MPI_INT,recvcount,1,MPI_INT,0,MPI_COMM_WORLD);
	/* MPI Gather the array */
	if(rank==0){
		displs[0]=0;
		for (int i=0;i<numtasks;i++)	MatC.nnz += recvcount[i];
		for (int i=0;i<numtasks-1;i++) 	displs[i+1]=displs[i]+recvcount[i];	
		MatC.coo_r = (int *) malloc(MatC.nnz*sizeof(int));
		MatC.coo_c = (int *) malloc(MatC.nnz*sizeof(int));
	}
	MPI_Gatherv(MatrixCOOArrC.coo_r,MatrixCOOArrC.nnz,MPI_INT,
		    MatC.coo_r,recvcount,displs,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Gatherv(MatrixCOOArrC.coo_c,MatrixCOOArrC.nnz,MPI_INT,
		    MatC.coo_c,recvcount,displs,MPI_INT,0,MPI_COMM_WORLD);
	/* MPI END */ 
	MPI_Finalize();
	clock_gettime(CLOCK_MONOTONIC,&ts_end);
	double ts_sec,ts_nsec;
	ts_sec = ts_end.tv_sec - ts_start.tv_sec;
	ts_nsec = ts_end.tv_nsec - ts_start.tv_nsec;
	printf("Execution Time %f ms\n",ts_sec*1000+ts_nsec/1000000);
	/* Sort MatC by column for tests  */
	if(rank==0){
		/* MatC = sortMat(MatC); */
		/* for (int i=0 ; i<MatC.nnz;i++) printf("Num %d i %d j %d \n",i+1,MatC.coo_r[i]+1,MatC.coo_c[i]+1); */
		printf("NNZ %d\n",MatC.nnz);
	}

	/* DEBUGGING PRINTOUTS */
	/* 
	   puts("MATRIX A\n");
	   for (int i =0 ;i<MatA.nnz;i++) printf("%d ",MatA.csc_r[i]);
	   puts("\n");
	   for (int i =0 ;i<MatA.rows+1;i++) printf("%d ",MatA.csc_c[i]);
	   puts("\nMATRIX B\n");
	   for (int i =0 ;i<MatB.nnz;i++) printf("%d ",MatB.csc_r[i]);
	   puts("\n");
	   for (int i =0 ;i<MatB.rows+1;i++) printf("%d ",MatB.csc_c[i]);
	   puts("\n");

*/

	  /* puts("\nMATRIX C\n"); */
	  /* for (int i =0 ;i<MatC.nnz;i++) printf("%d %d\n",MatC.coo_r[i],MatC.coo_c[i]); */


	/* Deallocate Matrices */
	for(int i=0;i<blockalloc*blockalloc;i++){
		free(MatrixArrA[i].csc_r);
		free(MatrixArrA[i].csc_c);
	}
	free(MatrixArrA);
	free(privMatrixArrA);
	free(MatrixArrB);
	free(MatrixCOOArrC.coo_c);
	free(MatrixCOOArrC.coo_r);
	/* Free csc arrays allocated from readFile(...)  */
	free(MatA.csc_r);
	free(MatA.csc_c);
	free(MatB.csc_r);
	free(MatB.csc_c);
	free(MatC.coo_r);
	free(MatC.coo_c);
}
/* Function that returns the min array */
int
imin(int a,int b){
	if(a>=b) {return b;}
	else {return a;}
}
/* This routine Multiplies 2 csc matrixes */
void 
matrixMult(cscMat MatA, cscMat MatB, int **totalHits, int *totalHitsSize, int firstFlag){
	/* This currently works for matrices RowsxRows */
	#pragma omp parallel 
	{
	#pragma omp parallel for schecdule(dynamic,256) nowait 
	for(int j=0;j<MatA.rows;j++){
		/* Check if is the first pass and the totalHit is empty */
		if(firstFlag==1){
			/* set finalHits to NULL */
			totalHits[j]=NULL;
			totalHitsSize[j]=0;
		}
		/* Iterate B column */
		for (int c=MatB.csc_c[j];c<MatB.csc_c[j+1];c++){
			int MatArow=MatB.csc_r[c];
			for (int r=MatA.csc_c[MatArow];r<MatA.csc_c[MatArow+1];r++){
				int i = MatA.csc_r[r];
				/* Set a flag that tracks if a same value (hit) is the answers */
				int hitflag=0;
				/* Check in current hits if this one exists */
				for(int hit=0;hit<totalHitsSize[j];hit++){
					if(i==totalHits[j][hit]) {
						/* set the flag and break */
						hitflag=1;
						break;
					}
				}
				if(hitflag==0){
					totalHitsSize[j]++;
					totalHits[j] = (int *)realloc(totalHits[j],totalHitsSize[j]*sizeof(int));
					totalHits[j][totalHitsSize[j]-1] = i;
				}
			}
		}
	}
	}
}
/* This function sorts a COO Matrix sorted by row to COO Matrix sorted by column */
/* Insertion Sort */
cooMat
sortMat(cooMat Matrix){
	/*Check Element if is smaller than the previous */
	for (int i =1 ;i<Matrix.nnz;i++){
		/* If it's smaller rollback by checking each element and swap with the i element until sorted in place */
		if (Matrix.coo_c[i]<Matrix.coo_c[i-1]){
			for (int j =i-1;j>=0;j--){
				if (Matrix.coo_c[j]>Matrix.coo_c[j+1]){
					int temp;
					temp = Matrix.coo_c[j+1];
					Matrix.coo_c[j+1] = Matrix.coo_c[j];
					Matrix.coo_c[j] = temp;
					temp = Matrix.coo_r[j+1];
					Matrix.coo_r[j+1] = Matrix.coo_r[j];
					Matrix.coo_r[j] = temp;
				}
			}
		}
	}
	int end,start,offset;
	end =0;
	while(1){
		start = end;
		offset=1;
		while(1){
			if (Matrix.coo_c[start]!=Matrix.coo_c[start+offset] || start+offset == Matrix.nnz ){
				end=start+offset;
				break;
			}
			else offset++;
		}
		/*Check Element if is smaller than the previous */
		for (int i =start+1 ;i<end;i++){
			/* If it's smaller rollback by checking each element and swap with the i element until sorted in place */
			if (Matrix.coo_r[i]<Matrix.coo_r[i-1]){
				for (int j =i-1;j>=start;j--){
					if (Matrix.coo_r[j]>Matrix.coo_r[j+1]){
						int temp;
						temp = Matrix.coo_c[j+1];
						Matrix.coo_c[j+1] = Matrix.coo_c[j];
						Matrix.coo_c[j] = temp;
						temp = Matrix.coo_r[j+1];
						Matrix.coo_r[j+1] = Matrix.coo_r[j];
						Matrix.coo_r[j] = temp;
					}
				}
			}
		}
		/* Check if the end is the end of the array */
		if(end==Matrix.nnz) break;
	}
	return Matrix;
}

cscMat
readFile(char *filename){
	FILE *f;
	int M,N,nz,ret;
	MM_typecode matcode;
	/* Open File and Handle errors */
	f= fopen (filename,"r");
	if (f == NULL){
		fprintf(stderr, "Line %d: Error opening file %s\n",__LINE__,strerror(errno));
		exit (EXIT_FAILURE);
	}
	/* Read Banner*/
	if(mm_read_banner(f,&matcode) != 0){
		printf("Matrix banner read not succesful\n");
		exit(1);
	}
	/* If banner is wrong format (array,dense) exit */
	if ((mm_is_matrix(matcode)==0)||(mm_is_sparse(matcode)==0)){
		printf("Check the format from the banner\n");printf("Must be matrix, sparse \n ");
		printf("Banner\n %s",mm_typecode_to_str(matcode));
	}	

	/* Read the matrix sizes for coordinate format*/
	if ((ret=mm_read_mtx_crd_size(f,&M,&N,&nz))!=0){
		printf("Cannot read coordinate");
		exit(0);
	}
	/* Allocating Memory for the COO format */
	int *coo_rows, *coo_cols;
	coo_rows = (int *)malloc(nz * sizeof(int));
	coo_cols = (int *)malloc(nz * sizeof(int));
	if (coo_rows==NULL || coo_cols==NULL) {
		printf("Memory allocation at line %d Failed\n",__LINE__);
		exit(EXIT_FAILURE);
	}
	/* Iterate .mtx and scan the values  */
	for (int i=0; i<nz; i++){
		ret = fscanf(f,"%d %d\n",&coo_rows[i],&coo_cols[i]);
		coo_rows[i]--;
		coo_cols[i]--;
	}
	/* Close the file */ 
	fclose(f);
	/* DEBUGGING TEMPORARY PRINTOUTS */
	/*printf("\tCOO\t\nRows - Cols - Non Zero Elements\n");
	  printf("%d - %d - %d\n",M,N,nz);
	  for (int i=0; i<10;i++)printf("Row %d Col %d\n",coo_rows[i],coo_cols[i]);i*/
	/* Convert to CSC  */
	/* Declaration of CSC format arrays  */ 
	int const row = M;
	int const nnz = nz;
	int isOneBased =0; 
	/*Allocate memory for csc  matrices */
	int * csc_r = (int *)malloc(nnz     * sizeof(int));
	int * csc_c = (int *)malloc((row+1) * sizeof(int));
	if (csc_r==NULL || csc_c==NULL) {
		printf("Memory allocation at line %d Failed\n",__LINE__);
		exit(EXIT_FAILURE);
	}
	/* Conversion */	
	coo2csc(csc_r,csc_c,coo_rows,coo_cols,nnz,row,isOneBased);
	/*Allocate Memory for Matix Stuct */
	cscMat Matrix; 
	Matrix.csc_r = csc_r;
	Matrix.csc_c = csc_c;
	Matrix.nnz = nnz;
	Matrix.rows = row;
	/* Freeing space */
	/*
	   free(csc_r);
	   free(csc_c);
	   */
	free(coo_rows);
	free(coo_cols);
	/* Return Function */
	return Matrix;
}

void
checkArgsNum(int argc){
	switch(argc){
		case 1:
			fprintf(stderr, "Line:%d No argument was given\n",__LINE__);
			exit(EXIT_FAILURE);
		case 2:
			fprintf(stderr, "Line:%d One argument is missing\n",__LINE__);
			exit(EXIT_FAILURE);
		case 3:
			break;
		default:
			fprintf(stderr, "Line:%d More args given\n",__LINE__);
			exit(EXIT_FAILURE);
	}
}
void coo2csc(
		int       * const row,       /*!< CSC row start indices */
		int       * const col,       /*!< CSC column indices */
		int const * const row_coo,   /*!< COO row indices */
		int const * const col_coo,   /*!< COO column indices */
		int const nnz,       /*!< Number of nonzero elements */
		int const n,         /*!< Number of rows/columns */
		int const         isOneBased /*!< Whether COO is 0- or 1-based */
	    ) {

	// ----- cannot assume that input is already 0!
	for (int l = 0; l < n+1; l++) col[l] = 0;


	// ----- find the correct column sizes
	for (int l = 0; l < nnz; l++)
		col[col_coo[l] - isOneBased]++;

	// ----- cumulative sum
	for (int i = 0, cumsum = 0; i < n; i++) {
		int temp = col[i];
		col[i] = cumsum;
		cumsum += temp;
	}
	col[n] = nnz;
	// ----- copy the row indices to the correct place
	for (int l = 0; l < nnz; l++) {
		int col_l;
		col_l = col_coo[l] - isOneBased;

		int dst = col[col_l];
		row[dst] = row_coo[l] - isOneBased;

		col[col_l]++;
	}
	// ----- revert the column pointers
	for (int i = 0, last = 0; i < n; i++) {
		int temp = col[i];
		col[i] = last;
		last = temp;
	}

}
