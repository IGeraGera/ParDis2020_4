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
void matrixMult(cscMat MatA, cscMat MatB, cooMat *MatC);

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
	double ts_sec,ts_nsec;
	clock_gettime(CLOCK_MONOTONIC,&ts_start);


	/* MPI */
	int numtasks,rank;
	struct cooMatrix *MatrixCOOArrA, *MatrixCOOArrB, MatrixCOOArrC;
	struct cscMatrix *MatrixArrA, *MatrixArrB;

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
	printf("rank %d start %d end %d\n",rank,workStart,workEnd);

	clock_gettime(CLOCK_MONOTONIC,&ts_end);
	ts_sec = ts_end.tv_sec - ts_start.tv_sec;
	ts_nsec = ts_end.tv_nsec - ts_start.tv_nsec;
	printf("rank %d Time till workload %f ms\n",rank,ts_sec*1000+ts_nsec/1000000);

	/* Allocate an arrays that will contain the  blocks COO and CSC  */
	MatrixCOOArrA = malloc(blockalloc*blockalloc*sizeof(struct cooMatrix));
	MatrixCOOArrB = malloc(blockalloc*blockalloc*sizeof(struct cooMatrix));
	if (MatrixCOOArrA==NULL || MatrixCOOArrB==NULL) {
		printf("Memory allocation at line %d Failed\n",__LINE__);
		exit(EXIT_FAILURE);
	}
	for (int i=0;i<blockalloc*blockalloc;i++){
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
	for(int j=0;j<MatA.rows;j++){
		for(int r=MatA.csc_c[j];r<MatA.csc_c[j+1];r++){
			int i = MatA.csc_r[r];
			/* Find Arri and Arrj  */
			int Arri = (int) floor(i/MatrixCOOArrA[0].rows);
			int Arrj = (int) floor(j/MatrixCOOArrA[0].rows);
			//printf("i %d j %d Arri %d Arrj %d \n",i,j,Arri,Arrj);

			/* Check if structure exist in MatrixCOOArrA and allocate it */	  
			int ptr = Arri +blockalloc*Arrj;
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
	for(int j=0;j<MatB.rows;j++){
		for(int r=MatB.csc_c[j];r<MatB.csc_c[j+1];r++){
			int i = MatB.csc_r[r];
			/* Find Arri and Arrj  */
			int Arri = (int) floor(i/MatrixCOOArrB[0].rows);
			int Arrj = (int) floor(j/MatrixCOOArrB[0].rows);
			//printf("Arri %d Arrj %d for i %d j %d \n",Arri,Arrj,i,j);

			/* Check if structure exist in MatrixCOOArrB and allocate it */	  
			int ptr = Arri +blockalloc*Arrj;
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
	clock_gettime(CLOCK_MONOTONIC,&ts_end);
	ts_sec = ts_end.tv_sec - ts_start.tv_sec;
	ts_nsec = ts_end.tv_nsec - ts_start.tv_nsec;
	printf("Rank %d Time to disect %f ms\n",rank,ts_sec*1000+ts_nsec/1000000);
	/* Covert All COO to CSC */
	//TODO: Remember to free
	MatrixArrA =  malloc(blockalloc*blockalloc*sizeof(struct cscMatrix));
	MatrixArrB = malloc(blockalloc*blockalloc*sizeof(struct cscMatrix));
	if (MatrixArrA==NULL || MatrixArrB==NULL) {
		printf("Memory allocation at line %d Failed\n",__LINE__);
		exit(EXIT_FAILURE);
	}
	for(int i=0;i<blockalloc;i++){
		for(int j=0;j<blockalloc;j++){
			int ptr = i + blockalloc*j;
			/* Declaration of CSC format arrays  */ 
			int row = MatrixCOOArrA[ptr].rows;
			int nnz = MatrixCOOArrA[ptr].nnz;
			int isOneBased =0; 
			/*Allocate memory for csc  matrices */
			MatrixArrA[ptr].nnz = nnz;
			MatrixArrA[ptr].rows =row;
			MatrixArrA[ptr].i=i*row;
			MatrixArrA[ptr].j=j*row;
			MatrixArrA[ptr].csc_r = (int *)malloc(nnz     * sizeof(int));
			MatrixArrA[ptr].csc_c = (int *)malloc((row+1) * sizeof(int));
			if (MatrixArrA[ptr].csc_r==NULL || MatrixArrA[ptr].csc_c==NULL) {
				printf("Memory allocation at line %d Failed\n",__LINE__);
				exit(EXIT_FAILURE);
			}
			/* Conversion */	
			if(MatrixCOOArrA[ptr].nnz==0){	// If the array is empty nnz=0 then the arrays are NULL
				MatrixArrA[ptr].csc_r = NULL;
				MatrixArrA[ptr].csc_c = NULL;
			}
			else{
				coo2csc(MatrixArrA[ptr].csc_r,MatrixArrA[ptr].csc_c,
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
			MatrixArrB[ptr].j=j*row;
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
	clock_gettime(CLOCK_MONOTONIC,&ts_end);
	ts_sec = ts_end.tv_sec - ts_start.tv_sec;
	ts_nsec = ts_end.tv_nsec - ts_start.tv_nsec;
	printf("Rank %d Time till COO to CSC %f ms\n",rank,ts_sec*1000+ts_nsec/1000000);
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
			/* if (rank==0) printf("Arri %d Arrj %d rows %d\n",Crow,Ccol,MatrixArrA[Cptr].rows); */
			/* Iterate the blocks in respect to row and column */
			for(int Ck=0;Ck<blockalloc;Ck++){
				int Aptr = Crow+blockalloc*Ck;
				int Bptr = Ck+blockalloc*Ccol;
				/* Temporary Matrix with result */
				cooMat tempC;
				tempC.coo_r = NULL;
				tempC.coo_c = NULL;
				tempC.nnz = 0;
				// If block is empty skip
				if (MatrixArrA[Aptr].nnz>0 && MatrixArrB[Bptr].nnz>0){
					/* Multiply the arrays A x B = Ctemp */
					matrixMult(MatrixArrA[Aptr],MatrixArrB[Bptr],&tempC);
					// DEBUG
					/* if (rank ==0){ */
					/* for (int i=0;i<tempC.nnz;i++){ */
					/* 	if(rank==0) printf("Arri %d Arrj %d i %d j %d\n",Crow,Ccol,tempC.coo_r[i]+MatrixArrA[Cptr].i+1,tempC.coo_c[i]+MatrixArrA[Cptr].j+1); */
					/* } */
					/* } */
					/* Merge the results of temporary matrix to MatrixCOOArrC*/
					/* Check if tempSubC is empty and the tempC loaded and pass the results to it */
					if (tempSubC.nnz == 0 && tempC.nnz!=0){
						tempSubC.coo_r = tempC.coo_r;
						tempSubC.coo_c = tempC.coo_c;
						tempSubC.nnz   = tempC.nnz;
						/* Sort the Matrix */
						sortMat(tempSubC);

					} 
					/* Check if both MatrixCOOArrC and tempC are empty and continue */
					else if ( tempC.nnz == 0 ) continue;
					/* If tempC and tempSubC are loaded merge the 2 of them */
					else {
						/* Sort the tempC matrix by columns */
						sortMat(tempC);
					        int *rows =NULL;
						int *cols = NULL;
						int nnz = 0;	
						int tempCin=0;
						int tempSubCin=0;

						/* The matrix is sorted by columns and rows */
						while(1){
							/* Check if any mat finished and copy the remaining elements */
							if (tempCin==tempC.nnz){
								int start = nnz;
								nnz+=tempSubC.nnz-tempSubCin;
								rows = realloc(rows,nnz*sizeof(int));
								cols = realloc(cols,nnz*sizeof(int));
								for (int i =0;i<tempSubC.nnz-tempSubCin;i++){
									rows[start+i] = tempSubC.coo_r[tempSubCin+i];
									cols[start+i] = tempSubC.coo_c[tempSubCin+i];
								}
								break;
							}
							if (tempSubCin==tempSubC.nnz){
								int start = nnz;
								nnz+=tempC.nnz-tempCin;
								rows = realloc(rows,nnz*sizeof(int));
								cols = realloc(cols,nnz*sizeof(int));
								for (int i =0;i<tempC.nnz-tempCin;i++){
									rows[start+i] = tempC.coo_r[tempCin+i];
									cols[start+i] = tempC.coo_c[tempCin+i];
								}
								break;
							}
							if (tempC.coo_c[tempCin]<tempSubC.coo_c[tempSubCin]){
								nnz++;
								rows = realloc(rows,nnz*sizeof(int));
								cols = realloc(cols,nnz*sizeof(int));
								rows[nnz-1] = tempC.coo_r[tempCin];
								cols[nnz-1] = tempC.coo_c[tempCin];
								tempCin++;
							}
							else if (tempC.coo_c[tempCin]>tempSubC.coo_c[tempSubCin]){
								nnz++;
								rows = realloc(rows,nnz*sizeof(int));
								cols = realloc(cols,nnz*sizeof(int));
								rows[nnz-1] = tempSubC.coo_r[tempSubCin];
								cols[nnz-1] = tempSubC.coo_c[tempSubCin];
								tempSubCin++;
							}
							else if (tempC.coo_c[tempCin]==tempSubC.coo_c[tempSubCin]){
								/* Check the row  */
								if (tempC.coo_r[tempCin]<tempSubC.coo_r[tempSubCin]){
									nnz++;
									rows = realloc(rows,nnz*sizeof(int));
									cols = realloc(cols,nnz*sizeof(int));
									rows[nnz-1] = tempC.coo_r[tempCin];
									cols[nnz-1] = tempC.coo_c[tempCin];
									tempCin++;
								}
								else if (tempC.coo_r[tempCin]==tempSubC.coo_r[tempSubCin]){
									nnz++;
									rows = realloc(rows,nnz*sizeof(int));
									cols = realloc(cols,nnz*sizeof(int));
									rows[nnz-1] = tempSubC.coo_r[tempSubCin];
									cols[nnz-1] = tempSubC.coo_c[tempSubCin];
									tempSubCin++;
									tempCin++;}
								else{
									nnz++;
									rows = realloc(rows,nnz*sizeof(int));
									cols = realloc(cols,nnz*sizeof(int));
									rows[nnz-1] = tempSubC.coo_r[tempSubCin];
									cols[nnz-1] = tempSubC.coo_c[tempSubCin];
									tempSubCin++;
								}
							}
						} 
						/* Free tempC coo matrices */
						free(tempC.coo_r);
						free(tempC.coo_c);
						free(tempSubC.coo_r);
						free(tempSubC.coo_c);
						tempSubC.coo_r =rows;
						tempSubC.coo_c =cols;
						tempSubC.nnz=nnz;
					}
				} 
			}
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
	clock_gettime(CLOCK_MONOTONIC,&ts_end);
	ts_sec = ts_end.tv_sec - ts_start.tv_sec;
	ts_nsec = ts_end.tv_nsec - ts_start.tv_nsec;
	printf("Rank %d Time till calc %f ms\n",rank,ts_sec*1000+ts_nsec/1000000);
	//sortMat(MatrixCOOArrC);
	clock_gettime(CLOCK_MONOTONIC,&ts_end);
	ts_sec = ts_end.tv_sec - ts_start.tv_sec;
	ts_nsec = ts_end.tv_nsec - ts_start.tv_nsec;
	printf("Rank %d Time till COO sort %f ms\n",rank,ts_sec*1000+ts_nsec/1000000);
	/* if(rank==2){ */
	/* 	puts("\nMATRIX C\n"); */
	/* 	for (int i =0 ;i<MatrixCOOArrC.nnz;i++) printf("%d %d\n",MatrixCOOArrC.coo_r[i]+1,MatrixCOOArrC.coo_c[i]+1); */
	/* } */ 
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
	ts_sec = ts_end.tv_sec - ts_start.tv_sec;
	ts_nsec = ts_end.tv_nsec - ts_start.tv_nsec;
	printf("Rank %d Execution Time %f ms\n",rank,ts_sec*1000+ts_nsec/1000000);
	/* Sort MatC by column for tests  */
	/* if(rank==0){ */
	/* 	MatC = sortMat(MatC); */
	/* 	for (int i=0 ; i<MatC.nnz;i++) printf("Num %d i %d j %d \n",i+1,MatC.coo_r[i]+1,MatC.coo_c[i]+1); */
	/* } */

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

/* This routine Multiplies 2 csc matrixes */
void 
matrixMult(cscMat MatA, cscMat MatB, cooMat *MatC){
	/* This currently works for matrices RowsxRows */
	for(int j=0;j<MatA.rows;j++){
		//if (j%1000==0) printf("%d\n",j);
		int *hits = (int *)calloc(MatA.rows,sizeof(int));
		if(hits==NULL){
			fprintf(stderr,"Line %d: Error Alocating Matrix hits",__LINE__);
			exit(EXIT_FAILURE);
		}
		/* Iterate B column */
		for (int c=MatB.csc_c[j];c<MatB.csc_c[j+1];c++){
			int MatArow=MatB.csc_r[c];
			for (int r=MatA.csc_c[MatArow];r<MatA.csc_c[MatArow+1];r++){
				/* check if exist */
				int i = MatA.csc_r[r];
				if (hits[i]==0){
					hits[i]=1;
					MatC->nnz ++;
					MatC->coo_r = (int *)realloc(MatC->coo_r,MatC->nnz*sizeof(int));
					MatC->coo_c = (int *)realloc(MatC->coo_c,MatC->nnz*sizeof(int));
					if (MatC->coo_c == NULL || MatC->coo_r==NULL)
					{fprintf(stderr,"Line %d: Error Alocating Matrix C",__LINE__);
						exit(EXIT_FAILURE);}
					MatC->coo_r[MatC->nnz-1] = i;
					MatC->coo_c[MatC->nnz-1] = j;

				}

			}
		}
		free(hits);
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
