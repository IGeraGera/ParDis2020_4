/* This is implementaion uses MPI without dissecting the matrix with F mask matrix */
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
	int size;
}cooMat;
/* Function Predeclarations */
void checkArgsNum(int argc);
void coo2csc(int * const row,int * const col,int const * const row_coo,
		int const * const col_coo,int const nnz,int const n,int const isOneBased);
cscMat readFile(char *filename);
cooMat sortMat(cooMat Matrix);
void matrixMult(cscMat MatA, cscMat MatB, cscMat MatF, cooMat *MatC, int **totalHits, int *totalHitsSize, int firstFlag);
void checkCOOMat(cooMat *Mat,int idx);
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
	/* Get MatF */
	cscMat MatF;
	MatF = readFile(argv[3]);
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
	MPI_Init(NULL,NULL);
	MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	/* Block size for each axis */
	int block =(int)ceil(((float)MatC.rows)/((float)numtasks)) ;
	/* Assign the workload to every thread */
	/* Every process gets a block of columns to calculate from the matrix C */
	int workStart,workEnd;
	workStart = block*rank;
	workEnd = workStart+block;
	if(workEnd>MatC.rows) workEnd=MatC.rows;
	int totalWorkload=workEnd-workStart;
	printf("rank %d start %d end %d\n",rank,workStart,workEnd);

	struct cooMatrix MatrixCOOArrC;
	MatrixCOOArrC.nnz=0;
	MatrixCOOArrC.coo_c=NULL;
	MatrixCOOArrC.coo_r=NULL;
	/* Start Calculating the result */
	/* Assign columns of C to processes */
	for(int j=workStart;j<workEnd;j++){
		int * finalHits = NULL;
		int nnz=0;
		for (int c=MatB.csc_c[j];c<MatB.csc_c[j+1];c++){
			int MatArow=MatB.csc_r[c];
			/* Init F column pointer */
			int fpointer = MatF.csc_c[j];
			int fpointerend = MatF.csc_c[j+1];
			for (int r=MatA.csc_c[MatArow];r<MatA.csc_c[MatArow+1];r++){
				/* Get row of MatA */
				int i = MatA.csc_r[r];
				/* Check iff the pointer reached the end  */
				if (fpointer==fpointerend) break;	
				/* if the row is lesser than the row the pointer points continue */
				else if (i<MatF.csc_r[fpointer]) continue;
				/* if the row is greater than the row the pointer points increment the pointer */
				while(i>MatF.csc_r[fpointer] && fpointer<fpointerend) fpointer++;
				/* Check again if the pointer passed the end  */
				if (fpointer==fpointerend) break;	
				/* if there is no hit and the row is in F the add it  */
				if(MatF.csc_r[fpointer]==i){
					for(int hit=0;hit<nnz;hit++){
						if(i==finalHits[hit]) break;
					}
					nnz++;
					finalHits = (int *)realloc(finalHits,nnz*sizeof(int));
					finalHits[nnz-1] = i;
					MatrixCOOArrC.nnz++;
					MatrixCOOArrC.coo_c=(int *)realloc(MatrixCOOArrC.coo_c,MatrixCOOArrC.nnz*sizeof(int));
					MatrixCOOArrC.coo_r=(int *)realloc(MatrixCOOArrC.coo_r,MatrixCOOArrC.nnz*sizeof(int));
					if (MatrixCOOArrC.coo_c == NULL || MatrixCOOArrC.coo_r==NULL)
					{fprintf(stderr,"Line %d: Error Alocating MatCOOArrC",__LINE__);
						exit(EXIT_FAILURE);}
					MatrixCOOArrC.coo_c[MatrixCOOArrC.nnz-1]=j;
					MatrixCOOArrC.coo_r[MatrixCOOArrC.nnz-1]=i;

				}
			}
		}
		free(finalHits);
	}

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
	if(rank==0){
		MatC = sortMat(MatC);
		for (int i=0 ; i<MatC.nnz;i++) printf("Num %d i %d j %d \n",i+1,MatC.coo_r[i]+1,MatC.coo_c[i]+1);
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
	free(MatrixCOOArrC.coo_c);
	free(MatrixCOOArrC.coo_r);
	/* Free csc arrays allocated from readFile(...)  */
	free(MatA.csc_r);
	free(MatA.csc_c);
	free(MatB.csc_r);
	free(MatB.csc_c);
	free(MatF.csc_r);
	free(MatF.csc_c);
	free(MatC.coo_r);
	free(MatC.coo_c);
}
/* Function that returns the min array */
int
imin(int a,int b){
	if(a>=b) {return b;}
	else {return a;}
}
void
checkCOOMat(cooMat *Mat,int idx){
	/* Check If the Size is larger that the index and allocate more space*/
	if (Mat->size==(idx+1) ){
		Mat->size*=2;
		Mat->coo_r = (int *) realloc(Mat->coo_r,Mat->size*sizeof(int));
		Mat->coo_c = (int *) realloc(Mat->coo_c,Mat->size*sizeof(int));
		if (Mat->coo_c == NULL || Mat->coo_r==NULL)
			{fprintf(stderr,"Line %d: Error reallocating coo matrix",__LINE__);
				exit(EXIT_FAILURE);}

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
			fprintf(stderr, "Line:%d Two arguments are missing\n",__LINE__);
			exit(EXIT_FAILURE);
		case 3:
			fprintf(stderr, "Line:%d One argument is missing\n",__LINE__);
			exit(EXIT_FAILURE);
		case 4:
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
