/* This method uses 2 CSC matrices and uses a better method by iterating the CSC matrices */
#define _POSIX_C_SOURCE 199309L
#include<stdio.h>
#include<stdlib.h>
#include<errno.h>
#include<string.h>
#include<sys/time.h>
#include<time.h>
#include"inc/mmio.h"
#include<omp.h>
#include<math.h>
extern int errno ;
/* Structs Declaration */
typedef struct cscMatrix{
  int *csc_r;
  int *csc_c;
  int nnz;
  int rows;
}cscMat;
typedef struct cooMatrix{
  int *coo_r;
  int *coo_c;
  int nnz;
  int rows;
}cooMat;
/* Function Predeclarations */
void checkArgsNum(int argc);
void coo2csc(int * const row,int * const col,int const * const row_coo,
             int const * const col_coo,int const nnz,int const n,int const isOneBased);
cscMat readFile(char *filename);
cooMat sortMat(cooMat Matrix);

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
  /* This currently works for matrices RowsxRows */
  int block = 30000;
  for(int j=0;j<MatC.rows;j+=block){
    //if (j%1000==0) printf("%d\n",j);
    //struct timespec start;
    //struct timespec end;
    //clock_gettime(CLOCK_MONOTONIC,&start);
    int *hits = (int *)calloc(MatC.rows*block,sizeof(int));
    if(hits==NULL){
      fprintf(stderr,"Line %d: Error Alocating Matrix hits",__LINE__);
      exit(EXIT_FAILURE);
    }
	printf("%d\n",j);
    /* Iterate B column */
    int currj=0;
    #pragma omp parallel for schedule(dynamic) private(currj) reduction(+:hits[:MatC.rows])
    for (int c=MatB.csc_c[j];c<MatB.csc_c[j+block];c++){
      int MatArow=MatB.csc_r[c];
      for (int r=MatA.csc_c[MatArow];r<MatA.csc_c[MatArow+1];r++){
	int i = MatA.csc_r[r];
	hits[i+MatC.rows*currj]+=1;
      }
      currj++;
    }
    // TODO: Fix the exit memory allocation
    for(int i=0;i<MatC.rows*block;i++){
      if(hits[i]>0){
      MatC.nnz ++;
      MatC.coo_r = (int *)realloc(MatC.coo_r,MatC.nnz*sizeof(int));
      MatC.coo_c = (int *)realloc(MatC.coo_c,MatC.nnz*sizeof(int));
      if (MatC.coo_c == NULL || MatC.coo_r==NULL)
        {fprintf(stderr,"Line %d: Error Alocating Matrix C",__LINE__);
          exit(EXIT_FAILURE);}
      MatC.coo_c[MatC.nnz-1] = j+j%block;
      MatC.coo_r[MatC.nnz-1] = i-(j%block)*MatC.rows;
      }
    }
    free(hits);
  //clock_gettime(CLOCK_MONOTONIC,&end);
  //double sec,nsec;
  //sec = end.tv_sec - start.tv_sec;
  //nsec = end.tv_nsec - start.tv_nsec;
  //printf("Execution Time %f ms\n",sec*1000+nsec/1000000);
  }
  clock_gettime(CLOCK_MONOTONIC,&ts_end);
  double ts_sec,ts_nsec;
  ts_sec = ts_end.tv_sec - ts_start.tv_sec;
  ts_nsec = ts_end.tv_nsec - ts_start.tv_nsec;
  printf("Execution Time %f ms\n",ts_sec*1000+ts_nsec/1000000);
  /* Sort MatC by column for tests  */
  MatC = sortMat(MatC);

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

//  puts("\nMATRIX C\n");
//  for (int i =0 ;i<MatC.nnz;i++) printf("%d %d\n",MatC.coo_r[i],MatC.coo_c[i]);



  /* Free csc arrays allocated from readFile(...)  */
  free(MatA.csc_r);
  free(MatA.csc_c);
  free(MatB.csc_r);
  free(MatB.csc_c);
  free(MatC.coo_r);
  free(MatC.coo_c);
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
  int const         nnz,       /*!< Number of nonzero elements */
  int const         n,         /*!< Number of rows/columns */
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
