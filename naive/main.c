#include<stdio.h>
#include<stdlib.h>
#include<errno.h>
#include<string.h>
#include"inc/mmio.h"
extern int errno ;
/* Function Predeclarations */
void
checkArgsNum(int argc);
/* Temporary placed function */
/* Convert COO to CSC  */
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
 
int
main(int argc, char *argv[]){
  checkArgsNum(argc);
  FILE *f;
  int M,N,nz,ret;
  MM_typecode matcode;
  /* Open File and Handle errors */
  f= fopen (argv[1],"r");
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
  /* Conversion */	
  coo2csc(csc_r,csc_c,coo_rows,coo_cols,nnz,row,isOneBased);
  /*for (int i =0 ;i<nnz;i++) printf("%d ",csc_r[i]);
  puts("\n");
  for (int i =0 ;i<row+1;i++) printf("%d ",csc_c[i]);*/


}



void
checkArgsNum(int argc){
  switch(argc){
    case 1:
      fprintf(stderr, "Line:%d No argument was given\n",__LINE__);
      exit(EXIT_FAILURE);
    case 2:
      break;
    default:
      fprintf(stderr, "Line:%d More args given\n",__LINE__);
      exit(EXIT_FAILURE);
  }
}
