/* SpMV: this file contains the I/O functions used to read and write matrices
 * in the Market Matrix format(see
 * https://math.nist.gov/MatrixMarket/formats.html#MMformat), using the
 * functions in mmio.c in its turn.
 * A vector of the appropriate size is generated and multiplied with the matrix.
 * There are also functions to generate your own matrices.
 * The reading also supports the "is_pattern" flag from the Matrix Market format. 
 */
#include <iostream>
#include <chrono>
#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include "mmio.h"

#define check(error) checkCudaCall(error, __LINE__)

//#define VERBOSE

#define N  512
#define M  512

#define REP 100

void checkCudaCall(cudaError_t result, uint32_t line)
{
    if (result != cudaSuccess)
    {
        printf("cuda error \n");
        printf("Line %u: %s\n", line, cudaGetErrorString(result));
        fflush(stdout);
        exit(1);
    }
}

__global__ void csc(const int *A_cols, const int *A_rows_idx, const float *A_values, const float *B, float *C) {
    int i,j;
    int col_start, col_end; 

    i = blockIdx.x;  

        col_start = A_cols[i];
        col_end = A_cols[i + 1];
        
        for (j = col_start; j < col_end; j++) {
            C[A_rows_idx[j]] += A_values[j] * B[i];
        }
}

void gpu_memory_init(int m, int n, int nzA, int *A_cols, int *A_rows_idx, float *A_values, float *B, float *C,
    int **d_A_cols, int **d_A_rows_idx, float **d_A_values, float **d_B, float **d_C)
{
    const uint32_t A_cols_size = sizeof(int) * (n + 1);
    const uint32_t A_rows_idx_size = sizeof(int) * nzA;
    const uint32_t A_values_size = sizeof(float) * nzA;
    const uint32_t B_size = sizeof(float) * n;
    const uint32_t C_size = sizeof(float) * m; 

    check(cudaMalloc(d_A_cols, A_cols_size));
    check(cudaMalloc(d_A_rows_idx, A_rows_idx_size));
    check(cudaMalloc(d_A_values, A_values_size));
    check(cudaMalloc(d_B, B_size));
    check(cudaMalloc(d_C, C_size));

    check(cudaMemcpy(*d_A_cols, A_cols, A_cols_size, cudaMemcpyHostToDevice));
    check(cudaMemcpy(*d_A_rows_idx, A_rows_idx, A_rows_idx_size, cudaMemcpyHostToDevice));
    check(cudaMemcpy(*d_A_values, A_values, A_values_size, cudaMemcpyHostToDevice));
    check(cudaMemcpy(*d_B, B, B_size, cudaMemcpyHostToDevice));
    check(cudaMemcpy(*d_C, C, C_size, cudaMemcpyHostToDevice));
}

void gpu_memory_free(int *d_A_cols, int *d_A_rows_idx, float *d_A_values, float *d_B, float *d_C){
    check(cudaFree(d_A_cols));
    check(cudaFree(d_A_rows_idx));
    check(cudaFree(d_A_values));
    check(cudaFree(d_B));
    check(cudaFree(d_C));
}

void csc_spmv(int m, int n, int nzA, int *A_cols, int *A_rows_idx, float *A_values, float *B, float *C) {
    int r;

    int *d_A_cols, *d_A_rows_idx;
    float *d_A_values;
    float *d_B;
    float *d_C;

    gpu_memory_init(m, n, nzA, A_cols, A_rows_idx, A_values, B, C, &d_A_cols, &d_A_rows_idx, &d_A_values, &d_B, &d_C);

    const auto start = std::chrono::system_clock::now();
    for (r=0; r<REP; r++) 
    {
        csc<<<n, 1>>>(d_A_cols, d_A_rows_idx, d_A_values, d_B, d_C);
        cudaDeviceSynchronize();
    }
    const auto end = std::chrono::system_clock::now();
    cudaMemcpy(C, d_C, sizeof(float) * m, cudaMemcpyDeviceToHost);

    using std::chrono::duration_cast;
    using std::chrono::milliseconds;
    const double avg_execution_time = duration_cast<milliseconds>(end - start).count() / 1000.0 / REP;
    printf("Seconds: %f\n", avg_execution_time);
    printf("REP: %i\n", REP);
    
    gpu_memory_free(d_A_cols, d_A_rows_idx, d_A_values, d_B, d_C);
}

/* 
 * for nz=1 - dense matrix; for nz>1 - sparse matrix with every nz element non-0; for nz=0: error
 */
void generate_mat(int m, int n, float *A, int nz) {
  int i;

  for (i=0; i<(m*n); i++) 
	A[i] = i/nz; //i/10; 
      	
}

void read_sparse(FILE *f, int m, int n, int nz, float *A, int is_pattern) {
  int i, row, col;
  float val;  
 
    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    if (is_pattern) {
        for (i=0; i<nz; i++)
        {
            fscanf(f, "%d %d\n", &row, &col);
            A[(row-1)*n+col-1] = 1.0f;   /* adjust from 1-based to 0-based */
        }
    }
    else {
      for (i=0; i<nz; i++)
      {
          fscanf(f, "%d %d %f\n", &row, &col, &val);
          A[(row-1)*n+col-1] = val;   /* adjust from 1-based to 0-based */
      }
    }
}

void write_sparse(FILE *f, int m, int p, const float *C) {
   int i, nz=0; 
   MM_typecode matcode;

   for (i=0; i<m*p; i++) if (C[i] != 0.0) nz++; 

    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    mm_set_coordinate(&matcode);
    mm_set_real(&matcode);

    mm_write_banner(f, matcode); 
    mm_write_mtx_crd_size(f, m, p, nz);

    for (i=0; i<m*p; i++) {
	if (C[i] != 0.0) 
          fprintf(f, "%d %d %f\n", i/p+1, i%p+1, C[i]);
    }

}

void write_vector(FILE* f, int m, const float* C) {
  int i;

  for (i=0; i<m; i++) {
      fprintf(f, "%f\n", C[i]);
  }
}

void read_dense(FILE *f, int m, int n, float *A) {
  int row, col;

  for(row=0; row<m; row++) { 
     for (col=0; col<n; col++) {
        fscanf(f, "%f ", &A[row*n+col]); 
     }
  } 
}

void print_mat(int m, int n, float *A) {
  int row, col; 

  for(row=0; row<m; row++) {
     for (col=0; col<n; col++) {
     	printf("%10.5f", A[row*n+col]);
     }
     printf("\n"); 
  }
}


int read_mat(int *m, int *n, int *nzA, FILE* fa, int *is_pattern) {
  MM_typecode ta;
  int ret_code = 0; 

  if (mm_read_banner(fa, &ta) != 0)
    {
        printf("Could not process Matrix Market banner for A.\n");
        return -3;
    }

  if (mm_is_pattern(ta)) 
	*is_pattern = 1;
  else
	*is_pattern = 0;


  if (mm_is_complex(ta)) return -6;

  if (mm_is_matrix(ta) && mm_is_sparse(ta))
    {
        if ((ret_code = mm_read_mtx_crd_size(fa, m, n, nzA)) !=0)
           return -10;
    }
  else if (mm_is_matrix(ta) && mm_is_array(ta)) {
	*nzA = 0;
        if ((ret_code = mm_read_mtx_array_size(fa, m, n)) !=0)
           return -11;

    }
  else return -8; 

  return ret_code;
}

/*
 *  * Converts matrix to CSC format. 
 *   * returns the number of nonZero's found during conversion
 *    */
int convert_to_csc(int m, int n, float *A, int *sA_cols, int *sA_row_idx, float *sA_vals) {
  int i,j;
  int checkNZ=0;
  float tmp;

  sA_cols[0]=0;
  checkNZ=0;
  for (i=0; i<n; i++) {
    for (j=0; j<m; j++) {
        tmp = A[j*n+i];
        if (tmp != 0) {
           sA_row_idx[checkNZ]=j;
           sA_vals[checkNZ]=tmp;
           checkNZ++;
        }
    }
    sA_cols[i+1]=checkNZ;
  }
  return checkNZ;
}

void print_mat_csc(int n, int nzA, int *sA_cols, int *sA_row_idx, float *sA_vals)
{
  int i;
  for (i=0; i<n+1; i++)
    printf("%d ", sA_cols[i]);
  printf("\n");

  for (i=0; i<nzA; i++)
    printf("%d ", sA_row_idx[i]);

  printf("\n");

  for (i=0; i<nzA; i++)
    printf("%10.5f ", sA_vals[i]);
  printf("\n");

}


int main (int argc, char** argv) {
 float *A, *B, *C;
 float *sA_vals;
 int *sA_cols, *sA_rows_idx;

 int m, n, err;
 int nzA=0, is_pattern = 1;
 FILE *fa, *fc;
  
#ifdef GENERATE 
 m=M; n=N; nzA=M*N/10; 
#else 
 if (argc < 3) {
    fprintf(stderr, "Usage: %s [matrix-market-filename] [result-vector-filename]\n", argv[0]);
    exit(1);
 }
 else {
    if ((fa = fopen(argv[1], "rt")) == NULL) exit(1);
    err = read_mat(&m, &n, &nzA, fa, &is_pattern);    
    if (err == -15) {
	printf("Matrices are incompatible! \n");
	fclose(fa); 
	exit(1);
    }
 }
#endif

 A = (float *)calloc(m*n,sizeof(float));
 if (A==NULL) {printf("Out of memory A! \n"); exit(1);}

#ifdef GENERATE
   generate_mat(m,n,A,nzA);
#else 
   if (nzA>0) {
	read_sparse(fa, m,n,nzA, A, is_pattern);
    }	
   else 
	read_dense(fa, m,n, A);
   fclose(fa);

        sA_cols = (int *)calloc(n+1,sizeof(int));
        sA_rows_idx = (int *)calloc(nzA,sizeof(int));
        sA_vals = (float *)calloc(nzA,sizeof(float));

        convert_to_csc(m,n, A, sA_cols, sA_rows_idx, sA_vals);
#ifdef VERBOSE
        print_mat(m,n,A);
        print_mat_csc(m, nzA, sA_cols, sA_rows_idx, sA_vals);
#endif

#endif

 B = (float *)calloc(n,sizeof(float));
 if (B==NULL) {printf("Out of memory B! \n"); exit(1);}

   generate_mat(n,1,B,1);
#ifdef VERBOSE
   print_mat(n,1,B);
#endif 

 C = (float *)calloc(m,sizeof(float));
 if (C==NULL) {printf("Out of memory C! \n"); exit(1);}

 /* Call the SpMV kernel. */
  csc_spmv(m, n, nzA, sA_cols, sA_rows_idx, sA_vals, B, C);

 if ((fc = fopen(argv[2], "wt")) == NULL) exit(3);    
// write_sparse(fc,n,m,C);
 write_vector(fc,m,C);
 fclose(fc);  

 free(A);
 free(B); 
 free(C);

}

