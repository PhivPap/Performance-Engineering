#ifndef MATRIX_PARSER
#define MATRIX_PARSER
#include <stdio.h>

void generate_mat(int m, int n, int p, float *A, float *B);
void read_sparse(FILE *f, int m, int n, int nz, float *A);
void write_sparse(FILE *f, int m, int p, const float *C);
void read_dense(FILE *f, int m, int n, float *A);
int read_mat(int *m, int *n, int *p, int *nzA, int *nzB, FILE *fa, FILE *fb);

#endif