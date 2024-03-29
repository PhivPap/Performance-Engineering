#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "mmio.h"
#include "MatrixParser/matrixparser.h"

#define N 512
#define M 512
#define P 512

#define REP 5


// m = rows of C and rows of A
// p = colums of C and colums of B
// n = colums of A and rows of B
// MATRIX = [ROWS][COLUMS]
// A = [m][n]
// B = [n][p]
// C = [m][p]
void matrix_mult(int m, int n, int p, float *A, float *B, float *C)
{
    int i, j, k;

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < p; j++)
        {
            C[i * p + j] = 0;
            for (k = 0; k < n; k++)
            {
                C[i * p + j] += A[i * n + k] * B[k * p + j];
            }
        }
    }
}


double get_operation_count(int m, int n, int p)
{
    return (double)m * n * p * 2; // exclude overhead 9 + (double)m * p * 3 + m;
}

int main(int argc, char **argv)
{
    float *A, *B, *C;
#ifdef TIMING
    struct timeval before, after;
#endif
    int m, n, p, r, err;
    int nzA = 0, nzB = 0;
    FILE *fa, *fb, *fc;

#ifdef GENERATE
    m = M;
    n = N;
    p = P;
#else
    if (argc < 3)
    {
        fprintf(stderr, "Usage: %s [martix1] [matrix2] [resultmatrix] \n", argv[0]);
        exit(1);
    }
    else
    {
        if ((fa = fopen(argv[1], "rt")) == NULL)
            exit(1);
        if ((fb = fopen(argv[2], "rt")) == NULL)
            exit(2);
        err = read_mat(&m, &n, &p, &nzA, &nzB, fa, fb);
        if (err == -15)
        {
            printf("Matrices are incompatible! \n");
            fclose(fa);
            fclose(fb);
            exit(1);
        }
    }
#endif

    A = (float *)calloc(m * n, sizeof(float));
    if (A == NULL)
    {
        printf("Out of memory A! \n");
        exit(1);
    }
    B = (float *)calloc(n * p, sizeof(float));
    if (B == NULL)
    {
        printf("Out of memory B! \n");
        exit(1);
    }

#ifdef GENERATE
    generate_mat(m, n, p, A, B);
#else
    if (nzA > 0)
        read_sparse(fa, m, n, nzA, A);
    else
        read_dense(fa, m, n, A);
    if (nzB > 0)
        read_sparse(fb, n, p, nzB, B);
    else
        read_dense(fb, n, p, B);
    fclose(fa);
    fclose(fb);
#endif

    C = (float *)calloc(m * p, sizeof(float));
    if (C == NULL)
    {
        printf("Out of memory C1! \n");
        exit(1);
    }
    // C2 = (float *)calloc(N*P,sizeof(float));
    // if (C2==NULL) {printf("Out of memory C2! \n"); exit(1);}

    // naive implementation
#ifdef TIMING
    double flops = get_operation_count(m, n, p);
    gettimeofday(&before, NULL);
#endif

    for (r = 0; r < REP; r++)
        matrix_mult(m, n, p, A, B, C);

        // for (r = 0; r < REP; r++)
        //     cuda_do_compute(m, n, p, A, B, C);

#ifdef TIMING
    gettimeofday(&after, NULL);
    double avg_execution_time = ((after.tv_sec + (after.tv_usec / 1000000.0)) -
                                 (before.tv_sec + (before.tv_usec / 1000000.0))) /
                                REP;
    printf("GFLOP/s: %f\n", (flops / (avg_execution_time)) / 1e9);
    printf("Seconds: %f\n", avg_execution_time);

#endif

#ifdef GENERATE
    if ((fc = fopen("gen_result.mtx", "wt")) == NULL)
        exit(3);
#else
    if ((fc = fopen(argv[3], "wt")) == NULL)
        exit(3);
#endif
    write_sparse(fc, m, p, C);
    fclose(fc);
    free(A);
    free(B);
    free(C);
    // free(C2);
}
