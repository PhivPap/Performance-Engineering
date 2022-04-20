#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <assert.h>
#include "mmio.h"
#include "MatrixParser/matrixparser.h"

#include <immintrin.h>

#define N 512
#define M 512
#define P 512

#define REP 5

enum Version { OPT, OPT_SIMD, LOOP_TILING};
enum Version version = OPT_SIMD;

void matrix_mult_optimised(int m, int n, int p, float *A, float *B, float *C){
    int i, j, k, in, jn;
    float sum;
    for (i = 0; i < m; i++){
        int ip = i * p;
        for (j = 0; j < p; j++){
            sum = 0;
            in = i * n;
            jn = j * n;

            for (k = 0; k < n; k++)
                sum += A[in + k] * B[jn + k];
                
            C[ip + j] = sum;
        }
    }
}


void matrix_mult_optimised_simd(int m, int n, int p, float *A, float *B, float *C)
{
    int REG_SIZE = 16; 
    int i, j, k, in, jn;
    float sum;
    float buffer[16];
    __m512 reg_a, reg_b, reg_res, reg_temp;

    int first_remaining = n - (n % REG_SIZE);

    for (i = 0; i < m; i++)
    {
        int ip = i * p;
        for (j = 0; j < p; j++)
        {
            sum = 0;
            reg_res = _mm512_set1_ps(0.0);
            in = i * n;
            jn = j * n;

            for (k = 0; k < n / REG_SIZE; k++)
            {
                // Assumes A is stored row major while B is stored column major
                reg_a = _mm512_loadu_ps(&A[in + k * REG_SIZE]);
                reg_b = _mm512_loadu_ps(&B[jn + k * REG_SIZE]);
                reg_temp = _mm512_mul_ps(reg_a, reg_b);
                reg_res = _mm512_add_ps(reg_res, reg_temp);
            }

            //remaining part
            for(int l = first_remaining; l < n; ++l){
                sum += A[in + l] * B[jn + l];
            }

            _mm512_store_ps(buffer, reg_res);
            sum += buffer[0] + buffer[1] + buffer[2] + buffer[3] + buffer[4] + buffer[5] + buffer[6] + buffer[7] +
                buffer[8] + buffer[9] + buffer[10] + buffer[11] + buffer[12] + buffer[13] + buffer[14] + buffer[15];

            C[ip + j] = sum;
        }
    }
}

// m = rows of C and rows of A
// p = colums of C and colums of B
// n = colums of A and rows of B
// MATRIX = [ROWS][COLUMS]
// A = [m][n]
// B = [n][p]
// C = [m][p]
void matrix_mult_optimised_loop_tiling(int m, int n, int p, float *A, float *B, float *C)
{
    int i, j, k, ii, kk;
    int ib = 2, kb = 2;
    float acc00, acc01, acc10, acc11;

    for (ii = 0; ii < m; ii += ib)
    {
        for (kk = 0; kk < p; kk += kb)
        {
            for (j=0; j < n; j += 2)
            {
                for (i = ii; i < ii + ib; i += 2)
                {
                    if (kk == 0)
                        acc00 = acc01 = acc10 = acc11 = 0;
                    else
                    {
                        acc00 = C[i * p + j];
                        acc01 = C[i * p + j + 1];
                        acc10 = C[(i + 1) * p + j];
                        acc11 = C[(i + 1) * p + j + 1];
                    }
                    for (k = kk; k < kk + kb; k++)
                    {
                        acc00 += B[k * p + j]       * A[i * n + k];
                        acc01 += B[k * p + j + 1]   * A[i * n + k];
                        acc10 += B[k * p + j]       * A[(i + 1) * n + k];
                        acc11 += B[k * p + j + 1]   * A[(i + 1) * n + k];
                    }
                    C[i * p + j] = acc00;
                    C[i * p + j + 1] = acc01;
                    C[(i + 1) * p + j] = acc10;
                    C[(i + 1) * p + j + 1] = acc11;
                }
            }
        }
    }
}

double get_operation_count(int m, int n, int p)
{
    switch (version){
        case OPT:           return (double)m * n * p * 2;
        case OPT_SIMD:      return (double)m * n * p * 2;
        case LOOP_TILING:   return (double)m * n * p * 2;
        default:            assert(0); 
    }
    return 0;
}

void to_column_major(int n, int p, float *B)
{
    float *Btmp = (float *)calloc(n * p, sizeof(float));
    if (Btmp == NULL)
    {
        printf("Out of memory B! \n");
        exit(1);
    }

    for (int row = 0; row < n; row++)
        for (int column = 0; column < p; column++)
            Btmp[column * n + row] = B[row * n + column];

    for (int i = 0; i < n * p; i++)
        B[i] = Btmp[i];
}

int main(int argc, char **argv)
{
    float *A, *B, *B_tmp, *C;
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
    B_tmp = (float *)calloc(n * p, sizeof(float));
    if (B_tmp == NULL)
    {
        printf("Out of memory B! \n");
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
        read_sparse(fb, n, p, nzB, B_tmp);
    else
        read_dense(fb, n, p, B_tmp);
    fclose(fa);
    fclose(fb);
#endif

    switch (version){
        case OPT:
        case OPT_SIMD: {
            transpose(n, p, B_tmp, B);
            break;
        }
        default: assert(0);
    }


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

    switch (version){
        case OPT: {
            for (r = 0; r < REP; r++)
                matrix_mult_optimised(m, n, p, A, B, C);
            break;
        }

        case OPT_SIMD: {
            for (r = 0; r < REP; r++)
                matrix_mult_optimised_simd(m, n, p, A, B, C);
            break;
        }

        case LOOP_TILING: {
            for (r = 0; r < REP; r++)
                matrix_mult_optimised_loop_tiling(m, n, p, A, B_tmp, C);
            break;
        }
    
        default: assert(0);
    }

#ifdef TIMING
    gettimeofday(&after, NULL);
    double avg_execution_time = ((after.tv_sec + (after.tv_usec / 1000000.0)) -
                                 (before.tv_sec + (before.tv_usec / 1000000.0))) /
                                REP;
    printf("GFLOP/s: %f\n", flops / (1e9 * avg_execution_time));
    printf("Seconds: %f\n", avg_execution_time);

#endif

#ifdef GENERATE
    if ((fc = fopen("gen_result.mtx", "wt")) == NULL)
        exit(3);
#else
    if ((fc = fopen(argv[3], "wt")) == NULL)
        exit(3);


#endif
    int nz = write_sparse(fc, m, p, C);
    fclose(fc);
    free(A);
    free(B);
    free(C);
}