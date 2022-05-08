#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "mmio.h"
#include <omp.h>
#include <immintrin.h>
#include <assert.h>
#include "MatrixParser/matrixparser.h"

#define N 512
#define M 512
#define P 512

#define REP 5

const int THREAD_COUNT = 4;

double matrix_mult_omp_simd(int m, int n, int p, float *A, float *B, float *C)
{
    int REG_SIZE = 16; 
    int i, j, k, in, jn;
    float sum;
    float buffers[THREAD_COUNT][16];
    __m256 reg_a, reg_b, reg_res, reg_temp;

    int first_remaining = n - (n % REG_SIZE);

    double elapsed = 0.0;

    #pragma omp parallel \
    reduction(max: elapsed) \
    default(none) \
    shared(m, n, p, A, B, C, REG_SIZE, buffers, first_remaining) \
    private(sum, j, k, in, jn, reg_a, reg_b, reg_res, reg_temp)
    {
        float *buffer = buffers[omp_get_thread_num()];
        struct timeval before, after;
        gettimeofday(&before, NULL);

        #pragma omp for  
        for (i = 0; i < m; i++)
        {
            int ip = i * p;
            for (j = 0; j < p; j++)
            {
                sum = 0;
                reg_res = _mm256_set1_ps(0.0);
                in = i * n;
                jn = j * n;

                for (k = 0; k < n / REG_SIZE; k++)
                {
                    // Assumes A is stored row major while B is stored column major
                    reg_a = _mm256_loadu_ps(&A[in + k * REG_SIZE]);
                    reg_b = _mm256_loadu_ps(&B[jn + k * REG_SIZE]);
                    reg_res = _mm256_fmadd_ps(reg_a, reg_b, reg_res);
                }

                //remaining part
                for(int l = first_remaining; l < n; ++l){
                    sum += A[in + l] * B[jn + l];
                }

                _mm256_store_ps(buffer, reg_res);
                sum += buffer[0] + buffer[1] + buffer[2] + buffer[3] + buffer[4] + buffer[5] + buffer[6] + buffer[7] +
                    buffer[8] + buffer[9] + buffer[10] + buffer[11] + buffer[12] + buffer[13] + buffer[14] + buffer[15];

                C[ip + j] = sum;
            }
        }

        gettimeofday(&after, NULL);
        elapsed = ((after.tv_sec + (after.tv_usec / 1000000.0)) -
                                    (before.tv_sec + (before.tv_usec / 1000000.0)));
    }
    return elapsed;
    
}

double get_operation_count(int m, int n, int p)
{
    return (double)m * n * p * 2;
}

int main(int argc, char **argv)
{
    omp_set_num_threads(THREAD_COUNT);
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

    transpose(n, p, B_tmp, B);

    C = (float *)calloc(m * p, sizeof(float));
    if (C == NULL)
    {
        printf("Out of memory C1! \n");
        exit(1);
    }


    double avg_execution_time = 0.0;
    for (r = 0; r < REP; r++)
        avg_execution_time += matrix_mult_omp_simd(m, n, p, A, B, C);
    avg_execution_time /= REP;

    
    printf("GFLOP/s: %f\n", get_operation_count(m, n, p) / (1e9 * avg_execution_time));
    printf("Seconds: %f\n", avg_execution_time);


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
