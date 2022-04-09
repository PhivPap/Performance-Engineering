#include <iostream>
#include <chrono>
#include <cuda.h>
#include "MatrixParser/matrixparser.h"
#include "mmio.h"

#define check(error) checkCudaCall(error, __LINE__)

const uint32_t REP = 10;
const uint32_t BLOCK_WIDTH = 16;
const uint32_t BLOCK_HEIGHT = 16;

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

__global__ void matrix_mul(float *A, float *B, float *C, int n, int C_width, int C_height)
{
    const uint32_t x = blockDim.x * blockIdx.x + threadIdx.x;
    const uint32_t y = blockDim.y * blockIdx.y + threadIdx.y;

    if (x >= C_width || y >= C_height)
        return;

    const uint32_t xn = x * n;
    float sum = 0;

    for (int i = 0; i < n; i++)
        sum += A[xn + i] * B[y + C_width * i];

    C[x * C_width + y] = sum;
}

void gpu_memory_init(int m, int n, int p,
                     float *A, float *B, float *C,
                     float **d_A, float **d_B, float **d_C)
{
    const uint32_t A_size = sizeof(float) * n * m;
    const uint32_t B_size = sizeof(float) * p * n;
    const uint32_t C_size = sizeof(float) * p * m; 

    check(cudaMalloc(d_A, A_size));
    check(cudaMalloc(d_B, B_size));
    check(cudaMalloc(d_C, C_size));

    check(cudaMemcpy(*d_A, A, A_size, cudaMemcpyHostToDevice));
    check(cudaMemcpy(*d_B, B, B_size, cudaMemcpyHostToDevice));
    check(cudaMemcpy(*d_C, C, C_size, cudaMemcpyHostToDevice));
}

void gpu_memory_free(float *d_A, float *d_B, float *d_C){
    check(cudaFree(d_A));
    check(cudaFree(d_B));
    check(cudaFree(d_C));
}

double get_operation_count(int m, int n, int p){
    return (double)m * n * p * 9;
}

void cuda_do_compute(int m, int n, int p, float *A, float *B, float *C)
{
    float *d_A, *d_B, *d_C;
    gpu_memory_init(m, n, p, A, B, C, &d_A, &d_B, &d_C);
    const uint32_t grid_height = m % BLOCK_HEIGHT == 0 ? m / BLOCK_HEIGHT : (m / BLOCK_HEIGHT) + 1;
    const uint32_t grid_width = p % BLOCK_WIDTH == 0 ? p / BLOCK_WIDTH : (p / BLOCK_WIDTH) + 1;
    dim3 grid(grid_width, grid_height);
    dim3 block(BLOCK_WIDTH, BLOCK_HEIGHT);

    const auto start = std::chrono::system_clock::now();
    matrix_mul<<<grid, block>>>(d_A, d_B, d_C, n, p, m);
    cudaDeviceSynchronize();
    const auto end = std::chrono::system_clock::now();

    const double avg_execution_time = (double)(end - start).count() / (double)10e9;
    printf("Reference code: %10.2f GFLOP/s \n", get_operation_count(m, n, p) / (avg_execution_time * 10e9));
    printf("Reference code: %10.2f seconds \n", avg_execution_time);
    cudaMemcpy(C, d_C, sizeof(float) * p * m, cudaMemcpyDeviceToHost);
    gpu_memory_free(d_A, d_B, d_C);
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
        cuda_do_compute(m, n, p, A, B, C);

#ifdef TIMING
    gettimeofday(&after, NULL);
    double avg_execution_time = ((after.tv_sec + (after.tv_usec / 1000000.0)) -
                                 (before.tv_sec + (before.tv_usec / 1000000.0))) /
                                REP;
    printf("Reference code: %10.2f GFLOP/s \n", flops / (10e9 * avg_execution_time));
    printf("Reference code: %10.2f seconds \n", avg_execution_time);

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