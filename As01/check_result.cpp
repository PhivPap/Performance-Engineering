#include <cassert>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "MatrixParser/matrixparser.h"

void check_result(char* file_expected, char* file_actual, int m, int n, int nz)
{
    printf("test passed1\n");
    fflush(stdout);
    FILE *fa, *fb;
    
    if ((fa = fopen(file_expected, "rt")) == NULL)
            exit(1);
    printf("test passed2\n");
    fflush(stdout);
    if ((fb = fopen(file_actual, "rt")) == NULL)
        exit(2);
    printf("test passed3\n");
    fflush(stdout);
    int size = m * n;
    float *expected = (float *)calloc(size, sizeof(float));
    float *actual = (float *)calloc(size, sizeof(float));
    printf("test passed4\n");
    fflush(stdout);
    
    read_sparse(fa, m, n, nz, expected);
    printf("test passed4.5\n");
    fflush(stdout);
    read_sparse(fb, m, n, nz, actual);
    printf("test passed5\n");
    fflush(stdout);

    for (int i = 0; i < size; ++i)
    {
        assert(("There's an error on index " + std::to_string(i), expected[i] == actual[i]));
    }
    printf("test passed6\n");
    fflush(stdout);
}