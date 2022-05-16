#include <stdlib.h>
// #include <sys/time.h>
#include "matrixparser.h"
#include "../mmio.h"

#define N 512
#define M 512
#define P 512

void transpose(int m, int n, float *A, float *B) {
    int i, j;
    for(i=0; i<m; i++) {
        for(j=0; j<n; j++) {
            B[i+j*m] = A[i*n+j];
        }
    }
}

void generate_mat(int m, int n, int p, float *A, float *B)
{
    int i;

    for (i = 0; i < (m * n); i++)
        A[i] = 1; // i/10;
    for (i = 0; i < (n * p); i++)
        B[i] = 1; // i/5;
}

void read_sparse(FILE *f, int m, int n, int nz, float *A)
{
    int i, row, col;
    float val;

    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    for (i = 0; i < nz; i++)
    {
        fscanf(f, "%d %d %f\n", &row, &col, &val);
        A[(row - 1) * n + col - 1] = val; /* adjust from 1-based to 0-based */
    }
}

int write_sparse(FILE *f, int m, int p, const float *C)
{
    int i, nz = 0;
    MM_typecode matcode;

    for (i = 0; i < m * p; i++)
        if (C[i] != 0.0)
            nz++;

    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    mm_set_coordinate(&matcode);
    mm_set_real(&matcode);

    mm_write_banner(f, matcode);
    mm_write_mtx_crd_size(f, m, p, nz);

    for (i = 0; i < m * p; i++)
    {
        if (C[i] != 0.0)
            fprintf(f, "%d %d %f\n", i / p + 1, i % p + 1, C[i]);
    }

    return nz;
}

void read_dense(FILE *f, int m, int n, float *A)
{
    int row, col;

    for (row = 0; row < m; row++)
    {
        for (col = 0; col < n; col++)
        {
            fscanf(f, "%f ", &A[row * n + col]);
            //	printf("%20.19f \n", A[row*(*n)+col]);
        }
    }
}

int read_mat(int *m, int *n, int *p, int *nzA, int *nzB, FILE *fa, FILE *fb)
{
    MM_typecode ta, tb;
    int ret_code;
    int n1;

    if (mm_read_banner(fa, &ta) != 0)
    {
        printf("Could not process Matrix Market banneri for A.\n");
        return -3;
    }
    if (mm_read_banner(fb, &tb) != 0)
    {
        printf("Could not process Matrix Market banner for B.\n");
        return -4;
    }

    if (mm_is_complex(ta))
        return -6;
    if (mm_is_complex(tb))
        return -7;

    if (mm_is_matrix(ta) && mm_is_sparse(ta))
    {
        if ((ret_code = mm_read_mtx_crd_size(fa, m, n, nzA)) != 0)
            return -10;
    }
    else if (mm_is_matrix(ta) && mm_is_array(ta))
    {
        *nzA = 0;
        if ((ret_code = mm_read_mtx_array_size(fa, m, n)) != 0)
            return -11;
    }
    else
        return -8;

    if (mm_is_matrix(tb) && mm_is_sparse(tb))
    {
        if ((ret_code = mm_read_mtx_crd_size(fb, &n1, p, nzB)) != 0)
            return -10;
    }
    else if (mm_is_matrix(tb) && mm_is_array(tb))
    {
        *nzB = 0;
        if ((ret_code = mm_read_mtx_array_size(fb, &n1, p)) != 0)
            return -11;
    }
    else
        return -9;

    if (*n != n1)
        return -15;

    return 0;
    /* find out size of sparse matrix .... */
}