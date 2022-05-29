#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

const float SCALAR = 7.777;

const int THREAD_COUNT = 4;

void generate_mat(int m, int n, float *A) {
  int i;

  for (i=0; i<(m*n); i++) 
	  A[i] = i;
}
     

int main(int argc, char **argv)
{
  omp_set_num_threads(THREAD_COUNT);
  
  float *A;
  int y, x;
  int m = 1000;
  A = (float *)calloc(m * m, sizeof(float));
  
  generate_mat(m, m, A);

#ifdef FALSESHARING
  #pragma omp parallel for \
  shared(m, A) \
  private(y, x) \
  schedule(static, 1)
#else
  #pragma omp parallel for \
  shared(m, A) \
  private(y, x) \
  schedule(static)
#endif
  for (y = 0; y < m; y++)
    for (x = 0; x < m; x++) 
      A[x + y * x] *= SCALAR;
}
