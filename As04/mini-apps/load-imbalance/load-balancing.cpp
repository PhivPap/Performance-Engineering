#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

const int THREAD_COUNT = 4;

void generate_mat(int m, int n, float *A) {
  int rows, cols, r, c;

  for(r = 0; r < m; r++) 
    for(c = 0; c < m; c++)
      if(r >= c)
        A[r + c * m] = r + c;
}
     

int main(int argc, char **argv)
{
  omp_set_num_threads(THREAD_COUNT);
  
  float *A, *B;
  int y, x;
  int m = 1000;
  A = (float *)calloc(m * m, sizeof(float));
  B = (float *)calloc(m, sizeof(float));

  generate_mat(m, m, A);

#ifdef IMBALANCE
  #pragma omp parallel for \
  shared(m, A, B) \
  private(y, x) \
  schedule(static)
#else
  #pragma omp parallel for \
  shared(m, A, B) \
  private(y, x) \
  schedule(static, 1)
#endif
  for (x = 0; x < m; x++)
    for (y = 0; y < m; y++) 
      B[y] += A[x + y * m];
}
