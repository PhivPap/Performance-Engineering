#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>

const int THREAD_COUNT = 16;

void generate_mat(int m, int n, float *A) {
  int rows, cols, r, c;

  for(r = 0; r < m; r++) 
  {
    int column_max;
    if (r < 15000)
      column_max = 100;
    else
      column_max = m;
    for(c = 0; c < column_max; c++)
      A[r + c * m] = (r + c) * 0.95;
  }
}
     

int main(int argc, char **argv)
{
  omp_set_num_threads(THREAD_COUNT);
  
  float *A, *B;
  int r, c;
  int m = 16000;
  A = (float *)calloc(m * m, sizeof(float));
  B = (float *)calloc(m, sizeof(float));

  generate_mat(m, m, A);

  struct timeval before, after;
  gettimeofday(&before, NULL); 

#ifdef IMBALANCE
  #pragma omp parallel for \
  shared(m, A, B) \
  private(r, c) \
  schedule(static)
#else
  #pragma omp parallel for \
  shared(m, A, B) \
  private(r, c) \
  schedule(static, 1)
#endif
  for(r = 0; r < m; r++) 
  {
    int column_max;
    if (r < 15000)
      column_max = 100;
    else
      column_max = m;
    for(c = 0; c < column_max; c++)
      B[r] += A[c + r * m];
  }
  
  gettimeofday(&after, NULL);
  printf("Execution time: %10.6f seconds \n", ((after.tv_sec + (after.tv_usec / 1000000.0)) -
            (before.tv_sec + (before.tv_usec / 1000000.0))));
}
