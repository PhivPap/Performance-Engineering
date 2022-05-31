#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>

const float SCALAR = 7.777;

const int THREAD_COUNT = 16;

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
  int m = 10000;
  A = (float *)calloc(m * m, sizeof(float));
  
  generate_mat(m, m, A);

  struct timeval before, after;
  gettimeofday(&before, NULL); 

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
  for (x = 0; x < m; x++)
    for (y = 0; y < m; y++) 
      A[x + y * m] *= SCALAR;
  
  gettimeofday(&after, NULL);
  printf("Execution time: %10.6f seconds \n", ((after.tv_sec + (after.tv_usec / 1000000.0)) -
            (before.tv_sec + (before.tv_usec / 1000000.0))));
}
