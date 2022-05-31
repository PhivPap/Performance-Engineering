#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>

const int THREAD_COUNT = 16;
const int SCALAR = 0.99;

void generate_mat(int m, int n, float *A) {
  int i;

  for (i=0; i<(m*n); i++) 
	  A[i] = i * SCALAR;
}
     

int main(int argc, char **argv)
{
  omp_set_num_threads(THREAD_COUNT);
  
  float *A;
  float B, val;
  int y, x;
  int m = 1000;
  A = (float *)calloc(m * m, sizeof(float));
  B = 0;
  
  generate_mat(m, m, A);

  struct timeval before, after;
  gettimeofday(&before, NULL); 
  
  #pragma omp parallel \
  shared(m, A, B) \
  private(y, x, val)
  {
#ifndef OVERHEAD
    float local;
#endif
    #pragma omp for \
    schedule(static)
    for (y = 0; y < m; y++)
        for (x = 0; x < m; x++) {
        val = A[x + y * m];
#ifdef OVERHEAD
        #pragma omp critical(B)
        {
            B += val;
        }
#else
        local += val;
#endif
      }
#ifndef OVERHEAD
    #pragma omp critical(B)
    {
      B += local;
    }
#endif
  }
  
  gettimeofday(&after, NULL);
  printf("Execution time: %10.6f seconds \n", ((after.tv_sec + (after.tv_usec / 1000000.0)) -
            (before.tv_sec + (before.tv_usec / 1000000.0))));
}
