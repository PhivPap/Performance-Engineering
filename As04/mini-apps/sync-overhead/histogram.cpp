#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

const int THREAD_COUNT = 4;
const int COLORS = 8;

void generate_mat(int m, int n, float *A) {
  int i;

  for (i=0; i<(m*n); i++) 
	  A[i] = i % COLORS;
}
     

int main(int argc, char **argv)
{
  omp_set_num_threads(THREAD_COUNT);
  
  float *A, *B;
  int y, x, color;
  int m = 1000;
  A = (float *)calloc(m * m, sizeof(float));
  B = (float *)calloc(COLORS, sizeof(float));
  
  generate_mat(m, m, A);
  
  #pragma omp parallel \
  shared(m, A, B) \
  private(y, x, color)
  {
#ifndef OVERHEAD
    float *local;
    local = (float *)calloc(COLORS, sizeof(float));
#endif
    #pragma omp for \
    schedule(static)
    for (y = 0; y < m; y++)
        for (x = 0; x < m; x++) {
        color = A[x + y * x];
#ifdef OVERHEAD
        #pragma omp critical(B)
        {
            B[color]++;
        }
#else
        local[color]++;
#endif
      }
#ifndef OVERHEAD
    #pragma omp critical(B)
    {
      for (int i = 0; i < COLORS; i++)
        B[i] += local[i];
    }
#endif
  }
}
