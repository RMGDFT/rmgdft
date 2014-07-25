#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"


void qgemv(int n, double alpha, double *A, double *X, double *Y)
{

  int ix, iy;
  __float128 *qA, *qX, *qY, *qYY, qalpha, temp;

  
  my_malloc(qA, 2*n*n, __float128);
  my_malloc(qX, 2*n, __float128);
  my_malloc(qY, 2*n, __float128);
  my_malloc(qYY, 2*n, __float128);

  qalpha = (__float128)alpha;

  for(ix = 0;ix < n*n;ix++) {
      qA[ix] = (__float128)A[ix];
  }
  for(ix = 0;ix < n;ix++) {
      qX[ix] = (__float128)X[ix];
  }
  for(ix = 0;ix < n;ix++) {
      qY[ix] = (__float128)0.0;
      qYY[ix] = (__float128)0.0;
  }

  for(ix = 0;ix < n;ix++) {
      temp = qalpha * qX[ix];
      for(iy = 0;iy < n;iy++) {
          qY[iy] += temp*qA[ix*n + iy];
      }
  }  
  for(ix = 0;ix < n;ix++) {
      temp = qalpha * qY[ix];
      for(iy = 0;iy < n;iy++) {
          qYY[iy] += temp*qA[ix*n + iy];
      }
  }  

  for(ix = 0;ix < n;ix++) {
      Y[ix] = (double)qYY[ix];
  }

  my_free(qYY);
  my_free(qY);
  my_free(qX);
  my_free(qA);
}
