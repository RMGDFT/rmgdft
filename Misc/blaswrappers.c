/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/


#include "main.h"

void QMD_daxpy (int n, double alpha, double * x, int incx, double * y, int incy)
{
//    daxpy (&n, &alpha, x, &incx, y, &incy);
    int i, iy=0, ix=0;

    if((incx == 1) && (incy == 1)) {

       for(i = 0;i < n;i++) {
           y[i] = alpha * x[i] + y[i];
       }
       return;
    }

    for(i = 0;i < n;i++) {
        y[iy] = alpha * x[ix] + y[iy];
        ix += incx;
        iy += incy;
    }

}

void QMD_dscal (int n, double alpha, double * x, int incx)
{
//    sscal (&n, &alpha, x, &incx);
      int i, ix=0;
      if(incx == 1) {
          for(i = 0;i < n;i++) {
              x[i] = alpha * x[i];
          }
          return;
      }
      for(i = 0;i < n;i++) {
          x[ix] = alpha * x[ix];
          ix += incx;
      }
}

void QMD_dcopy (int n, double * x, int incx, double * y, int incy)
{
//    scopy (&n, x, &incx, y, &incy);
      int i, ix=0, iy=0;
      if((incx == 1) && (incy == 1)) {
         for(i = 0;i < n;i++) {
             y[i] = x[i];
         }
         return;
      }

      for(i = 0;i < n;i++) {
          y[iy] = x[ix];
          ix += incx;
          iy += incy;
      }
}

double QMD_ddot (int n, double * x, int incx, double * y, int incy)
{
//    return sdot (&n, x, &incx, y, &incy);
      int i, ix = 0, iy = 0;
      double stemp = 0.0;

      if((incx == 1) && (incy == 1)) {
          for(i = 0;i < n;i++) {
              stemp += y[i] * x[i];
          }
          return stemp;
      }
      for(i = 0;i < n;i++) {
          stemp += y[iy] * x[ix];
          ix += incx;
          iy += incy;
      }
      return stemp;

}

void QMD_dswap(int n, double * x, int incx, double *y, int incy)
{
      int i, ix=0, iy=0;
      double temp;

      for(i = 0;i < n;i++) {
          temp = y[iy];
          y[iy] = x[ix];
          x[ix] = temp;
          ix += incx;
          iy += incy;
      }


}

void my_copy(double *in, double *out, int length) {
    int ione = 1;
    QMD_dcopy(length, in, ione, out, ione);
}
void my_scal(double alpha, double *vect, int length) {
    int ione = 1;
    QMD_dscal(length, alpha, vect, ione); 
}
void my_axpy(double alpha, double *in, double *out, int length) {
    int ione = 1;
    QMD_daxpy(length, alpha, in, ione, out, ione);
}
void my_swap(double *vec1, double *vec2, int length) {
    int ione =1;
    QMD_dswap(length, vec1, ione, vec2, ione);
}

/******/
void QMD_saxpy (int n, float alpha, float * x, int incx, float * y, int incy)
{
//    daxpy (&n, &alpha, x, &incx, y, &incy);
    int i, iy=0, ix=0;

    if((incx == 1) && (incy == 1)) {

       for(i = 0;i < n;i++) {
           y[i] = alpha * x[i] + y[i];
       }
       return;
    }

    for(i = 0;i < n;i++) {
        y[iy] = alpha * x[ix] + y[iy];
        ix += incx;
        iy += incy;
    }

}

void QMD_sscal (int n, float alpha, float * x, int incx)
{
//    sscal (&n, &alpha, x, &incx);
      int i, ix=0;
      if(incx == 1) {
          for(i = 0;i < n;i++) {
              x[i] = alpha * x[i];
          }
          return;
      }
      for(i = 0;i < n;i++) {
          x[ix] = alpha * x[ix];
          ix += incx;
      }
}

void QMD_scopy (int n, float * x, int incx, float * y, int incy)
{
//    scopy (&n, x, &incx, y, &incy);
      int i, ix=0, iy=0;
      if((incx == 1) && (incy == 1)) {
         for(i = 0;i < n;i++) {
             y[i] = x[i];
         }
         return;
      }

      for(i = 0;i < n;i++) {
          y[iy] = x[ix];
          ix += incx;
          iy += incy;
      }
}

float QMD_sdot (int n, float * x, int incx, float * y, int incy)
{
      int i, ix = 0, iy = 0;
      double stemp = 0.0;

      if((incx == 1) && (incy == 1)) {
          for(i = 0;i < n;i++) {
              stemp += y[i] * x[i];
          }
          return stemp;
      }
      for(i = 0;i < n;i++) {
          stemp += (double)y[iy] * (double)x[ix];
          ix += incx;
          iy += incy;
      }
      return (float)stemp;

}

void QMD_sswap(int n, float * x, int incx, float *y, int incy)
{
      int i, ix=0, iy=0;
      float temp;

      for(i = 0;i < n;i++) {
          temp = y[iy];
          y[iy] = x[ix];
          x[ix] = temp;
          ix += incx;
          iy += incy;
      }


}
