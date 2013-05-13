/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/blaswrappers.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5
 * COPYRIGHT
 *   Copyright (C) 1995  Emil Briggs
 *   Copyright (C) 1998  Emil Briggs, Charles Brabec, Mark Wensell, 
 *                       Dan Sullivan, Chris Rapcewicz, Jerzy Bernholc
 *   Copyright (C) 2001  Emil Briggs, Wenchang Lu,
 *                       Marco Buongiorno Nardelli,Charles Brabec, 
 *                       Mark Wensell,Dan Sullivan, Chris Rapcewicz,
 *                       Jerzy Bernholc
 * FUNCTION
 *   Wrapper routines for blas calls.
 *   saxpy, scopy, sscal, sdot 
 * SOURCE
 */


#include "main.h"

void QMD_daxpy (int n, rmg_double_t alpha, rmg_double_t * x, int incx, rmg_double_t * y, int incy)
{
//    saxpy (&n, &alpha, x, &incx, y, &incy);
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

void QMD_dscal (int n, rmg_double_t alpha, rmg_double_t * x, int incx)
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

void QMD_dcopy (int n, rmg_double_t * x, int incx, rmg_double_t * y, int incy)
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

rmg_double_t QMD_ddot (int n, rmg_double_t * x, int incx, rmg_double_t * y, int incy)
{
//    return sdot (&n, x, &incx, y, &incy);
      int i, ix = 0, iy = 0;
      rmg_double_t stemp = 0.0;

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

void QMD_dswap(int n, rmg_double_t * x, int incx, rmg_double_t *y, int incy)
{
      int i, ix=0, iy=0;
      rmg_double_t temp;

      for(i = 0;i < n;i++) {
          temp = y[iy];
          y[iy] = x[ix];
          x[ix] = temp;
          ix += incx;
          iy += incy;
      }


}

void my_copy(rmg_double_t *in, rmg_double_t *out, int length) {
    int ione = 1;
    QMD_dcopy(length, in, ione, out, ione);
}
void my_scal(rmg_double_t alpha, rmg_double_t *vect, int length) {
    int ione = 1;
    QMD_dscal(length, alpha, vect, ione); 
}
void my_axpy(rmg_double_t alpha, rmg_double_t *in, rmg_double_t *out, int length) {
    int ione = 1;
    QMD_daxpy(length, alpha, in, ione, out, ione);
}
void my_swap(rmg_double_t *vec1, rmg_double_t *vec2, int length) {
    int ione =1;
    QMD_dswap(length, vec1, ione, vec2, ione);
}

/******/
void QMD_saxpy (int n, rmg_float_t alpha, rmg_float_t * x, int incx, rmg_float_t * y, int incy)
{
//    saxpy (&n, &alpha, x, &incx, y, &incy);
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

void QMD_sscal (int n, rmg_float_t alpha, rmg_float_t * x, int incx)
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

void QMD_scopy (int n, rmg_float_t * x, int incx, rmg_float_t * y, int incy)
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

rmg_float_t QMD_sdot (int n, rmg_float_t * x, int incx, rmg_float_t * y, int incy)
{
      int i, ix = 0, iy = 0;
      rmg_double_t stemp = 0.0;

      if((incx == 1) && (incy == 1)) {
          for(i = 0;i < n;i++) {
              stemp += y[i] * x[i];
          }
          return stemp;
      }
      for(i = 0;i < n;i++) {
          stemp += (rmg_double_t)y[iy] * (rmg_double_t)x[ix];
          ix += incx;
          iy += incy;
      }
      return (rmg_float_t)stemp;

}

void QMD_sswap(int n, rmg_float_t * x, int incx, rmg_float_t *y, int incy)
{
      int i, ix=0, iy=0;
      rmg_float_t temp;

      for(i = 0;i < n;i++) {
          temp = y[iy];
          y[iy] = x[ix];
          x[ix] = temp;
          ix += incx;
          iy += incy;
      }


}
