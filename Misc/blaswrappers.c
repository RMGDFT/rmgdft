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

void QMD_saxpy (int n, REAL alpha, REAL * x, int incx, REAL * y, int incy)
{
//    saxpy (&n, &alpha, x, &incx, y, &incy);
    int i, iy=0, ix=0;
    for(i = 0;i < n;i++) {
        y[iy] = alpha * x[ix] + y[iy];
        ix += incx;
        iy += incy;
    }

}

void QMD_sscal (int n, REAL alpha, REAL * x, int incx)
{
//    sscal (&n, &alpha, x, &incx);
      int i, ix=0;
      for(i = 0;i < n;i++) {
          x[ix] = alpha * x[ix];
          ix += incx;
      }
}

void QMD_scopy (int n, REAL * x, int incx, REAL * y, int incy)
{
//    scopy (&n, x, &incx, y, &incy);
      int i, ix=0, iy=0;
      for(i = 0;i < n;i++) {
          y[iy] = x[ix];
          ix += incx;
          iy += incy;
      }
}

REAL QMD_sdot (int n, REAL * x, int incx, REAL * y, int incy)
{
//    return sdot (&n, x, &incx, y, &incy);
      int i, ix = 0, iy = 0;
      REAL stemp = 0.0;
      for(i = 0;i < n;i++) {
          stemp += y[iy] * x[ix];
          ix += incx;
          iy += incy;
      }
      return stemp;

}

void QMD_sswap(int n, REAL * x, int incx, REAL *y, int incy)
{
      int i, ix=0, iy=0;
      REAL temp;

      for(i = 0;i < n;i++) {
          temp = y[iy];
          y[iy] = x[ix];
          x[ix] = temp;
          ix += incx;
          iy += incy;
      }


}

void my_copy(REAL *in, REAL *out, int length) {
    int ione = 1;
    QMD_scopy(length, in, ione, out, ione);
}
void my_scal(REAL alpha, REAL *vect, int length) {
    int ione = 1;
    QMD_sscal(length, alpha, vect, ione); 
}
void my_axpy(REAL alpha, REAL *in, REAL *out, int length) {
    int ione = 1;
    QMD_saxpy(length, alpha, in, ione, out, ione);
}
void my_swap(REAL *vec1, REAL *vec2, int length) {
    int ione;
    QMD_sswap(length, vec1, ione, vec2, ione);
}

/******/
