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


#include "md.h"
#include "blas.h"

void QMD_saxpy(int n, REAL alpha, REAL * x, int incx, REAL * y, int incy)
{
    saxpy(&n, &alpha, x, &incx, y, &incy);
}

void QMD_sscal(int n, REAL alpha, REAL * x, int incx)
{
    sscal(&n, &alpha, x, &incx);
}

void QMD_scopy(int n, REAL * x, int incx, REAL * y, int incy)
{
    scopy(&n, x, &incx, y, &incy);
}

REAL QMD_sdot(int n, REAL * x, int incx, REAL * y, int incy)
{
    return sdot(&n, x, &incx, y, &incy);
}

/******/
