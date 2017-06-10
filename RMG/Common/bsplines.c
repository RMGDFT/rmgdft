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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "grid.h"
#include "common_prototypes.h"
#include "main.h"

static void bsplvb (double *t, int jhigh, int index, double x, int left, double *biatx);
static void banfac (double *w, int nroww, int nrow, int nbandl, int nbandu, int *iflag);
static void banslv (double *w, int nroww, int nrow, int nbandl, int nbandu, double *b);
static void huntn (double *xx, int n, int kord, double x, int *jlo);
static void dbsnak (int nx, double *xvec, int kxord, double *xknot);
static double dbsdca (int iderx, double x, int kx, double *xknot, double *bcoef, int leftx);


static void spli3d (double *xyzvec, int ldf, int mdf, int zdf,
             double *xyzdata, double *xyzknot, int n, int k,
             int m, int l, double *work2, double *work3, double *bcoef, int nxmax, int nymax);
static void spli2d (double *xyvec, int ld, double *xydata, double *xyknot, int n, int k, int m,
             double *work2, double *work3, double *bcoef);
static int dbs2in (int nx, double *xvec, int ny,
                             double *yvec, double *xydata, int ldf, int kx,
                             int ky, double *xknot, double *yknot, double *bcoef);
static double dbs2vl (double x, double y, int kx, int ky,
               double *xknot, double *yknot, int nx, int ny, double *bcoef);
static void dbs3in (int nx, double *xvec, int ny,
             double *yvec, int nz, double *zvec, double *xyzdata,
             int ldf, int mdf, int zdf,
             int kx, int ky, int kz, double *xknot, double *yknot, double *zknot, double *bcoef);
static double dbs3vl (double x, double y, double z__, int kx,
               int ky, int kz, double *xknot, double *yknot,
               double *zknot, int nx, int ny, int nz, double *bcoef);
static void get_biats (int nxvec, double *xvec, int nyvec, double *yvec,
                int nzvec, double *zvec, int kx, int ky, int kz,
                double *xknot, double *yknot, double *zknot,
                int nx, int ny, int nz,
                double *biatx, double *biaty, double *biatz, int *leftx, int *lefty, int *leftz);
static void dbs3gd2 (int nxvec, int nyvec, int nzvec,
              int kx, int ky, int kz,
              int nx, int ny, int nz,
              double *bcoef, double *value,
              int ldvalue, int mdvalue, int zdvalue,
              double *biatx, double *biaty, double *biatz, int *leftx, int *lefty, int *leftz);


/* Table of constant values */

static int c__1 = 1;
static int c__0 = 0;

/* -*-fortran-*- */


/*   VERSION 2.2 */

/*   This library contains routines for B-spline interpolation in */
/*   one, two, and three dimensions. Part of the routines are based */
/*   on the book by Carl de Boor: A practical guide to Splines (Springer, */
/*   New-York 1978) and have the same calling sequence and names as */
/*   the corresponding routines from the IMSL library. For documen- */
/*   tation see the additional files. NOTE: The results in the demo */
/*   routines may vary slightly on different architectures. */

/*   by W. Schadow 07/19/98 */
/*   last changed by W. Schadow 07/28/2000 */


/*   Wolfgang Schadow */
/*   TRIUMF */
/*   4004 Wesbrook Mall */
/*   Vancouver, B.C. V6T 2A3 */
/*   Canada */

/*   email: schadow@triumf.ca  or  schadow@physik.uni-bonn.de */

/*   www  : http://www.triumf.ca/people/schadow */


/*  ------------------------------------------------------------------ */


/*   Copyright (C) 2000 Wolfgang Schadow */

/*   This library is free software; you can redistribute it and/or */
/*   modify it under the terms of the GNU Library General Public */
/*   License as published by the Free Software Foundation; either */
/*   version 2 of the License, or (at your option) any later version. */

/*   This library is distributed in the hope that it will be useful, */
/*   but WITHOUT ANY WARRANTY; without even the implied warranty of */
/*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU */
/*   Library General Public License for more details. */

/*   You should have received a copy of the GNU Library General Public */
/*   License along with this library; if not, write to the */
/*   Free Software Foundation, Inc., 59 Temple Place - Suite 330, */
/*   Boston, MA  02111-1307, USA. */


/*  ------------------------------------------------------------------ */


/*   The following routines are included: */

/*            dbsnak */

/*            dbsint */
/*            dbsval */
/*            dbsder */
/*            dbs1gd */

/*            dbs2in */
/*            dbs2dr */
/*            dbs2vl */
/*            dbs2gd */

/*            dbs3in */
/*            dbs3vl */
/*            dbs3dr */
/*            dbs3gd */

/*  ------------------------------------------------------------------ */

/*  NEW: corrected some error messages */
/*       some changes in the checks of dbs3dg to find a possible */
/*       error earlier. */

/*  ------------------------------------------------------------------ */

/*  NEW: documentation included, changed some comments */

/*  ------------------------------------------------------------------ */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */






static void dbsnak (int nx, double *xvec, int kxord, double *xknot)
{

    /*These two variables should be known all the time, therefore they are
     * static*/
    static int first = 1;
    static double eps;


    /* Local variables */
    int i;
    extern double dlamch (char *, int);

    /* Function Body */
    if (first)
    {

        first = 0;

        /*dlamch determines double precision machine parameters
         * It should be LAPACK routine*/
        eps = dlamch ("precision", (int) 9);

        /*Print out eps */
        if (pct.gridpe == 0)
            printf ("\n\n subroutine dbsnak: eps = %e", eps);

    }                           /*end if (first) */

    if (kxord < 0 || kxord > nx)
    {

        printf ("\n\n subroutine dbsnak: error");
        printf ("\n 0 <= kxord <= nx is required.");
        printf ("\n kxord = %d and nx = %d is given.", kxord, nx);

        exit (0);
    }                           /*end if (kxord < 0 || kxord > nx) */

    for (i = 0; i < kxord; ++i)
        xknot[i] = xvec[0];

    if (kxord % 2 == 0)
        for (i = kxord; i < nx; i++)
            xknot[i] = xvec[i - kxord / 2];


    else
    {

        for (i = kxord; i < nx; ++i)
            xknot[i] = 0.5 * (xvec[i - kxord / 2] + xvec[i - kxord / 2 - 1]);
    }

    for (i = nx; i < nx + kxord; i++)
        xknot[i] = xvec[nx - 1] * (1.0 + eps);

}                               /* End dbsnak */


void dbsint (int nx, double *xvec, double *xdata, int kx, double *xknot, double *bcoef);

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* Subroutine */ void dbsint (int nx, double *xvec, double *xdata,
                              int kx, double *xknot, double *bcoef)
{

    /* Local variables */
    int i, jj, ik, ix, kxm1, lenq;
    double *work;
    int kpkm2, iflag;
    double xveci;
    int leftx;

    /* Function Body */
    if (kx > 8)
    {

        printf ("\n\n subroutine dbsint: error");
        printf ("\n kx > 8");
        printf ("\n kx = %d", kx);

        exit (0);
    }


    kxm1 = kx - 1;
    kpkm2 = kxm1 << 1;
    leftx = kx;
    lenq = nx * (kx + kxm1);

    /*Memory for work */
    my_malloc (work, lenq, double);

    for (i = 0; i < lenq; ++i)
        work[i] = 0.;

    /*At the beginning of this loop, leftx should be index of the
     * last zero knot*/
    for (ix = 0; ix < nx; ++ix)
    {
        xveci = xvec[ix];
        leftx = rmg_max (leftx, ix + 1);

        if (xveci < xknot[leftx - 1])
        {
            /*L998: */
            printf ("\n\n subroutine dbsint:");
            printf ("\n xknot(ix) <= xknot(ix+1) required.");
            printf ("\n ix %d, xknot[ix] %f xknot[ix+1] %f", ix, xknot[ix], xknot[ix + 1]);
            exit (0);
        }

        /*I am not 100% sure this is exactly true,
         * if this does not work, uncomment the while loop below*/
        while (xknot[leftx] <= xveci)
            leftx++;


#if 0
        int i2 = ix + kx + 1;
        int nxp1 = nx + 1;
        int ilp1mx = rmg_min (i2, nxp1);
        /*This should probably find leftx so that xvec[ix] < xknot[leftx] */
        while (1)
        {
            if (xveci < xknot[leftx])
                break;
            ++leftx;

            if (leftx >= ilp1mx)
            {

                --leftx;
                if (xveci > xknot[leftx])
                {

                    /*This is L998 */
                    printf ("\n\n subroutine dbsint:");
                    printf ("\n xknot(ix) <= xknot(ix+1) required.");
                    printf ("\n ix %d, xknot[ix] %f xknot[ix+1] %f", ix, xknot[ix], xknot[ix + 1]);
                    exit (0);
                }

                break;
            }


        }                       /*end  while (1) */
#endif


        bsplvb (xknot, kx, 1, xveci, leftx, bcoef);

        jj = ix - leftx + 2 + (leftx - kx) * (kx + kxm1);

        for (ik = 0; ik < kx; ++ik)
        {
            jj += kpkm2;
            work[jj - 1] = bcoef[ik];
        }

    }                           /*end for (ix = 0; ix < nx; ++ix) */


    //i1 = kx + kxm1;
    //banfac_(work, &i1, &nx, &kxm1, &kxm1, &iflag);
    banfac (work, kx + kxm1, nx, kxm1, kxm1, &iflag);

    if (iflag == 2)
    {
        printf ("\n\n subroutine dbsint: error: no solution of linear equation system !!!");
        exit (0);
    }

    for (ix = 0; ix < nx; ++ix)
        bcoef[ix] = xdata[ix];

    //i1 = kx + kxm1;
    //banslv_(work, &i1, &nx, &kxm1, &kxm1, bcoef);
    banslv (work, kx + kxm1, nx, kxm1, kxm1, bcoef);

    my_free (work);

}                               /* dbsint_ */



double dbsval (double x, int kx, double *xknot, int nx, double *bcoef);
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
double dbsval (double x, int kx, double *xknot, int nx, double *bcoef)
{
    /* System generated locals */
    int i1, i2;
    double ret_val;


    /* Local variables */
    double dl[8];
    int ik, il;
    double dr[8];
    int ix;
    double work[8], save1, save2;
    int leftx;


    /* Function Body */
    /*8 is maximum order */
    if (kx > 8)
    {

        printf ("\n subroutine dbsval:");
        printf ("\n kx <= kxmax required.");

        exit (0);
    }

/*     check if xknot(i) <= xknot(i+1) and calculation of i so that */
/*     xknot(i) <= x < xknot(i+1) */

    i1 = nx + kx - 1;
    //ix = 0;
    ix = kx;
    while (xknot[ix] <= x)
    {
        ix++;

        if (ix == i1)
        {
            printf (" Matching interval cannot be found");
            exit (0);
        }

    }                           /*end while (xknot[ix] < x) */

    leftx = ix;

#if 0
    leftx = 0;
    for (ix = 0; ix < i1; ++ix)
    {

        if (xknot[ix] > xknot[ix + 1])
        {
            printf ("\n subroutine dbsval:");
            printf ("\n xknot(ix) <= xknot(ix+1) required.");
            printf ("\n  ix %d xknot[ix] %f xknot[ix + 1] %f", ix, xknot[ix], xknot[ix + 1]);
            exit (0);
        }
        if (xknot[ix] <= x && x < xknot[ix + 1])
        {
            leftx = ix + 1;
        }
    }
#endif

    /*This should not happen */
#if 0
    if (leftx == 0)
    {

        printf ("\n subroutine dbsval:");
        printf ("\n ix with xknot(ix) <= x < xknot(ix+1) required.");
        printf ("\n x = %f", x);

        exit (0);
    }
#endif

    i2 = leftx - kx;
    i1 = kx - 1;

    for (ik = 0; ik < i1; ++ik)
    {
        work[ik] = bcoef[i2 + ik];
        dl[ik] = x - xknot[i2 + ik];
        dr[ik] = xknot[leftx + ik] - x;
    }

    i1 = kx - 1;
    i2 = leftx - 1;
    work[i1] = bcoef[i2];
    dl[i1] = x - xknot[i2];

    for (ik = 0; ik < i1; ++ik)
    {

        save2 = work[ik];

        for (il = ik + 1; il < kx; ++il)
        {
            save1 = work[il];
            work[il] = (dl[il] * work[il] + dr[il - ik - 1] * save2) / (dl[il] + dr[il - ik - 1]);
            save2 = save1;
        }

    }

    ret_val = work[kx - 1];
    return ret_val;

}                               /* dbsval_ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* Subroutine */ static void bsplvb (double *t, int jhigh,
                                     int index, double x, int left, double *biatx)
{
    /* Initialized data */

    int j = 1;


    /* Local variables */
    int i;
    double dl[8], dr[8];
    int jp1;
    double term, saved;


    /* Function Body */

    if (index != 2)
    {
        /*L10: */
        j = 1;
        biatx[0] = 1.;
        if (j >= jhigh)
            return;
    }

/*L20:*/
    do
    {
        jp1 = j + 1;
        dr[j - 1] = t[left + j - 1] - x;
        dl[j - 1] = x - t[left - j];
        saved = 0.;

        for (i = 0; i < j; ++i)
        {
            term = biatx[i] / (dr[i] + dl[jp1 - i - 2]);
            biatx[i] = saved + dr[i] * term;
            saved = dl[jp1 - i - 2] * term;
        }

        biatx[jp1 - 1] = saved;
        j = jp1;

    }
    while (j < jhigh);
}                               /* bsplvb_ */







/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* Subroutine */ static void banfac (double *w, int nroww, int nrow,
                                     int nbandl, int nbandu, int *iflag)
{
    /* System generated locals */
    int w_dim1, w_offset, i2, i3;

    /* Local variables */
    int i, j, k, ipk, jmax, kmax, midmk;
    double pivot;
    int nrowm1, middle;
    double factor;

    /* Parameter adjustments */
    w_dim1 = nroww;
    w_offset = 1 + w_dim1;
    //w -= w_offset;

    /* Function Body */
    *iflag = 1;
    middle = nbandu + 1;
    nrowm1 = nrow - 1;

    /*This should not happen, number of grid points should be reasonable number */
#if 0
    if (nrowm1 < 0)
    {
        /*L999 */
        *iflag = 2;
        return;
    }
#endif

    if (!nrowm1)
    {

        /*L900 */
        if (w[middle + nrow * w_dim1 - w_offset] != 0.)
            return;

    }

    /*Also nbandl and nbandu should be positive, otherwise we are using 
     * order 1 or less*/
#if 0
    if (nbandl <= 0)
    {

        i1 = nrowm1;
        for (i = 0; i < i1; ++i)
        {
            if (w[middle + i * w_dim1 - 1] == 0.)
            {
                /*L999 */
                *iflag = 2;
                return;
            }
/* L20: */
        }

        /*L900 */
        if (w[middle + nrow * w_dim1 - w_offset] != 0.)
            return;

    }

/*L30:*/
    if (nbandu <= 0)
    {

        i1 = nrowm1;
        for (i = 0; i < i1; ++i)
        {
            pivot = w[middle + i * w_dim1 - 1];
            if (pivot == 0.)
            {
                /*L999 */
                *iflag = 2;
                return;
            }
            i2 = nbandl, i3 = nrow - i - 1;
            jmax = min (i2, i3);
            i2 = jmax;

            for (j = 0; j < i2; ++j)
                w[middle + j + i * w_dim1] /= pivot;
            /* L50: */
        }

        return;
    }
#endif


    /*L60: */
    for (i = 0; i < nrowm1; ++i)
    {

        pivot = w[middle + i * w_dim1 - 1];
        if (pivot == 0.)
        {
            /*L999 */
            *iflag = 2;
            return;
        }
        i2 = nbandl, i3 = nrow - i - 1;
        jmax = rmg_min (i2, i3);

        for (j = 0; j < jmax; ++j)
            w[middle + j + i * w_dim1] /= pivot;

        i2 = nbandu, i3 = nrow - i - 1;
        kmax = rmg_min (i2, i3);

        for (k = 0; k < kmax; ++k)
        {
            ipk = i + 2 + k;
            midmk = middle - k - 1;
            factor = w[midmk + ipk * w_dim1 - w_offset];
            i3 = jmax;

            for (j = 0; j < i3; ++j)
                w[midmk + j + (ipk - 1) * w_dim1] -= w[middle + j + i * w_dim1] * factor;
        }

    }                           /*end for (i = 0; i < nrowm1; ++i) */
/*L900:*/
    if (w[middle + nrow * w_dim1 - w_offset] != 0.)
        return;

    /*L999: */
    *iflag = 2;
    return;
}                               /* banfac */




/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* Subroutine */ static void banslv (double *w, int nroww, int nrow,
                                     int nbandl, int nbandu, double *b)
{
    /* System generated locals */
    int w_dim1, w_offset, i1, i2, i3;

    /* Local variables */
    int i, j, jmax, nrowm1, middle;

    w_dim1 = nroww;
    w_offset = 1 + w_dim1;
    //w -= w_offset;

    /* Function Body */
    middle = nbandu + 1;
    if (nrow == 1)
    {

        /*L99 */
        b[0] /= w[middle + w_dim1 - w_offset];
        return;
    }

    nrowm1 = nrow - 1;

    if (nbandl != 0)
    {

        i1 = nrowm1;
        for (i = 0; i < i1; ++i)
        {
            /* Computing MIN */
            i2 = nbandl, i3 = nrow - i - 1;
            jmax = rmg_min (i2, i3);
            i2 = jmax;

            for (j = 0; j < i2; ++j)
                b[i + j + 1] -= b[i] * w[middle + j + i * w_dim1];

        }

    }

/*L30:*/
    if (nbandu <= 0)
    {


        for (i = 0; i < nrow; ++i)
            b[i] /= w[i * w_dim1];

        return;
    }

    for (i = nrow - 1; i >= 1; --i)
    {
        b[i] /= w[middle + i * w_dim1 - 1];
        i1 = nbandu, i2 = i;

        jmax = rmg_min (i1, i2);
        i1 = jmax;

        for (j = 0; j < i1; ++j)
            b[i - j - 1] -= b[i] * w[middle - j + i * w_dim1 - 2];
    }


    /*L99: */
    b[0] /= w[middle + w_dim1 - w_offset];
    return;
}                               /* banslv_ */



/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
static double dbsdca (int iderx, double x, int kx, double *xknot, double *bcoef, int leftx)
{
    /* System generated locals */
    int i1, i2;
    double ret_val;

    /* Local variables */
    int i;
    double y, dl[8];
    int ik, il;
    double dr[8], dik, bsp[8], sum, save, work[8], save1, save2;


/* This routine is equivalent to the routine dbsder, but it does not */
/* check the parameters!!c */

/* Evaluates the derivative of a spline, given its B-spline representation. */


/*   iderx  - order of the derivative to be evaluated.  (input) */
/*            in particular, iderx = 0 returns the value of the */
/*            spline. */
/*   x      - point at which the spline is to be evaluated.  (input) */
/*   kx     - order of the spline.  (input) */
/*   xknot  - array of length nx+kx containing the knot */
/*            sequence.  (input) */
/*            xknot must be nondecreasing. */
/*   nx     - number of B-spline coefficients.  (input) */
/*   bcoef  - array of length nx containing the B-spline */
/*            coefficients.  (input) */
/*   leftx  - number of the intervall of xknot that includes x */
/*   dbsdca - value of the ideriv-th derivative of the spline at x. */
/*            (output) */


    /* Function Body */
    if (iderx == 0)
    {
        i1 = kx - 1;
        for (ik = 0; ik < i1; ++ik)
        {
            work[ik] = bcoef[leftx + ik - kx];
            dl[ik] = x - xknot[leftx + ik - kx];
            dr[ik] = xknot[leftx + ik] - x;
/* L20: */
        }
        work[kx - 1] = bcoef[leftx - 1];
        dl[kx - 1] = x - xknot[leftx - 1];
        i1 = kx - 1;
        for (ik = 0; ik < i1; ++ik)
        {
            save2 = work[ik];
            i2 = kx;
            for (il = ik + 1; il < i2; ++il)
            {
                save1 = work[il];
                work[il] = (dl[il] * work[il] + dr[il - ik - 1] *
                            save2) / (dl[il] + dr[il - ik - 1]);
                save2 = save1;
/* L40: */
            }
/* L30: */
        }
        ret_val = work[kx - 1];
    }
    else if (iderx >= 1 && iderx < kx)
    {
        bsp[0] = 1.;
        i1 = kx - iderx - 1;
        for (ik = 0; ik < i1; ++ik)
        {
            dr[ik] = xknot[leftx + ik] - x;
            dl[ik] = x - xknot[leftx - ik - 1];
            save = bsp[0];
            bsp[0] = 0.;
            i2 = ik + 1;
            for (il = 0; il < i2; ++il)
            {
                y = save / (dr[il - 1] + dl[ik - il]);
                bsp[il] += dr[il] * y;
                save = bsp[il + 1];
                bsp[il + 1] = dl[ik - il] * y;
/* L60: */
            }
/* L50: */
        }
        i1 = kx;
        for (ik = 0; ik < i1; ++ik)
        {
            work[ik] = bcoef[leftx + ik - kx];
            dr[ik] = xknot[leftx + ik] - x;
            dl[ik] = x - xknot[leftx + ik - kx];
/* L70: */
        }
        i1 = iderx;
        for (ik = 0; ik < i1; ++ik)
        {
            dik = (double) (kx - ik - 1);
            save2 = work[ik];
            i2 = kx;
            for (il = ik + 1; il < i2; ++il)
            {
                save1 = work[il];
                work[il] = dik * (work[il] - save2) / (dl[il] + dr[il - ik - 1]);
                save2 = save1;
/* L90: */
            }
/* L80: */
        }
        sum = 0.;
        i1 = kx - iderx;
        for (i = 1; i <= i1; ++i)
        {
            sum += bsp[i - 1] * work[iderx + i - 1];
/* L100: */
        }
        ret_val = sum;
    }
    else
    {
        ret_val = 0.;
    }
    return ret_val;
}                               /* dbsdca */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* Subroutine */ static int dbs2in (int nx, double *xvec, int ny,
                             double *yvec, double *xydata, int ldf, int kx,
                             int ky, double *xknot, double *yknot, double *bcoef)
{


    /* Local variables */
#if 0
    double work1[16900] /* was [130][130] */ , work2[130], work3[1950];
#endif
    double *work1, *work2, *work3;
    int maxnxny, max2;

/*  Computes a two-dimensional tensor-product spline interpolant, */
/*  returning the tensor-product B-spline coefficients. */

/*    nx     - number of data points in the x-direction.  (input) */
/*    xvec   - array of length nx containing the data points in */
/*             the x-direction.  (input) */
/*             xdata must be strictly increasing. */
/*    ny     - number of data points in the y-direction.  (input) */
/*    yvec   - array of length ny containing the data points in */
/*             the y-direction.  (input) */
/*             ydata must be strictly increasing. */
/*    xydata - array of size nx by nydata containing the values to */
/*             be interpolated.  (input) */
/*             fdata(i,j) is the value at (xdata(i),ydata(j)). */
/*    ldf    - the leading dimension of fdata exactly as specified in */
/*             the dimension statement of the calling program. */
/*             (input) */
/*    kx     - order of the spline in the x-direction.  (input) */
/*             kxord must be less than or equal to nxdata. */
/*    ky     - order of the spline in the y-direction.  (input) */
/*             kyord must be less than or equal to nydata. */
/*    xknot  - array of length nx+kx containing the knot */
/*             sequence in the x-direction.  (input) */
/*             xknot must be nondecreasing. */
/*    yknot  - array of length ny+ky containing the knot */
/*             sequence in the y-direction.  (input) */
/*             yknot must be nondecreasing. */
/*    bcoef  - array of length nx*ny containing the */
/*             tensor-product B-spline coefficients.  (output) */
/*             bscoef is treated internally as a matrix of size nxdata */
/*             by nydata. */


/*     dimensions should be */
/*                  work1(max(nx,ny),max(nx,ny)) */
/*                  work2(max(nx,ny)) */
/*                  work3(max((2*kx-1)*nx,(2*ky-1)*ny)) */
    maxnxny = rmg_max (nx, ny);
    max2 = rmg_max ((2 * kx - 1) * nx, (2 * ky - 1) * ny);

    my_malloc (work1, maxnxny * maxnxny, double);
    my_malloc (work2, maxnxny, double);
    my_malloc (work3, max2, double);


    /* Function Body */
#if 0
    if (*kx > 8)
    {
        printf ("\n\n subroutine dbs2in: error");
        printf ("\n kx > kxmax");
        printf ("\n kx = %d  kxmax = %d", *kx, 8);

        exit (0);
    }
    if (*ky > 8)
    {

        printf ("\n\n subroutine dbs2in: error");
        printf ("\n ky > kymax");
        printf ("\n ky = %d  kymax = %d", *ky, 8);

        exit (0);
    }
#endif
    spli2d (xvec, ldf, xydata, xknot, nx, kx, ny, work2, work3, work1);
    spli2d (yvec, ny, work1, yknot, ny, ky, nx, work2, work3, bcoef);

    my_free (work1);
    my_free (work2);
    my_free (work3);
    return 0;
}                               /* dbs2in */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
static void spli2d (double *xyvec, int ld, double *xydata, double *xyknot, int n, int k, int m,
             double *work2, double *work3, double *bcoef)
{
    /* System generated locals */
    int xydata_dim1, bcoef_dim1, i1, i2;


    /* Local variables */
    int i, j, jj, km1, left, lenq, kpkm2, iflag;
    double xyveci;

    /* Parameter adjustments */
    bcoef_dim1 = m;
    xydata_dim1 = ld;

    /* Function Body */
    km1 = k - 1;
    kpkm2 = km1 << 1;
    left = k;
    lenq = n * (k + km1);
    i1 = lenq;
    for (i = 0; i < i1; ++i)
    {
        work3[i] = 0.;
/* L10: */
    }
    i1 = n;
    for (i = 0; i < i1; ++i)
    {
        xyveci = xyvec[i];
/* Computing MIN */
        i2 = i + 1 + k;
        left = rmg_max (left, i + 1);

#if 0
        int np1 = n + 1;
        int ilp1mx = rmg_min (i2, np1);
        if (xyveci < xyknot[left - 1])
        {
            /*L998: */
            printf ("\n\n subroutine db2in:");
            printf ("\n i with knot(i) <= x/y < knot(i+1) required.");
            printf ("\n knot(1)     =  %f", xyknot[0]);
            printf ("\n knot(n+k)   =  %f", xyknot[n + k - 1]);
            printf ("\n        x/y  =  %f", xyveci);

            exit (0);
        }
#endif

        /*I am not 100% sure about this, in case of problems,
         * comment this while loop and uncomment the commented code below*/
        while (xyveci >= xyknot[left])
            ++left;
#if 0
      L30:
        if (xyveci < xyknot[left])
        {
            goto L40;
        }
        ++left;
        if (left < ilp1mx)
        {
            goto L30;
        }
        --left;
        if (xyveci > xyknot[left])
        {
            printf ("\n\n subroutine db2in:");
            printf ("\n i with knot(i) <= x/y < knot(i+1) required.");
            printf ("\n knot(1)     =  %f", xyknot[0]);
            printf ("\n knot(n+k)   =  %f", xyknot[n + *k - 1]);
            printf ("\n        x/y  =  %f", xyveci);

            exit (0);
        }
      L40:
#endif
        i2 = n + k;
        bsplvb (xyknot, k, c__1, xyveci, left, work2);
        jj = i - left + 1 + (left - k) * (k + km1);
        i2 = k;
        for (j = 0; j < i2; ++j)
        {
            jj += kpkm2;
            work3[jj] = work2[j];
/* L50: */
        }
/* L20: */
    }
    i1 = k + km1;
    banfac (work3, i1, n, km1, km1, &iflag);

    if (iflag == 2)
    {
        /*L999: */
        printf ("\n\n subroutine db2in:");
        printf ("\n no solution of linear equation system !!!");
        exit (0);
    }

/*L60:*/
    i1 = m;
    for (j = 0; j < i1; ++j)
    {
        i2 = n;
        for (i = 0; i < i2; ++i)
        {
            work2[i] = xydata[i + j * xydata_dim1];
/* L80: */
        }
        i2 = k + km1;
        banslv (work3, i2, n, km1, km1, work2);
        i2 = n;
        for (i = 0; i < i2; ++i)
        {
            bcoef[j + i * bcoef_dim1] = work2[i];
/* L90: */
        }
/* L70: */
    }
}                               /* spli2d */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
static double dbs2vl (double x, double y, int kx, int ky,
               double *xknot, double *yknot, int nx, int ny, double *bcoef)
{
    /* System generated locals */
    int bcoef_dim1, i1;
    double ret_val;


    /* Local variables */
    int i, iky;
    double work[8];
    int leftx, lefty;

/*  evaluates a two-dimensional tensor-product spline, given its */
/*  tensor-product B-spline representation.    use numeric */

/*   x      - x-coordinate of the point at which the spline is to be */
/*            evaluated.  (input) */
/*   y      - y-coordinate of the point at which the spline is to be */
/*            evaluated.  (input) */
/*   kx     - order of the spline in the x-direction.  (input) */
/*   ky     - order of the spline in the y-direction.  (input) */
/*   xknot  - array of length nx+kx containing the knot */
/*            sequence in the x-direction.  (input) */
/*            xknot must be nondecreasing. */
/*   yknot  - array of length ny+ky containing the knot */
/*            sequence in the y-direction.  (input) */
/*            yknot must be nondecreasing. */
/*   nx     - number of B-spline coefficients in the x-direction. */
/*            (input) */
/*   ny     - number of B-spline coefficients in the y-direction. */
/*            (input) */
/*   bcoef  - array of length nx*ny containing the */
/*            tensor-product B-spline coefficients.  (input) */
/*            bscoef is treated internally as a matrix of size nx */
/*            by ny. */
/*   dbs2vl - value of the spline at (x,y).  (output) */


/*     check if k <= kmax */

    /* Parameter adjustments */
    bcoef_dim1 = nx;
    /*bcoef_offset = 1 + bcoef_dim1 * 1;
       bcoef -= bcoef_offset; */

    /* Function Body */
    if (kx > 8)
    {
        printf ("\n\n subroutine dbs2vl:");
        printf ("\n kx <= kxmax required.");

        exit (0);
    }
    if (ky > 8)
    {
        printf ("\n\n subroutine dbs2vl:");
        printf ("\n kx <= kxmax required.");

        exit (0);
    }

/*     check if knot(i) <= knot(i+1) and calculation of i so that */
/*     knot(i) <= x < knot(i+1) */

    leftx = 0;
    i1 = nx + kx - 1;
    for (i = 0; i < i1; ++i)
    {
        if (xknot[i] > xknot[i + 1])
        {

            printf ("\n\n subroutine dbs2vl:");
            printf ("\n xknot(i) <= xknot(i+1) required.");
            printf ("\n i %d xknot(i) %f xknot(i+1) %f", i + 1, xknot[i], xknot[i + 1]);

            exit (0);
        }
        if (xknot[i] <= x && x < xknot[i + 1])
        {
            leftx = i + 1;
        }
/* L10: */
    }
    if (leftx == 0)
    {
        printf ("\n\n subroutine dbs2vl:");
        printf ("\n i with xknot(i) <= x < xknot(i+1) required.");
        printf ("\n xknot(i)   = %f", xknot[i]);
        printf ("\n x = %f", x);
        printf ("\n xknot(i+1) = %f", xknot[i + 1]);

        exit (0);
    }
    lefty = 0;
    i1 = ny + ky - 1;
    for (i = 0; i < i1; ++i)
    {
        if (yknot[i] > yknot[i + 1])
        {
            printf ("\n\n subroutine dbs2vl:");
            printf ("\n yknot(i) <= yknot(i+1) required.");
            printf ("\n i %d yknot(i) %f yknot(i+1) %f", i + 1, yknot[i], yknot[i + 1]);

            exit (0);
        }
        if (yknot[i] <= y && y < yknot[i + 1])
        {
            lefty = i + 1;
        }
/* L20: */
    }
    if (lefty == 0)
    {
        printf ("\n\n subroutine dbs2vl:");
        printf ("\n i with yknot(i) <= y < yknot(i+1) required.");
        printf ("\n yknot(i)   = %f", yknot[i]);
        printf ("\n y = %f", y);
        printf ("\n yknot(i+1) = %f", yknot[i + 1]);
    }
    i1 = ky;
    for (iky = 0; iky < i1; ++iky)
    {
        work[iky] = dbsdca (c__0, x, kx, xknot, &bcoef[(lefty - ky + iky) * bcoef_dim1], leftx);
/* L30: */
    }
    ret_val = dbsval (y, ky, &yknot[lefty - ky], ky, work);
    return ret_val;
}                               /* dbs2vl */



/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
static void dbs3in (int nx, double *xvec, int ny,
             double *yvec, int nz, double *zvec, double *xyzdata,
             int ldf, int mdf, int zdf,
             int kx, int ky, int kz, double *xknot, double *yknot, double *zknot, double *bcoef)
{
    /* System generated locals */
    int bcoef_dim1, bcoef_dim2, i1;


    /* Local variables */
    int iz;
    //double work1[2197000];   /* was [130][130][130] */ 
    /*double work2[130]; */
    /*double        work3[1950]; */
    double *work1;              /* was [130][130][130] */
    double *work2;
    double *work3;

/*  Computes a three-dimensional tensor-product spline interpolant, */
/*  returning the tensor-product B-spline coefficients. */

/*   nx      - number of data points in the x-direction.  (input) */
/*   xvec    - array of length nxdata containing the data points in */
/*             the x-direction.  (input) */
/*             xdata must be increasing. */
/*   ny      - number of data points in the y-direction.  (input) */
/*   yvec    - array of length nydata containing the data points in */
/*             the y-direction.  (input) */
/*             ydata must be increasing. */
/*   nz      - number of data points in the z-direction.  (input) */
/*   zvec    - array of length nzdata containing the data points in */
/*             the z-direction.  (input) */
/*             zdata must be increasing. */
/*   xyzdata - array of size nx by ny by nz containing the */
/*             values to be interpolated.  (input) */
/*             xyzdata(i,j,k) contains the value at */
/*             (xvec(i),yvec(j),zvec(k)). */
/*   ldf     - leading dimension of fdata exactly as specified in the */
/*             dimension statement of the calling program.  (input) */
/*   mdf     - middle dimension of fdata exactly as specified in the */
/*             dimension statement of the calling program.  (input) */
/*   kx      - order of the spline in the x-direction.  (input) */
/*             kxord must be less than or equal to nxdata. */
/*   ky      - order of the spline in the y-direction.  (input) */
/*             kyord must be less than or equal to nydata. */
/*   kz      - order of the spline in the z-direction.  (input) */
/*             kzord must be less than or equal to nzdata. */
/*   xknot   - array of length nx+kx containing the knot */
/*             sequence in the x-direction.  (input) */
/*             xknot must be nondecreasing. */
/*   yknot   - array of length ny+ky containing the knot */
/*             sequence in the y-direction.  (input) */
/*             yknot must be nondecreasing. */
/*   zknot   - array of length nz+kz containing the knot */
/*             sequence in the z-direction.  (input) */
/*             zknot must be nondecreasing. */
/*   bcoef   - array of length nx*ny*nz containing the */
/*             tensor-product B-spline coefficients.  (output) */
/*             bscoef is treated internally as a matrix of size nx */
/*             by ny by nz. */


/*     dimensions should be */
/*              work1(nx,ny,nz) */
/*              work2(nz) */
/*              work3((2*kz-1)*nz) */


    my_malloc (work1, nx * ny * nz, double);
    my_malloc (work2, nz, double);
    my_malloc (work3, (2 * kz - 1) * nz, double);
    /* Parameter adjustments */
    //--xvec;
    //--yvec;
    bcoef_dim1 = nx;
    bcoef_dim2 = ny;
    //bcoef_offset = 1 + bcoef_dim1 * (1 + bcoef_dim2 * 1);
    //bcoef -= bcoef_offset;
    //--zvec;
#if 0
    xyzdata_dim1 = ldf;
    xyzdata_dim2 = mdf;
    //xyzdata_offset = 1 + xyzdata_dim1 * (1 + xyzdata_dim2 * 1);
    xyzdata_offset = 1 + 1 * xyzdata_dim1 + 1 * xyzdata_dim1 * xyzdata_dim2;
    xyzdata -= xyzdata_offset;
#endif
    //--xknot;
    //--yknot;
    //--zknot;



    /* Function Body */
    if (kx > 8)
    {
        printf ("\n\n subroutine dbs3in: error");
        printf ("\n kx > kxmax");
        printf ("\n kx = %d kxmax = %d", kx, 8);

        exit (0);
    }
    if (ky > 8)
    {
        printf ("\n\n subroutine dbs3in: error");
        printf ("\n ky > kymax");
        printf ("\n ky = %d kymax = %d", ky, 8);

        exit (0);
    }
    if (kz > 8)
    {
        printf ("\n\n subroutine dbs3in: error");
        printf ("\n kz > kzmax");
        printf ("\n kz = %d kzmax = %d", kz, 8);

        exit (0);

    }
    spli3d (zvec, ldf, mdf, zdf, xyzdata, zknot, nz, kz, nx, ny, work2, work3, work1, nx, ny);

    i1 = nz;
    for (iz = 0; iz < i1; ++iz)
    {
        dbs2in (nx, xvec, ny, yvec,
                //&work1[(iz * 130 + 1) * 130 - 17030], &c__130, 
                /*for some reason they are passing (iz-1),0,0 and not iz,1,1 */
                &work1[iz * ny * nx], nx,
                //kx, ky, &xknot[1], &yknot[1], &bcoef[(iz * bcoef_dim2 + 1) * bcoef_dim1 + 1]);
                kx, ky, xknot, yknot, &bcoef[iz * bcoef_dim2 * bcoef_dim1]);
        //exit(0);
/* L10: */
    }
    my_free (work1);
    my_free (work2);
    my_free (work3);
}                               /* dbs3in */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
static void spli3d (double *xyzvec, int ldf, int mdf, int zdf,
             double *xyzdata, double *xyzknot, int n, int k,
             int m, int l, double *work2, double *work3, double *bcoef, int nxmax, int nymax)
{
    /* System generated locals */
    int xyzdata_dim2, xyzdata_dim3, bcoef_dim1, bcoef_dim2, i1, i2, i3;

    /* Local variables */
    int i, j, jj, in, km1, left, lenq, kpkm2, iflag;
    double xyzveci;

    /* Parameter adjustments */
    //xyzdata_dim1 = ldf;
    xyzdata_dim2 = mdf;
    xyzdata_dim3 = zdf;
    //xyzdata_offset = 1 + xyzdata_dim1 * (1 + xyzdata_dim2 * 1);
    //xyzdata_offset = 1 + 1*xyzdata_dim1  + 1*xyzdata_dim2*xyzdata_dim1;
    //xyzdata -= xyzdata_offset;
    //--work2;
    //--xyzvec;
    //--work3;
    //--xyzknot;
    bcoef_dim1 = nxmax;
    bcoef_dim2 = nymax;
    //bcoef_offset = 1 + bcoef_dim1 * (1 + bcoef_dim2 * 1);
    //bcoef -= bcoef_offset;


    /* Function Body */
    km1 = k - 1;
    kpkm2 = km1 << 1;
    left = k;
    lenq = n * (k + km1);
    i1 = lenq;
    for (i = 0; i < i1; ++i)
    {
        work3[i] = 0.;
/* L10: */
    }
    i1 = n;
    for (i = 0; i < i1; ++i)
    {
        xyzveci = xyzvec[i];
/* Computing MIN */
        i2 = i + 1 + k;
        left = rmg_max (left, i + 1);
#if 0
        int np1 = n + 1;
        int ilp1mx = rmg_min (i2, np1);
        if (xyzveci < xyzknot[left - 1])
        {

            /*L998: */
            printf ("\n\n subroutine db3in:");
            printf ("\n i with knot(i) <= x/y/z < knot(i+1) required.");
            printf ("\n knot(1)   = %f", xyzknot[1]);
            printf ("\n knot(n+k) =  %f", xyzknot[n + k - 1]);
            printf ("\n x/y/z = %f", xyzveci);
            exit (0);
        }
#endif

        /*I am not 100% sure about this, in case of problems,
         * comment this while loop and uncomment the commented code below*/
        while (xyzveci >= xyzknot[left])
            ++left;

#if 0
      L30:
        if (xyzveci < xyzknot[left])
        {
            goto L40;
        }
        ++left;
        if (left < ilp1mx)
        {
            goto L30;
        }
        --left;
        if (xyzveci > xyzknot[left])
        {

            /*L998: */
            printf ("\n\n subroutine db3in:");
            printf ("\n i with knot(i) <= x/y/z < knot(i+1) required.");
            printf ("\n knot(1)   = %f", xyzknot[1]);
            printf ("\n knot(n+k) =  %f", xyzknot[n + k - 1]);
            printf ("\n x/y/z = %f", xyzveci);
            exit (0);
        }
      L40:
#endif
        i2 = n + k;
        bsplvb (xyzknot, k, c__1, xyzveci, left, work2);
        jj = i - left + 2 + (left - k) * (k + km1);
        i2 = k;
        for (j = 0; j < i2; ++j)
        {
            jj += kpkm2;
            work3[jj - 1] = work2[j];
/* L50: */
        }
/* L20: */
    }
    i1 = k + km1;
    banfac (work3, i1, n, km1, km1, &iflag);

    if (iflag == 2)
    {
        /*L999: */
        printf ("\n\n subroutine db3in: error");
        printf ("\n no solution of linear equation system !!!");
        exit (0);
    }

    i1 = l;
    for (j = 0; j < i1; ++j)
    {
        i2 = m;
        for (i = 0; i < i2; ++i)
        {
            i3 = n;
            for (in = 0; in < i3; ++in)
            {
                //work2[in] = xyzdata[i + (j + in * xyzdata_dim2) * xyzdata_dim1];
                //work2[in] = xyzdata[i + j*xyzdata_dim1 + in * xyzdata_dim2*xyzdata_dim1];
                work2[in] = xyzdata[i * xyzdata_dim2 * xyzdata_dim3 + j * xyzdata_dim3 + in];
/* L90: */
            }
            i3 = k + km1;
            banslv (work3, i3, n, km1, km1, work2);
            i3 = n;
            for (in = 0; in < i3; ++in)
            {
                //bcoef[i + (j + in * bcoef_dim2) * bcoef_dim1] = work2[in];
                bcoef[i + j * bcoef_dim1 + in * bcoef_dim2 * bcoef_dim1] = work2[in];
/* L100: */
            }
/* L80: */
        }
/* L70: */
    }

}                               /* spli3d */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
static double dbs3vl (double x, double y, double z__, int kx,
               int ky, int kz, double *xknot, double *yknot,
               double *zknot, int nx, int ny, int nz, double *bcoef)
{
    /* System generated locals */
    int bcoef_dim1, bcoef_dim2, i1;
    double ret_val;


    /* Local variables */
    int i, jz;
    double work[8];
    int nintz;

/*  Evaluates a three-dimensional tensor-product spline, given its */
/*  tensor-product B-spline representation. */

/*   x      - x-coordinate of the point at which the spline is to be */
/*            evaluated.  (input) */
/*   y      - y-coordinate of the point at which the spline is to be */
/*            evaluated.  (input) */
/*   z      - z-coordinate of the point at which the spline is to be */
/*            evaluated.  (input) */
/*   kx     - order of the spline in the x-direction.  (input) */
/*   ky     - order of the spline in the y-direction.  (input) */
/*   kz     - order of the spline in the z-direction.  (input) */
/*   xknot  - array of length nx+kx containing the knot */
/*            sequence in the x-direction.  (input) */
/*            xknot must be nondecreasing. */
/*   yknot  - array of length ny+ky containing the knot */
/*            sequence in the y-direction.  (input) */
/*            yknot must be nondecreasing. */
/*   zknot  - array of length nz+kz containing the knot */
/*            sequence in the z-direction.  (input) */
/*            zknot must be nondecreasing. */
/*   nx     - number of B-spline coefficients in the x-direction. */
/*            (input) */
/*   ny     - number of B-spline coefficients in the y-direction. */
/*            (input) */
/*   nz     - number of B-spline coefficients in the z-direction. */
/*            (input) */
/*   bcoef  - array of length nx*ny*nz containing the */
/*            tensor-product B-spline coefficients.  (input) */
/*            bscoef is treated internally as a matrix of size nx */
/*            by ny by nz. */
/*   dbs3vl - value of the spline at (x,y,z).  (output) */


/*     dimension should be */
/*            dimension work(kz) */

/*     check if k <= kmax */
    /* Parameter adjustments */
    bcoef_dim1 = nx;
    bcoef_dim2 = ny;
    //bcoef_offset = 1 + bcoef_dim1 * (1 + bcoef_dim2 * 1);
    //bcoef -= bcoef_offset;

    /* Function Body */
    if (kx > 8)
    {

        printf ("\n\n subroutine dbs3vl:");
        printf ("\n kx <= kxmax required.");

        exit (0);
    }
    if (ky > 8)
    {
        printf ("\n\n subroutine dbs3vl:");
        printf ("\n ky <= kymax required.");

        exit (0);
    }
    if (kz > 8)
    {
        printf ("\n\n subroutine dbs3vl:");
        printf ("\n kz <= kzmax required.");

        exit (0);
    }

/*     check if knot(i) <= knot(i+1) and calculation of i so that */
/*     knot(i) <= x < knot(i+1) */

    nintz = 0;
    i1 = nz + kz - 1;
    for (i = 0; i < i1; ++i)
    {
        if (zknot[i] > zknot[i + 1])
        {
            printf ("\n\n subroutine dbs3vl:");
            printf ("\n zknot(i) <= zknot(i+1) required.");
            printf ("\n i %d zknot(i) %f zknot(i+1) %f", i, zknot[i], zknot[i + 1]);

            exit (0);
        }
        if (zknot[i] <= z__ && z__ < zknot[i + 1])
        {
            nintz = i + 1;
        }
/* L10: */
    }
    if (nintz == 0)
    {
        printf ("\n\n subroutine dbs3vl:");
        printf ("\n i with zknot(i) <= z < zknot(i+1) required.");
        printf ("\n zknot(i)   = %f", zknot[i]);
        printf ("\n z = %f", z__);
        printf ("\n zknot(i+1) = %f", zknot[i + 1]);

        exit (0);
    }
    i1 = kz;
    for (jz = 0; jz < i1; ++jz)
    {
        work[jz] =
            dbs2vl (x, y, kx, ky, xknot, yknot, nx, ny,
                    &bcoef[((nintz - kz + jz) * bcoef_dim2) * bcoef_dim1]);
/* L40: */
    }
    ret_val = dbsval (z__, kz, &zknot[nintz - kz], kz, work);
    return ret_val;
}                               /* dbs3vl */




/*This sets up biatx, biaty, biatz and leftx, lefty, leftz arrays*/
/*Input :nxvec - number of points in x direction for grid we want to inteprolate on
* Input :xvec - array containing x coordiantes for grid  on which we want to 
                interpolate, length of this array is given by nxvec
*Input :nyvec - number of points in y direction for grid we want to inteprolate on
* Input xvec - array containing y coordiantes for grid  on which we want to 
                interpolate, length of this array is given by nxvec
*Input :nzvec - number of points in z direction for grid we want to inteprolate on
* Input zvec - array containing z coordiantes for grid  on which we want to 
                interpolate, length of this array is given by nxvec
*Input: kx, ky, ky - order of splines in x, y and z directions
*Input:  xknot, yknot, zknot - knots in x, y and z directions, length should be
*                       nx+kx, ny+ky, nz+kz, respectivly
*Input: nx, ny, nz - number of bspline coefficients in x, y, and z directions
*Output: biatx, ,biaty, biatz - arrays (double) with lengths: nxvec*kx, nyvec*ky, nzvec*kz,
*Output: leftx, lefty, leftz, - int arrays with lengths nxvec, nyvec, nzvec*/


static void get_biats (int nxvec, double *xvec, int nyvec, double *yvec,
                int nzvec, double *zvec, int kx, int ky, int kz,
                double *xknot, double *yknot, double *zknot,
                int nx, int ny, int nz,
                double *biatx, double *biaty, double *biatz, int *leftx, int *lefty, int *leftz)
{
    /* System generated locals */
    int i1, i2, i3;

    /* Local variables */
    int i;
#if 0
    double dl[1600] /* was [200][8] */ , dr[1600] /* was [200][8] */ ;
    double biatx[1600] /* was [200][8] */ , biaty[1600] /* was [200][8] */ ,
        biatz[1600] /* was [200][8] */ ;
    double term[200];
    double save1[200];
    int leftx[200], lefty[200], leftz[200];
#endif
    double *dl, *dr, *term, *save1;
    int ik, il;
    int ix, iy, iz;
    int same;
    int next;
    int maxnvec, maxk;
    int dl_dim, dr_dim, biatx_dim, biaty_dim, biatz_dim;

    /*find max of nxvec,nyvec,nzvec */
    maxnvec = rmg_max (nxvec, nyvec);
    maxnvec = rmg_max (maxnvec, nzvec);

    /*find max of kx,ky,kz */
    maxk = rmg_max (kx, ky);
    maxk = rmg_max (maxk, kz);

    dl_dim = maxnvec;
    dr_dim = maxnvec;
    biatx_dim = nxvec;
    biaty_dim = nyvec;
    biatz_dim = nzvec;


    my_malloc (dl, maxnvec * maxk, double);
    my_malloc (dr, maxnvec * maxk, double);
    my_malloc (term, maxnvec, double);
    my_malloc (save1, maxnvec, double);


    /* Function Body */
    /*     check if kx <= kxmax */
    if (kx > 8)
    {
        printf ("\n\n subroutine dbs3gd:");
        printf ("\n kx <= kxmax required.");

        exit (0);
    }
    i1 = nx + kx - 1;
    for (i = 0; i < i1; ++i)
    {
        if (xknot[i] > xknot[i + 1])
        {
            printf ("\n\n subroutine dbs3gd:");
            printf ("\n xknot(i) <= xknot(i+1) required.");
            printf ("\n i %d xknot(i) %f xknot(i+1) %f", i + 1, xknot[i], xknot[i + 1]);

            exit (0);
        }
/* L10: */
    }
    i1 = nxvec;
    for (i = 0; i < i1; ++i)
    {
        if (xvec[i] < xknot[0] || xvec[i] > xknot[nx + kx - 1])
        {
            printf ("\n\n subroutine dbs3gd:");
            printf ("\n ix with xknot(ix) <= x < xknot(ix+1) required.");
            printf ("\n x = %f", xvec[i]);

            exit (0);
        }
/* L20: */
    }
    leftx[0] = 0;
    i1 = nx + kx;
    huntn (xknot, i1, kx, xvec[0], leftx);
    i1 = nxvec;
    for (ix = 1; ix < i1; ++ix)
    {
        leftx[ix] = leftx[ix - 1];
        same = xknot[leftx[ix] - 1] <= xvec[ix] && xvec[ix] <= xknot[leftx[ix]];
        if (!same)
        {
            ++leftx[ix];
            next = xknot[leftx[ix] - 1] <= xvec[ix] && xvec[ix] <= xknot[leftx[ix]];
            if (!next)
            {
                i2 = nx + kx;
                huntn (xknot, i2, kx, xvec[ix], &leftx[ix]);
            }
        }
/* L30: */
    }

/*     check if ky <= kymax */

    if (ky > 8)
    {
        printf ("\n\n subroutine dbs3gd:");
        printf ("\n ky <= kymax required.");

        exit (0);
    }
    i1 = ny + ky - 1;
    for (i = 0; i < i1; ++i)
    {
        if (yknot[i] > yknot[i + 1])
        {
            printf ("\n\n subroutine dbs3gd:");
            printf ("\n yknot(i) <= yknot(i+1) required.");
            printf ("\n i %d yknot(i) %f yknot(i+1) %f", i + 1, yknot[i], yknot[i + 1]);

            exit (0);
        }
/* L40: */
    }
    i1 = nyvec;
    for (i = 0; i < i1; ++i)
    {
        if (yvec[i] < yknot[0] || yvec[i] > yknot[ny + ky - 1])
        {
            printf ("\n\n subroutine dbs3gd:");
            printf ("\n iy with yknot(iy) <= y < yknot(iy+1) required.");
            printf ("\n y = %f", yvec[i + 1]);

            exit (0);
        }
/* L50: */
    }
    lefty[0] = 0;
    i1 = ny + ky;
    huntn (yknot, i1, ky, yvec[0], lefty);
    i1 = nyvec;
    for (iy = 1; iy < i1; ++iy)
    {
        lefty[iy] = lefty[iy - 1];
        same = yknot[lefty[iy] - 1] <= yvec[iy] && yvec[iy] <= yknot[lefty[iy]];
        if (!same)
        {
            ++lefty[iy];
            next = yknot[lefty[iy] - 1] <= yvec[iy] && yvec[iy] <= yknot[lefty[iy]];
            if (!next)
            {
                i2 = ny + ky;
                huntn (yknot, i2, ky, yvec[iy], &lefty[iy]);
            }
        }
/* L60: */
    }

/*     check if kz <= kzmax */

    if (kz > 8)
    {
        printf ("\n\n subroutine dbs3gd:");
        printf ("\n kz <= kzmax required.");

        exit (0);
    }
    i1 = nz + kz - 1;
    for (i = 0; i < i1; ++i)
    {
        if (zknot[i] > zknot[i + 1])
        {
            printf ("\n\n subroutine dbs3gd:");
            printf ("\n zknot(i) <= zknot(i+1) required.");
            printf ("\n i %d zknot(i) %f zknot(i+1) %f", i + 1, zknot[i], zknot[i + 1]);

            exit (0);
        }
/* L70: */
    }
    i1 = nzvec;
    for (i = 0; i < i1; ++i)
    {
        if (zvec[i] < zknot[0] || zvec[i] > zknot[nz + kz - 1])
        {
            printf ("\n\n subroutine dbs3gd:");
            printf ("\n iz with zknot(iz) <= z < zknot(iz+1) required.");
            printf ("\n z = %f", zvec[i + 1]);

            exit (0);
        }
/* L80: */
    }
    leftz[0] = 0;
    i1 = nz + kz;
    huntn (zknot, i1, kz, zvec[0], leftz);
    i1 = nzvec;
    for (iz = 1; iz < i1; ++iz)
    {
        leftz[iz] = leftz[iz - 1];
        same = zknot[leftz[iz] - 1] <= zvec[iz] && zvec[iz] <= zknot[leftz[iz]];
        if (!same)
        {
            ++leftz[iz];
            next = zknot[leftz[iz] - 1] <= zvec[iz] && zvec[iz] <= zknot[leftz[iz]];
            if (!next)
            {
                i2 = nz + kz;
                huntn (zknot, i2, kz, zvec[iz], &leftz[iz]);
            }
        }


/* L90: */
    }

    /*At this point we should have leftx, lefty and leftz setup */
    for (ix = 0; ix < nxvec; ++ix)
        biatx[ix] = 1.;

    for (ik = 0; ik < kx - 1; ++ik)
    {
        for (ix = 0; ix < nxvec; ++ix)
        {
            dr[ix + ik * dr_dim] = xknot[leftx[ix] + ik] - xvec[ix];
            dl[ix + ik * dl_dim] = xvec[ix] - xknot[leftx[ix] - ik - 1];
            save1[ix] = 0.;
/* L120: */
        }
        for (il = 0; il < ik + 1; ++il)
        {
            i3 = nxvec;
            for (ix = 0; ix < i3; ++ix)
            {
                term[ix] = biatx[ix + il * biatx_dim] /
                    (dr[ix + il * dr_dim] + dl[ix + (ik - il) * dl_dim]);
                biatx[ix + il * biatx_dim] = save1[ix] + dr[ix + il * dr_dim] * term[ix];
                save1[ix] = dl[ix + (ik - il) * dl_dim] * term[ix];
/* L140: */
            }
/* L130: */
        }
        for (ix = 0; ix < nxvec; ++ix)
            biatx[ix + (ik + 1) * biatx_dim] = save1[ix];

        /* L110: */
    }

    for (iy = 0; iy < nyvec; ++iy)
        biaty[iy] = 1.;

    for (ik = 0; ik < ky - 1; ++ik)
    {
        for (iy = 0; iy < nyvec; ++iy)
        {
            dr[iy + ik * dr_dim] = yknot[lefty[iy] + ik] - yvec[iy];
            dl[iy + ik * dl_dim] = yvec[iy] - yknot[lefty[iy] - ik - 1];
            save1[iy] = 0.;
/* L180: */
        }
        for (il = 0; il < ik + 1; ++il)
        {
            i3 = nyvec;
            for (iy = 0; iy < i3; ++iy)
            {
                term[iy] = biaty[iy + il * biaty_dim] / (dr[iy + il * dr_dim]
                                                         + dl[iy + (ik - il) * dl_dim]);
                biaty[iy + il * biaty_dim] = save1[iy] + dr[iy + il * dr_dim] * term[iy];
                save1[iy] = dl[iy + (ik - il) * dl_dim] * term[iy];
/* L200: */
            }
/* L190: */
        }
        i2 = nyvec;
        for (iy = 0; iy < i2; ++iy)
        {
            biaty[iy + (ik + 1) * biaty_dim] = save1[iy];
/* L210: */
        }
/* L170: */
    }
    for (iz = 0; iz < nzvec; ++iz)
        biatz[iz] = 1.;

    for (ik = 0; ik < kz - 1; ++ik)
    {
        for (iz = 0; iz < nzvec; ++iz)
        {
            dr[iz + ik * dr_dim] = zknot[leftz[iz] + ik] - zvec[iz];
            dl[iz + ik * dl_dim] = zvec[iz] - zknot[leftz[iz] - ik - 1];
            save1[iz] = 0.;
/* L240: */
        }
        for (il = 0; il < ik + 1; ++il)
        {

            for (iz = 0; iz < nzvec; ++iz)
            {
                term[iz] = biatz[iz + il * biatz_dim] /
                    (dr[iz + il * dr_dim] + dl[iz + (ik - il) * dl_dim]);
                biatz[iz + il * biatz_dim] = save1[iz] + dr[iz + il * dr_dim] * term[iz];
                save1[iz] = dl[iz + (ik - il) * dl_dim] * term[iz];
/* L260: */
            }
/* L250: */
        }
        i2 = nzvec;
        for (iz = 0; iz < nzvec; ++iz)
            biatz[iz + (ik + 1) * biatz_dim] = save1[iz];

        /* L230: */
    }


    /*Biats should be setup here */




    /*Release memory */
    my_free (dl);
    my_free (dr);
    my_free (term);
    my_free (save1);


}                               /* dbs3gd_ */

static void dbs3gd2 (int nxvec, int nyvec, int nzvec,
              int kx, int ky, int kz,
              int nx, int ny, int nz,
              double *bcoef, double *value,
              int ldvalue, int mdvalue, int zdvalue,
              double *biatx, double *biaty, double *biatz, int *leftx, int *lefty, int *leftz)
{
    /* System generated locals */
    int bcoef_dim1, bcoef_dim2, value_dim2, value_dim3;

    /* Local variables */
    int ix, iy, iz, ikx, iky, ikz;
    int biatx_dim, biaty_dim, biatz_dim;
    int v23, b12, lx, ly, lz, bcz, bcyz;
    double bz, byz, *ptr;

    bcoef_dim1 = nx;
    bcoef_dim2 = ny;
    value_dim2 = mdvalue;
    value_dim3 = zdvalue;

    biatx_dim = nxvec;
    biaty_dim = nyvec;
    biatz_dim = nzvec;


    v23 = value_dim3 * value_dim2;
    b12 = bcoef_dim2 * bcoef_dim1;



    for (iz = 0; iz < nzvec; ++iz)
    {
        lz = leftz[iz];

        for (iy = 0; iy < nyvec; ++iy)
        {
            ly = lefty[iy];

            for (ix = 0; ix < nxvec; ++ix)
            {
                lx = leftx[ix];
                //index = ix*v23 + iy*value_dim3 + iz;
                ptr = &value[ix * v23 + iy * value_dim3 + iz];
                *ptr = ZERO;

                for (ikz = 0; ikz < kz; ++ikz)
                {
                    bz = biatz[iz + ikz * biatz_dim];
                    bcz = (lz - kz + ikz) * b12;

                    for (iky = 0; iky < ky; ++iky)
                    {
                        byz = biaty[iy + iky * biaty_dim] * bz;
                        bcyz = (ly - ky + iky) * bcoef_dim1 + bcz;

                        for (ikx = 0; ikx < kx; ++ikx)
                        {
                            //value[index] += 
                            *ptr += biatx[ix + ikx * biatx_dim] * byz * bcoef[lx - kx + ikx + bcyz];

                        }
                    }
                }
            }
        }
    }

#if 0
    /*This was here originally instaed of above
     * While this is shorter, it takes much more computer time*/
    for (ikz = 0; ikz < kz; ++ikz)
        for (iky = 0; iky < ky; ++iky)
            for (ikx = 0; ikx < kx; ++ikx)
                for (iz = 0; iz < nzvec; ++iz)
                    for (iy = 0; iy < nyvec; ++iy)
                        for (ix = 0; ix < nxvec; ++ix)
                            value[ix * value_dim3 * value_dim2 + iy * value_dim3 + iz] +=
                                biatx[ix + ikx * biatx_dim] * biaty[iy + iky * biaty_dim] *
                                biatz[iz + ikz * biatz_dim] *
                                bcoef[leftx[ix] - kx + ikx
                                      + (lefty[iy] - ky + iky) * bcoef_dim1 + (leftz[iz] - kz +
                                                                               ikz) * bcoef_dim2 *
                                      bcoef_dim1];
#endif

}                               /* dbs3gd_ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
static void huntn (double *xx, int n, int kord, double x, int *jlo)
{
    int jm, inc, jhi, max__, null;


    /*     works only for B-Splines (order n) */

    /* Parameter adjustments */
    --xx;

    /* Function Body */
    max__ = n - kord;
    null = kord;
    if (*jlo <= null || *jlo > max__)
    {
        *jlo = null;
        jhi = max__ + 1;
        //goto L30;
    }

    else
    {
        inc = 1;
        if (x >= xx[*jlo])
        {

            do
            {
                jhi = *jlo + inc;

                if (jhi > max__)
                {
                    jhi = max__ + 1;
                    break;
                }
                if (x < xx[jhi])
                    break;
                *jlo = jhi;
                inc += inc;

            }
            while (x >= xx[jhi]);
#if 0
          L10:
            jhi = *jlo + inc;
            if (jhi > max__)
            {
                jhi = max__ + 1;
            }
            else if (x >= xx[jhi])
            {
                *jlo = jhi;
                inc += inc;
                goto L10;
            }
#endif
        }
        else
        {
            jhi = *jlo;

            do
            {
                *jlo = jhi - inc;

                if (*jlo <= null)
                {
                    *jlo = null;
                    break;
                }

                if (x >= xx[*jlo])
                    break;

                jhi = *jlo;
                inc += inc;
            }
            while (x < xx[*jlo]);
#if 0
          L20:
            *jlo = jhi - inc;
            if (*jlo <= null)
            {
                *jlo = null;
            }
            else if (x < xx[*jlo])
            {
                jhi = *jlo;
                inc += inc;
                goto L20;
            }
#endif
        }

    }

    while (jhi - *jlo != 1)
    {

        jm = (jhi + *jlo) / 2;

        if (x > xx[jm])
            *jlo = jm;
        else
            jhi = jm;

    }


#if 0
  L30:
    if (jhi - *jlo == 1)
    {
        return;
    }
    jm = (jhi + *jlo) / 2;
    if (x > xx[jm])
    {
        *jlo = jm;
    }
    else
    {
        jhi = jm;
    }
    goto L30;
#endif
}                               /* huntn_ */





/*This is a function that does interpolation*/
void bspline_interp_full (double * rho, double * rho_f)
{
    double *rho_c;

    int xorder, yorder, zorder, cmax_x, cmax_y, cmax_z, ix, iy, iz;
    double *bscoef, *crdsx, *crdsy, *crdsz, t1;
    double cspacing_x, cspacing_y, cspacing_z, fspacing_x, fspacing_y, fspacing_z;
    double offset_x, offset_y, offset_z;

    /*These variables should keep values between calls */
    static int flag = 0;
    static double *biatx, *biaty, *biatz, *knots_x, *knots_y, *knots_z;
    static int *leftx, *lefty, *leftz;





    /*Order of spline interpolation */
    xorder = ct.interp_order;
    yorder = ct.interp_order;
    zorder = ct.interp_order;

    /*Number of points in bigger coarse grid */
    cmax_x = (get_PX0_GRID() + 2 * ct.interp_trade);
    cmax_y = (get_PY0_GRID() + 2 * ct.interp_trade);
    cmax_z = (get_PZ0_GRID() + 2 * ct.interp_trade);

    my_malloc (bscoef, cmax_x * cmax_y * cmax_z, double);

    /*This will hold coordinates of the grid points */
    my_malloc (crdsx, get_FPX0_GRID(), double);
    my_malloc (crdsy, get_FPY0_GRID(), double);
    my_malloc (crdsz, get_FPZ0_GRID(), double);

    if (!flag)
    {

        /*Memory for knots and B-spline coefficients */
        my_malloc (knots_x, cmax_x + xorder, double);
        my_malloc (knots_y, cmax_y + yorder, double);
        my_malloc (knots_z, cmax_z + zorder, double);

        my_malloc (biatx, get_FPX0_GRID() * xorder, double);
        my_malloc (biaty, get_FPY0_GRID() * yorder, double);
        my_malloc (biatz, get_FPZ0_GRID() * zorder, double);

        my_malloc (leftx, get_FPX0_GRID(), int);
        my_malloc (lefty, get_FPY0_GRID(), int);
        my_malloc (leftz, get_FPZ0_GRID(), int);
    }



    /*Setup crds for knots */
    t1 = 1.0 / (double) (cmax_x - 1);
    for (ix = 0; ix < cmax_x; ix++)
        crdsx[ix] = (double) ix *t1;

    t1 = 1.0 / (double) (cmax_y - 1);
    for (ix = 0; ix < cmax_y; ix++)
        crdsy[ix] = (double) ix *t1;

    t1 = 1.0 / (double) (cmax_z - 1);
    for (ix = 0; ix < cmax_z; ix++)
        crdsz[ix] = (double) ix *t1;


    /*Here we do things that need to be done only once - 
     * first time we are going this function*/
    if (!flag)
    {
        /*Find knots */
        dbsnak (cmax_x, crdsx, xorder, knots_x);
        dbsnak (cmax_y, crdsy, yorder, knots_y);
        dbsnak (cmax_z, crdsz, zorder, knots_z);


        /*Make sure that knots are non-decreasing */
        for (ix = 0; ix < cmax_x + xorder - 1; ix++)
            if (knots_x[ix] > knots_x[ix + 1])
            {
                printf ("\n PE %d: knots_x[%d]=%f knots_x[%d]=%f",
                        pct.gridpe, ix, knots_x[ix], ix + 1, knots_x[ix + 1]);
                error_handler ("knots_x are not non-decreasing");
            }


        for (iy = 0; iy < cmax_y + yorder - 1; iy++)
            if (knots_y[iy] > knots_y[iy + 1])
            {
                printf ("\n PE %d: knots_y[%d]=%f knots_y[%d]=%f",
                        pct.gridpe, iy, knots_y[iy], iy + 1, knots_y[iy + 1]);
                error_handler ("knots_y are not non-decreasing");
            }


        for (iz = 0; iz < cmax_z + zorder - 1; iz++)
            if (knots_z[iz] > knots_z[iz + 1])
            {
                printf ("\n PE %d: knots_z[%d]=%f knots_z[%d]=%f",
                        pct.gridpe, iz, knots_z[iz], iz + 1, knots_z[iz + 1]);
                error_handler ("knots_z are not non-decreasing");
            }

    }                           /*end if (!flag) */




    /*Allocate memory for rho_c,  do required trade_image, get interpolation
     * coeffcients and release memory for rho_c
     * In short, do everything that depends on ct.interp_trade*/

    my_malloc (rho_c, cmax_x * cmax_y * cmax_z, double);
    trade_imagesx (rho, rho_c, get_PX0_GRID(), get_PY0_GRID(), get_PZ0_GRID(), ct.interp_trade, FULL_TRADE);
    /*Find B-spline coeffients */
    dbs3in (cmax_x, crdsx, cmax_y, crdsy, cmax_z, crdsz, rho_c,
            cmax_x, cmax_y, cmax_z, xorder, yorder, zorder, knots_x, knots_y, knots_z, bscoef);

    my_free (rho_c);



    /*if (pct.gridpe == 0) {
       printf("\n Value at 0,0,0 %e", dbs3vl(0,0,0,xorder,yorder,zorder,
       knots_x, knots_y, knots_z, cmax_x, cmax_y, cmax_z, bscoef)) ;
       printf("rho_c->b[0][0][0] %e", rho_c->b[0][0][0]);
       } */



    /*Here we do things that need to be done only once - 
     * first time we are going this function*/
    if (!flag)
    {

        /*Spacing on coarse grid */
        cspacing_x = 1.0 / (double) (cmax_x - 1);
        cspacing_y = 1.0 / (double) (cmax_y - 1);
        cspacing_z = 1.0 / (double) (cmax_z - 1);

        /*Fine grid spacing */
        fspacing_x = cspacing_x / get_FG_RATIO();
        fspacing_y = cspacing_y / get_FG_RATIO();
        fspacing_z = cspacing_z / get_FG_RATIO();

        offset_x = ct.interp_trade * cspacing_x;
        offset_y = ct.interp_trade * cspacing_y;
        offset_z = ct.interp_trade * cspacing_z;

        /*Setup coordinates that will be used to interpolate */
        for (ix = 0; ix < get_FPX0_GRID(); ix++)
            crdsx[ix] = offset_x + ix * fspacing_x;

        for (iy = 0; iy < get_FPY0_GRID(); iy++)
            crdsy[iy] = offset_y + iy * fspacing_y;

        for (iz = 0; iz < get_FPZ0_GRID(); iz++)
            crdsz[iz] = offset_z + iz * fspacing_z;


        /*Calculate arrays that depend on grid
         * These are needed when calculating interpolation values, 
         * but need to be setup only once*/
        get_biats (get_FPX0_GRID(), crdsx, get_FPY0_GRID(), crdsy, get_FPZ0_GRID(), crdsz,
                   xorder, yorder, zorder,
                   knots_x, knots_y, knots_z,
                   cmax_x, cmax_y, cmax_z, biatx, biaty, biatz, leftx, lefty, leftz);

        flag = 1;

    }                           /*end if (!flag) */


    dbs3gd2 (get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(),
             xorder, yorder, zorder,
             cmax_x, cmax_y, cmax_z, bscoef,
             rho_f, get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(), biatx, biaty, biatz, leftx, lefty, leftz);


    my_free (bscoef);
    my_free (crdsx);
    my_free (crdsy);
    my_free (crdsz);

}
