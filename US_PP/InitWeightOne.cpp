/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "make_conf.h"
#include "const.h"
#include "grid.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include <complex>
#include "Kpoint.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"
#include "AtomicInterpolate.h"

static double harmonic_func_00(double r, double *b);
static double harmonic_func_10(double r, double *b);
static double harmonic_func_11(double r, double *b);
static double harmonic_func_12(double r, double *b);
static double harmonic_func_20(double r, double *b);
static double harmonic_func_21(double r, double *b);
static double harmonic_func_22(double r, double *b);
static double harmonic_func_23(double r, double *b);
static double harmonic_func_24(double r, double *b);
static double harmonic_func_30(double r, double *b);
static double harmonic_func_31(double r, double *b);
static double harmonic_func_32(double r, double *b);
static double harmonic_func_33(double r, double *b);
static double harmonic_func_34(double r, double *b);
static double harmonic_func_35(double r, double *b);
static double harmonic_func_36(double r, double *b);
void InitWeightOne (SPECIES * sp, fftw_complex * rtptr, int ip, int l, int m, fftw_plan p1)
{

    int idx, ix, iy, iz, size, ibegin, iend;
    double r, ax[3], xc, yc, zc, bx[3];
    double t1, hxx, hyy, hzz;
    std::complex<double> *weptr, *gwptr;
    double (*func) (double, double *) = NULL;

    switch(l)
    {
        case S_STATE: 
            func = harmonic_func_00;
            break;
        case P_STATE:
            if(m == 0) func = harmonic_func_10;
            if(m == 1) func = harmonic_func_12;
            if(m == 2) func = harmonic_func_11;
            break;
        case D_STATE:
            if(m == 0) func = harmonic_func_20;
            if(m == 1) func = harmonic_func_21;
            if(m == 2) func = harmonic_func_22;
            if(m == 3) func = harmonic_func_23;
            if(m == 4) func = harmonic_func_24;
            break;
//        case F_STATE:
//            if(m == 0) func = harmonic_func_30;
//            if(m == 1) func = harmonic_func_31;
//            if(m == 2) func = harmonic_func_32;
//            if(m == 3) func = harmonic_func_33;
//            if(m == 4) func = harmonic_func_34;
//            if(m == 5) func = harmonic_func_35;
//            if(m == 6) func = harmonic_func_36;
//            break;
        default:
            printf("\n projecotr with ip = %d not programed \n", ip);
            exit(0);
    }




    // define functions to distiguish s, px, py, pz, ....



    /* nlfdim is size of the non-local box in the double grid */
    size = sp->nlfdim * sp->nlfdim * sp->nlfdim;

    weptr = new std::complex<double>[size];
    gwptr = new std::complex<double>[size];

    hxx = get_hxgrid() / (double) ct.nxfgrid;
    hyy = get_hygrid() / (double) ct.nyfgrid;
    hzz = get_hzgrid() / (double) ct.nzfgrid;

    /*We assume that ion is in the center of non-local box */
    //    ibegin = -(sp->nldim / 2) * ct.nxfgrid;

    ibegin = -sp->nlfdim/2;
    iend = ibegin + sp->nlfdim;
    int ixx, iyy, izz;
    for (ix = ibegin; ix < iend; ix++)
    {
        ixx = ix;
        if (ixx < 0) ixx = ix + sp->nlfdim;
        xc = (double) ix *hxx;

        for (iy = ibegin; iy < iend; iy++)
        {
            iyy = iy;
            if (iyy < 0) iyy = iy + sp->nlfdim;
            yc = (double) iy *hyy;

            for (iz = ibegin; iz < iend; iz++)
            {

                izz = iz;
                if (izz < 0) izz = iz + sp->nlfdim;
                zc = (double) iz *hzz;

                ax[0] = xc;
                ax[1] = yc;
                ax[2] = zc;

                r = metric (ax);
                to_cartesian (ax, bx);
                r += 1.0e-10;

                t1 = AtomicInterpolateInline(&sp->betalig[ip][0], r);
                idx = ixx *sp->nlfdim * sp->nlfdim + iyy * sp->nlfdim + izz;
                weptr[idx] = func(r, bx) * t1;

                if((ix*2 + sp->nlfdim) == 0 || (iy*2 + sp->nlfdim) == 0 || (iz*2 + sp->nlfdim) == 0 ) 
                    weptr[idx] = 0.0;

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

    fftw_execute_dft (p1, reinterpret_cast<fftw_complex*>(weptr), reinterpret_cast<fftw_complex*>(gwptr));

    pack_gftoc (sp, reinterpret_cast<fftw_complex*>(gwptr), reinterpret_cast<fftw_complex*>(rtptr));

    delete []gwptr;
    delete []weptr;
}

static double harmonic_func_00(double r, double *b)
{
    return sqrt(1.0/(4.0 * PI));
}

static double harmonic_func_10(double r, double *b)
{
    return sqrt(3.0/(4.0 * PI)) * b[0]/r;
}


static double harmonic_func_11(double r, double *b)
{
    return sqrt(3.0/(4.0 * PI)) * b[1]/r;
}

static double harmonic_func_12(double r, double *b)
{
    return sqrt(3.0/(4.0 * PI)) * b[2]/r;
}

static double harmonic_func_20(double r, double *b)
{
    return sqrt(15.0/(4.0 * PI)) * b[0] * b[1]/(r*r);
}

static double harmonic_func_21(double r, double *b)
{
    return sqrt(15.0/(4.0 * PI)) * b[0] * b[2]/(r*r);
}

static double harmonic_func_22(double r, double *b)
{
    return sqrt(15.0/(4.0 * PI)) * b[1] * b[2]/(r*r);
}

static double harmonic_func_23(double r, double *b)
{
    return sqrt(15.0/(4.0 * PI))* (b[0] * b[0] - b[1]*b[1])/(2.0*r*r);
}

static double harmonic_func_24(double r, double *b)
{
    return sqrt(5.0/(4.0 * PI))* (3.0 * b[2] * b[2] - r*r)/(2.0*r*r);
}


