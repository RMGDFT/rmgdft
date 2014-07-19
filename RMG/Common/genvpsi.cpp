/************************** SVN Revision Information **************************
 **    $Id: genvpsi.c 2012 2013-05-13 17:22:57Z ebriggs $    **
******************************************************************************/

/****f* QMD-MGDFT/genvpsi.c *****
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
 *   void genvpsi(double *psi, double *sg_twovpsi, double *vtot, double *vnl, double *kd,
 *                double kmag, int dimx, int dimy, int dimz)
 *   Used to generate the product 2*V*psi for the Mehrstellen B operator.
 * INPUTS
 *   psi: wave function
 *   vtot: total potential 
 *   vnl:  results of the non-local operator acting on wave function
 *   kd:   k . (del)psi
 *   kmag: |k|^2
 *   dimx, dimy, dimz: dimentions of the array
 * OUTPUT
 *   sg_twovpsi: see equation (10) of 
 *             Briggs, Sullivan and Bernholc, PRB54, 14362 (1996).
 * PARENTS
 *   mg_eig_state.c subdiag_mpi.c subdiag_smp.c
 * CHILDREN
 *   nothing
 * SOURCE
 */

#include <complex>
#include "make_conf.h"
#include "rmgtypes.h"
#include "common_prototypes.h"
#include "transition.h"

#define TWO 2.0


// Gamma point float version
void CPP_genvpsi (float * psi, float * sg_twovpsi, double * vtot, double * vnl, void * kdp,
              double kmag, int dimx, int dimy, int dimz)
{

    int ix, iy, iz;
    int incx, incy;

    double *kd = (double *)kdp;

    incy = dimz;
    incx = (dimy) * (dimz);

    /* Generate 2 * V * psi */
    for (ix = 0; ix < dimx; ix++)
    {

        for (iy = 0; iy < dimy; iy++)
        {

            for (iz = 0; iz < dimz; iz++)
            {
                sg_twovpsi[(ix) * incx + (iy) * incy + iz] =
                    TWO * psi[ix * incx + iy * incy + iz] *
                    vtot[ix * incx + iy * incy + iz] + TWO * vnl[ix * incx + iy * incy + iz];
            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


}                               /* end genvpsi */

// Gamma point double version
void CPP_genvpsi (double * psi, double * sg_twovpsi, double * vtot, double * vnl, double * kdp,
              double kmag, int dimx, int dimy, int dimz)
{

    int ix, iy, iz;
    int incx, incy;

    double *kd = (double *)kdp;

    incy = dimz;
    incx = (dimy) * (dimz);

    /* Generate 2 * V * psi */
    for (ix = 0; ix < dimx; ix++)
    {

        for (iy = 0; iy < dimy; iy++)
        {

            for (iz = 0; iz < dimz; iz++)
            {
                sg_twovpsi[ix * incx + iy * incy + iz] =
                    TWO * psi[ix * incx + iy * incy + iz] *
                    vtot[ix * incx + iy * incy + iz] + TWO * vnl[ix * incx + iy * incy + iz];
            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

}

// complex float version
void CPP_genvpsi (std::complex<float> * psi, std::complex<float> * sg_twovpsi, double * vtot, std::complex<double> * vnl, void * kdp,
              double kmag, int dimx, int dimy, int dimz)
{

    int ix, iy, iz;
    int incx, incy;

    std::complex<double> *kd = (std::complex<double> *)kdp;

    incy = dimz;
    incx = (dimy) * (dimz);

    /* Generate 2 * V * psi */
    for (ix = 0; ix < dimx; ix++)
    {

        for (iy = 0; iy < dimy; iy++)
        {

            for (iz = 0; iz < dimz; iz++)
            {
                sg_twovpsi[ix * incx + iy * incy + iz] = (std::complex<float>)
                    (
                    2.0 * (std::complex<double>)psi[ix * incx + iy * incy + iz] *
                    (vtot[ix * incx + iy * incy + iz] + 0.5 * kmag) +
                    2.0 * kd[ix * incx + iy * incy + iz] + TWO * vnl[ix * incx + iy * incy + iz]
                    );
            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

}                               /* end genvpsi */


// complex double version
void CPP_genvpsi (std::complex<double> * psi, std::complex<double> * sg_twovpsi, double * vtot, std::complex<double> * vnl, void * kdp,
              double kmag, int dimx, int dimy, int dimz)
{

    int ix, iy, iz;
    int incx, incy;

    std::complex<double> *kd = (std::complex<double> *)kdp;

    incy = dimz;
    incx = (dimy) * (dimz);

    /* Generate 2 * V * psi */
    for (ix = 0; ix < dimx; ix++)
    {

        for (iy = 0; iy < dimy; iy++)
        {

            for (iz = 0; iz < dimz; iz++)
            {
                sg_twovpsi[ix * incx + iy * incy + iz] = (std::complex<double>)
                    (
                    2.0 * (std::complex<double>)psi[ix * incx + iy * incy + iz] *
                    (vtot[ix * incx + iy * incy + iz] + 0.5 * kmag) +
                    2.0 * kd[ix * incx + iy * incy + iz] + TWO * vnl[ix * incx + iy * incy + iz]
                    );
            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

}                               /* end genvpsi */


extern "C" void genvpsi (double * psi, double * sg_twovpsi, double * vtot, double * vnl, double * kd,
              double kmag, int dimx, int dimy, int dimz)
{

    CPP_genvpsi(psi, sg_twovpsi, vtot, vnl, kd, kmag, dimx, dimy, dimz);

}
extern "C" void genvpsi_f (float * psi, float * sg_twovpsi, double * vtot, double * vnl, double * kd,
              double kmag, int dimx, int dimy, int dimz)
{

    CPP_genvpsi(psi, sg_twovpsi, vtot, vnl, kd, kmag, dimx, dimy, dimz);

}
