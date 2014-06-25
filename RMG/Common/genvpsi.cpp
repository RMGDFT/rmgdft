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
 *   void genvpsi(rmg_double_t *psi, rmg_double_t *sg_twovpsi, rmg_double_t *vtot, rmg_double_t *vnl, rmg_double_t *kd,
 *                rmg_double_t kmag, int dimx, int dimy, int dimz)
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


#include "make_conf.h"
#include "rmgtypes.h"
#include "common_prototypes.h"
#include "auxiliary.h"

#define TWO 2.0


template void CPP_genvpsi<double>(double*, double*, double*, double*, double*, double, int, int, int);
template void CPP_genvpsi<float>(float*, float*, double*, double*, double*, double, int, int, int);

template <typename RmgType>
void CPP_genvpsi (RmgType * psi, RmgType * sg_twovpsi, rmg_double_t * vtot, rmg_double_t * vnl, rmg_double_t * kd,
              rmg_double_t kmag, int dimx, int dimy, int dimz)
{

    int ix, iy, iz;
    int incx, incy;
    int incx1, incy1;

    incy = dimz;
    incx = (dimy) * (dimz);

    incy1 = dimz;
    incx1 = dimy * dimz;

    /* Generate 2 * V * psi */
    for (ix = 0; ix < dimx; ix++)
    {

        for (iy = 0; iy < dimy; iy++)
        {

            for (iz = 0; iz < dimz; iz++)
            {
#if GAMMA_PT
                sg_twovpsi[(ix) * incx + (iy) * incy + iz] =
                    TWO * psi[ix * incx1 + iy * incy1 + iz] *
                    vtot[ix * incx1 + iy * incy1 + iz];

#else
                sg_twovpsi[(ix) * incx + (iy) * incy + iz] =
                    TWO * psi[ix * incx1 + iy * incy1 + iz] *
                    (vtot[ix * incx1 + iy * incy1 + iz] + 0.5 * kmag) +
                    TWO * kd[ix * incx1 + iy * incy1 + iz] +
                    TWO * vnl[ix * incx1 + iy * incy1 + iz];
#endif

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


}                               /* end genvpsi */


extern "C" void genvpsi (rmg_double_t * psi, rmg_double_t * sg_twovpsi, rmg_double_t * vtot, rmg_double_t * vnl, rmg_double_t * kd,
              rmg_double_t kmag, int dimx, int dimy, int dimz)
{

    CPP_genvpsi<double>(psi, sg_twovpsi, vtot, vnl, kd, kmag, dimx, dimy, dimz);

}
extern "C" void genvpsi_f (rmg_float_t * psi, rmg_float_t * sg_twovpsi, rmg_double_t * vtot, rmg_double_t * vnl, rmg_double_t * kd,
              rmg_double_t kmag, int dimx, int dimy, int dimz)
{

    CPP_genvpsi<float>(psi, sg_twovpsi, vtot, vnl, kd, kmag, dimx, dimy, dimz);

}
