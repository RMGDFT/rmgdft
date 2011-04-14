/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

/****f* QMD-MGDFT/app_nl.c *****
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
 *   void app_nl(REAL *psiR, REAL *psiI, REAL *workR, REAL *workI, 
 *               int state, int flag, int kidx, int tid)
 *    Applies the non-local potential operator to an orbital.
 *    Also mixes the projections of the non-local operator on the
 *    orbitals from step to step.
 * INPUTS
 *   psiR: Points to the real part of the orbital
 *   psiI: Points to the imaginary part of the orbital
 *   state: Index of the orbital from zero
 *   kidx: K-point that this orbital belongs to
 *   flag: a flag to control whether we want to mix the projections.
 *   tid: Thread ID, only used for timing purposes when running 
 *        in shared memory mode.
 * OUTPUT
 *   workR: Real part of the Non-local potential applied to orbital
 *   workI: Imaginary part of the Non-local potential applied to orbital
 * PARENTS
 *   get_milliken.c mg_eig_state.c subdiag_mpi.c subdiag_smp.c
 * CHILDREN
 *   global_sums.c
 * SOURCE
 */


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"


void app_nl_eig (REAL * psiR, REAL * psiI, REAL * workR, REAL * workI, REAL *sintR, REAL *sintI, int state,
                 int kidx)
{

    int idx, ion, stop, ip, sindex, index2;
    int *pidx;
    int i, j, nh;
    int incx = 1, alloc, step, count;
    REAL *weiptr, *mptr, *dnmI, coeffR, coeffI;
    REAL *nworkR, *nworkI, *pR, *pI, *psintR;
    ION *iptr;

    #if !GAMMA_PT
      REAL *psintI;
    #endif



    alloc = P0_BASIS;
    if (alloc < ct.max_nlpoints)
        alloc = ct.max_nlpoints;
    my_calloc (nworkR, 2 * alloc, REAL);
    nworkI = nworkR + alloc;

    /*Base index for sintR and sintI */
    sindex = kidx * ct.num_ions * ct.num_states * ct.max_nl + state * ct.max_nl;



    /* Zero out the work array */
    for (idx = 0; idx < P0_BASIS; idx++)
        workR[idx] = 0.0;
#if !GAMMA_PT
    for (idx = 0; idx < P0_BASIS; idx++)
        workI[idx] = 0.0;
#endif



    /*Simplified version of the loop below */
#if 0
    /* Loop over ions once again */
    for (ion = 0; ion < ct.num_ions; ion++)
    {

        stop = pct.idxptrlen[ion];
        if (!stop)
            continue;

        psintR = &sintR[ion * ct.num_states * ct.max_nl + sindex];
#if !GAMMA_PT
        psintI = &sintI[ion * ct.num_states * ct.max_nl + sindex];
#endif

        weiptr = pct.weight[ion];
        nh = pct.prj_per_ion[ion];
        dnmI = pct.dnmI[ion];

        pR = pct.phaseptr[ion];
        pR += 2 * kidx * stop;
        pI = pR + stop;
        pidx = pct.nlindex[ion];


        for (i = 0; i < nh; i++)
        {
            mptr = weiptr + i * stop;
            coeffR = 0.0;
            coeffI = 0.0;
            inh = i * nh;
            for (j = 0; j < nh; j++)
            {
                coeffR += dnmI[inh + j] * psintR[j];
#if !GAMMA_PT
                coeffI += dnmI[inh + j] * psintI[j];
#endif
            }                   /* end for j */
            saxpy (&stop, &coeffR, mptr, &incx, nworkR, &incx);
#if !GAMMA_PT
            saxpy (&stop, &coeffI, mptr, &incx, nworkI, &incx);
#endif
        }                       /*end for i */

#if GAMMA_PT
        /* Write back the results */
        for (idx = 0; idx < stop; idx++)
        {
            workR[pidx[idx]] += nworkR[idx];
            nworkR[idx] = 0.0;

        }                       /* end for */
#else

        /* Write back the results */
        for (idx = 0; idx < stop; idx++)
            workR[pidx[idx]] += (nworkR[idx] * pR[idx] + nworkI[idx] * pI[idx]);

        for (idx = 0; idx < stop; idx++)
        {
            workI[pidx[idx]] += (-nworkR[idx] * pI[idx] + nworkI[idx] * pR[idx]);
            nworkR[idx] = nworkI[idx] = 0.0;

        }                       /* end for */
#endif



    }                           /* end for */
#else



    /*Atomized version of the loop above, it seems to give slightly better performance,
     * but it is quite difficult to read and understand*/

    /* Loop over ions once again */
    for (ion = 0; ion < ct.num_ions; ion++)
    {

        stop = pct.idxptrlen[ion];
        if (!stop)
            continue;

        iptr = &ct.ions[ion];
        weiptr = pct.weight[ion];
        nh = pct.prj_per_ion[ion];
        dnmI = pct.dnmI[ion];

        psintR = &sintR[ion * ct.num_states * ct.max_nl + sindex];
#if !GAMMA_PT
        psintI = &sintI[ion * ct.num_states * ct.max_nl + sindex];
#endif

        pR = pct.phaseptr[ion];
        pR += 2 * kidx * stop;
        pI = pR + stop;


        count = 0;
        step = 128;
        pidx = pct.nlindex[ion];
        while (count < stop)
        {

            if ((step + count) > stop)
                step = stop - count;
            /* Now apply the non-local operator to the wavefunctions */
            for (idx = 0; idx < step; idx++)
                nworkR[idx] = ZERO;
#if !GAMMA_PT
            for (idx = 0; idx < step; idx++)
                nworkI[idx] = ZERO;
#endif

            if (step)
            {
                for (i = 0; i < nh; i++)
                {
                    mptr = weiptr + i * stop;
                    coeffR = 0.0;
                    coeffI = 0.0;
                    for (j = 0; j < nh; j++)
                    {
                        coeffR += dnmI[i * nh + j] * psintR[j];
#if !GAMMA_PT
                        coeffI += dnmI[i * nh + j] * psintI[j];
#endif
                    }           /* end for j */
                    QMD_saxpy (step, coeffR, &mptr[count], incx, nworkR, incx);
#if !GAMMA_PT
                    QMD_saxpy (step, coeffI, &mptr[count], incx, nworkI, incx);
#endif
                }               /*end for i */
            }

#if GAMMA_PT
            /* Write back the results */
            for (idx = 0; idx < step; idx++)
                workR[pidx[idx + count]] += nworkR[idx];
#else

            /* Write back the results */
            for (idx = 0; idx < step; idx++)
                workR[pidx[idx + count]] += (nworkR[idx] * pR[idx + count] +
                                             nworkI[idx] * pI[idx + count]);

            for (idx = 0; idx < step; idx++)
                workI[pidx[idx + count]] += (-nworkR[idx] * pI[idx + count] +
                                             nworkI[idx] * pR[idx + count]);
#endif


            count += step;
        }                       /* end while */

    }                           /* end for */
#endif

    my_free (nworkR);


}                               /* end app_nl */

/******/
