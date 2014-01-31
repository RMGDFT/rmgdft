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
 *   void app_nl(rmg_double_t *psiR, rmg_double_t *psiI, rmg_double_t *workR, rmg_double_t *workI, 
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

#define FAST_NLS 1


void app_nls (rmg_double_t * psiR, rmg_double_t * psiI, rmg_double_t * workR, rmg_double_t * workI, rmg_double_t *work2R, rmg_double_t *work2I, rmg_double_t *sintR, rmg_double_t *sintI, int state,
        int kidx)
{

    int idx, ion, gion, stop, ip, sindex, index2, ione=1, itwo=2, P0_BASIS;
    int *pidx;
    int i, j, nh, inh;
    int incx = 1, alloc, step, count;
    rmg_double_t *weiptr, *mptr, *dnmI, coeffR, coeffI, coeff2R, coeff2I;
    rmg_double_t *nworkR, *nworkI, *nwork2R, *nwork2I, *pR, *pI, *psintR, *qqq;
    ION *iptr;
    SPECIES *sp;
    rmg_double_t coeffMatR[2*MAX_NL], coeffMatI[2*MAX_NL], rzero = 0.0, rone=1.0;
    char *transa = "n";

#if !GAMMA_PT
    rmg_double_t *psintI;
#endif

    P0_BASIS = get_P0_BASIS();
    alloc = P0_BASIS;
    if (alloc < ct.max_nlpoints)
        alloc = ct.max_nlpoints;
    my_calloc (nworkR, 4 * alloc, rmg_double_t);
    nworkI = nworkR + alloc;
    nwork2R = nworkI + alloc;
    nwork2I = nwork2R + alloc;

    /*Base index for sintR and sintI */
    sindex = kidx * pct.num_nonloc_ions * ct.num_states * ct.max_nl + state * ct.max_nl;



    /* Zero out the work array */
    for (idx = 0; idx < P0_BASIS; idx++)
        workR[idx] = 0.0;

    my_copy(psiR, work2R, P0_BASIS);

#if !GAMMA_PT
    for (idx = 0; idx < P0_BASIS; idx++)
        workI[idx] = 0.0;

    my_copy(psiI, work2I, P0_BASIS);
#endif
    stop = P0_BASIS;

    /* Loop over ions once again */
    weiptr = pct.weight;
    for (ion = 0; ion < pct.num_nonloc_ions; ion++)
    {

        /*Actual index of the ion under consideration*/
        gion = pct.nonloc_ions_list[ion];
        iptr = &ct.ions[gion];
        sp = &ct.sp[iptr->species];

        nh = sp->nh;

        /*This needs to be here, since nonlocal ions include those that have overlap due to either beta
         * or Q and here we only need those that overlap due to beta*/
        if (pct.idxptrlen[gion])
        {


            psintR = &sintR[ion * ct.num_states * ct.max_nl + sindex];
#if !GAMMA_PT
            psintI = &sintI[ion * ct.num_states * ct.max_nl + sindex];
#endif

            dnmI = pct.dnmI[gion];
            qqq = pct.qqq[gion];

            pR = pct.phaseptr[gion];
            pR += 2 * kidx * stop;
            pI = pR + stop;
            pidx = pct.nlindex[gion];


            // Set this up as a pair of blas level 2 calls
            for (i = 0; i < nh; i++)
            {
                coeffMatR[i] = 0.0;         // First row of left hand matrix
                coeffMatR[nh + i] = 0.0;    // Second row of left hand matrix
                inh = i * nh;
                for (j = 0; j < nh; j++)
                {
                    coeffMatR[i]  += dnmI[inh + j] * psintR[j];
                    coeffMatR[nh + i] += qqq[inh + j] * psintR[j];
#if !GAMMA_PT
                    coeffMatI[i]  += dnmI[inh + j] * psintR[j];
                    coeffMatI[nh + i] += qqq[inh + j] * psintR[j];
#endif
                }
            }

            dgemv(transa, &stop, &nh, &rone, weiptr, &stop, coeffMatR, &ione, &rzero, nworkR, &ione);
            dgemv(transa, &stop, &nh, &rone, weiptr, &stop, &coeffMatR[nh], &ione, &rzero, nwork2R, &ione);
#if !GAMMA_PT
            dgemv(transa, &stop, &nh, &rone, weiptr, &stop, coeffMatI, &ione, &rzero, nworkI, &ione);
            dgemv(transa, &stop, &nh, &rone, weiptr, &stop, &coeffMatI[nh], &ione, &rzero, nwork2I, &ione);
#endif



#if GAMMA_PT
            /* Write back the results */
            for (idx = 0; idx < stop; idx++)
            {
                workR[idx] += nworkR[idx];
                work2R[idx] += nwork2R[idx];
            }                       /* end for */
#else

            /* Write back the results */
            for (idx = 0; idx < stop; idx++)
            {
                workR[pidx[idx]] += (nworkR[idx] * pR[idx] + nworkI[idx] * pI[idx]);
                work2R[pidx[idx]] += (nwork2R[idx] * pR[idx] + nwork2I[idx] * pI[idx]);
            }

            for (idx = 0; idx < stop; idx++)
            {
                workI[pidx[idx]] += (-nworkR[idx] * pI[idx] + nworkI[idx] * pR[idx]);
                work2I[pidx[idx]] += (-nwork2R[idx] * pI[idx] + nwork2I[idx] * pR[idx]);
            }                       /* end for */
#endif

        }

        weiptr += ct.max_nl * P0_BASIS;


    }                           /* end for */

    my_free (nworkR);


}                               /* end app_nl */

/******/
