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



void app_nls_allstates (REAL * psiR, REAL * psiI, REAL * workR, REAL * workI, REAL *work2R, REAL *work2I, REAL *sintR, REAL *sintI, int kidx)
{

    int idx, ion, gion, stop, ip, sindex, index2, ione=1, itwo=2;
    int *pidx;
    int i, j, nh, inh;
    int incx = 1, alloc, step, count;
    REAL *weiptr, *mptr, *dnmI, coeffR, coeffI, coeff2R, coeff2I;
    REAL *nworkR, *nworkI, *nwork2R, *nwork2I, *pR, *pI, *psintR, *qqq;
    ION *iptr;
    SPECIES *sp;
    REAL coeffMatR[2*MAX_NL], coeffMatI[2*MAX_NL], rzero = 0.0, rone=1.0;
    char *transa = "n", *transt = "t";
    REAL *sintR_compack, *sintI_compack;
    int istate, proj_index;

#if !GAMMA_PT
    error_handler("\n app_nls_allstate is not programed for kpoint");
#endif

    if(pct.num_tot_proj == 0)
    {
        for(i = 0; i < ct.num_states * pct.P0_BASIS; i++)
            workR[i] = 0.0;
        my_copy(psiR, work2R, ct.num_states * pct.P0_BASIS);
        return;
    }
        

    alloc =pct.P0_BASIS;

    alloc = pct.num_tot_proj * ct.num_states;
    my_calloc (sintR_compack, alloc, REAL);
    my_calloc (nworkR, alloc, REAL);

    for(i = 0; i < ct.num_states * pct.num_tot_proj; i++)
        sintR_compack[i] = 0.0;
    /*Base index for sintR and sintI */


    for(istate = 0; istate < ct.num_states; istate++)
    {
        sindex = kidx * pct.num_nonloc_ions * ct.num_states * ct.max_nl + istate * ct.max_nl;
        for (ion = 0; ion < pct.num_nonloc_ions; ion++)
        {
            proj_index = ion * ct.max_nl;
            psintR = &sintR[ion * ct.num_states * ct.max_nl + sindex];
            /*Actual index of the ion under consideration*/
            gion = pct.nonloc_ions_list[ion];
            iptr = &ct.ions[gion];
            sp = &ct.sp[iptr->species];

            nh = sp->nh;
            for (i = 0; i < nh; i++)
            {
                sintR_compack[istate * pct.num_tot_proj + proj_index + i] =
                    psintR[i]; 
            }
        }
    }

    for (i = 0; i < pct.num_tot_proj * pct.num_tot_proj; i++)
    {
        pct.M_dnm[i] = 0.0;
        pct.M_qqq[i] = 0.0;
    }


    // set up pct.M_qqq and pct.M_dnm, this can be done outside in the
    // init.c or get_ddd get_qqq, we need to check the order
    proj_index = 0;
    for (ion = 0; ion < pct.num_nonloc_ions; ion++)
    {

        /*Actual index of the ion under consideration*/
        proj_index = ion * ct.max_nl;
        gion = pct.nonloc_ions_list[ion];
        iptr = &ct.ions[gion];
        sp = &ct.sp[iptr->species];

        nh = sp->nh;

        dnmI = pct.dnmI[gion];
        qqq = pct.qqq[gion];

        for (i = 0; i < nh; i++)
        {
            inh = i * nh;
            for (j = 0; j < nh; j++)
            {

                idx = (proj_index + i) * pct.num_tot_proj + proj_index + j;
                pct.M_dnm[idx] = dnmI[inh+j];
                pct.M_qqq[idx] = qqq[inh+j];
            }
        }
    }


    dgemm (transa, transa, &pct.num_tot_proj, &ct.num_states, &pct.num_tot_proj, 
            &rone, pct.M_dnm,  &pct.num_tot_proj, sintR_compack, &pct.num_tot_proj,
            &rzero,  nworkR, &pct.num_tot_proj);
    dgemm (transa, transa, &pct.P0_BASIS, &ct.num_states, &pct.num_tot_proj, 
            &rone, pct.weight,  &pct.P0_BASIS, nworkR, &pct.num_tot_proj,
            &rzero,  workR, &pct.P0_BASIS);



    my_copy(psiR, work2R, ct.num_states * pct.P0_BASIS);

    dgemm (transa, transa, &pct.num_tot_proj, &ct.num_states, &pct.num_tot_proj, 
            &rone, pct.M_qqq,  &pct.num_tot_proj, sintR_compack, &pct.num_tot_proj,
            &rzero,  nworkR, &pct.num_tot_proj);
    dgemm (transa, transa, &pct.P0_BASIS, &ct.num_states, &pct.num_tot_proj, 
            &rone, pct.weight,  &pct.P0_BASIS, nworkR, &pct.num_tot_proj,
            &rone,  work2R, &pct.P0_BASIS);



    my_free (nworkR);
    my_free (sintR_compack);


}                               /* end app_nl */

/******/
