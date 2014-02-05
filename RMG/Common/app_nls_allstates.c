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
#include "common_prototypes.h"


#if GPU_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif



void app_nls_allstates (rmg_double_t * psiR, rmg_double_t * psiI, rmg_double_t * workR, rmg_double_t * workI, 
rmg_double_t *work2R, rmg_double_t *work2I, rmg_double_t *Bns, rmg_double_t *BnsI, rmg_double_t *sintR, rmg_double_t *sintI, int kidx)
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
    char *transa = "n", *transt = "t";
    rmg_double_t *sintR_compack, *sintI_compack;
    int istate, proj_index;
#if GPU_ENABLED
    cublasStatus_t custat;
    cublasOperation_t cu_transT = CUBLAS_OP_T, cu_transN = CUBLAS_OP_N;
#endif

#if !GAMMA_PT
    error_handler("\n app_nls_allstate is not programed for kpoint");
#endif

    P0_BASIS = get_P0_BASIS();

    if(pct.num_tot_proj == 0)
    {
        for(i = 0; i < ct.num_states * P0_BASIS; i++)
        {
            workR[i] = 0.0;
            Bns[i] = 0.0;
        }
        my_copy(psiR, work2R, ct.num_states * P0_BASIS);
        return;
    }
        

    alloc =P0_BASIS;

    alloc = pct.num_tot_proj * ct.num_states;
#if GPU_ENABLED
    sintR_compack = ct.gpu_host_temp1;
#else
    my_calloc (sintR_compack, alloc, rmg_double_t);
#endif
    my_calloc (nworkR, alloc, rmg_double_t);

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

#if GPU_ENABLED
    cublasSetVector( ct.num_states * pct.num_tot_proj, sizeof( rmg_double_t ), sintR_compack, ione, ct.gpu_work3, ione );

    cublasSetVector( pct.num_tot_proj*pct.num_tot_proj, sizeof( rmg_double_t ), pct.M_dnm, ione, ct.gpu_work1, ione );
    cublasDgemm (ct.cublas_handle, cu_transN, cu_transN, pct.num_tot_proj, ct.num_states, pct.num_tot_proj, 
            &rone, ct.gpu_work1,  pct.num_tot_proj, ct.gpu_work3, pct.num_tot_proj,
            &rzero,  ct.gpu_work2, pct.num_tot_proj);
    cublasDgemm (ct.cublas_handle, cu_transN, cu_transN, P0_BASIS, ct.num_states, pct.num_tot_proj, 
            &rone, ct.gpu_Bweight,  P0_BASIS, ct.gpu_work2, pct.num_tot_proj,
            &rzero,  ct.gpu_temp, P0_BASIS);
    cublasGetVector( ct.num_states * P0_BASIS, sizeof( rmg_double_t ), ct.gpu_temp, ione, workR, ione);


    cublasSetVector( ct.num_states * P0_BASIS, sizeof( rmg_double_t ), psiR, ione, ct.gpu_temp, ione);
    cublasSetVector( pct.num_tot_proj*pct.num_tot_proj, sizeof( rmg_double_t ), pct.M_qqq, ione, ct.gpu_work1, ione );
    cublasDgemm (ct.cublas_handle, cu_transN, cu_transN, pct.num_tot_proj, ct.num_states, pct.num_tot_proj, 
            &rone, ct.gpu_work1,  pct.num_tot_proj, ct.gpu_work3, pct.num_tot_proj,
            &rzero,  ct.gpu_work2, pct.num_tot_proj);
    cublasDgemm (ct.cublas_handle, cu_transN, cu_transN, P0_BASIS, ct.num_states, pct.num_tot_proj, 
            &rone, ct.gpu_weight,  P0_BASIS, ct.gpu_work2, pct.num_tot_proj,
            &rone,  ct.gpu_temp, P0_BASIS);
    cublasGetVector( ct.num_states * P0_BASIS, sizeof( rmg_double_t ), ct.gpu_temp, ione, work2R, ione);

    cublasDgemm (ct.cublas_handle, cu_transN, cu_transN, P0_BASIS, ct.num_states, pct.num_tot_proj, 
            &rone, ct.gpu_Bweight,  P0_BASIS, ct.gpu_work2, pct.num_tot_proj,
            &rzero,  ct.gpu_temp, P0_BASIS);
    cublasGetVector( ct.num_states * P0_BASIS, sizeof( rmg_double_t ), ct.gpu_temp, ione, Bns, ione);

#else

    dgemm (transa, transa, &pct.num_tot_proj, &ct.num_states, &pct.num_tot_proj, 
            &rone, pct.M_dnm,  &pct.num_tot_proj, sintR_compack, &pct.num_tot_proj,
            &rzero,  nworkR, &pct.num_tot_proj);
    dgemm (transa, transa, &P0_BASIS, &ct.num_states, &pct.num_tot_proj, 
            &rone, pct.Bweight,  &P0_BASIS, nworkR, &pct.num_tot_proj,
            &rzero,  workR, &P0_BASIS);


    my_copy(psiR, work2R, ct.num_states * P0_BASIS);

    dgemm (transa, transa, &pct.num_tot_proj, &ct.num_states, &pct.num_tot_proj, 
            &rone, pct.M_qqq,  &pct.num_tot_proj, sintR_compack, &pct.num_tot_proj,
            &rzero,  nworkR, &pct.num_tot_proj);
    dgemm (transa, transa, &P0_BASIS, &ct.num_states, &pct.num_tot_proj, 
            &rone, pct.weight,  &P0_BASIS, nworkR, &pct.num_tot_proj,
            &rone,  work2R, &P0_BASIS);
    dgemm (transa, transa, &P0_BASIS, &ct.num_states, &pct.num_tot_proj, 
            &rone, pct.Bweight,  &P0_BASIS, nworkR, &pct.num_tot_proj,
            &rzero,  Bns, &P0_BASIS);
#endif


    my_free (nworkR);
#if !GPU_ENABLED
    my_free (sintR_compack);
#endif


}                               /* end app_nl */

/******/
