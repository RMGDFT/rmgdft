#include "const.h"
#include "common_prototypes.h"
#include "State.h"
#include "Kpoint.h"
#include "BaseThread.h"
#include "TradeImages.h"
#include "RmgTimer.h"
#include "RmgThread.h"
#include "GlobalSums.h"
#include "vhartree.h"
#include "packfuncs.h"
#include "blas.h"

// Transitional headers
#include "transition.h"
#include "typedefs.h"


template <typename RmgType> void AppNls(Kpoint<RmgType> *kpoint, double *work, double *work2, double *sintR, double *sintI, double *Bns)
{

    int num_states = kpoint->get_nstates();
    int P0_BASIS = Rmg_G->get_P0_BASIS(1);

    char *transa = "n", *transt = "t";

    RmgType psiR = kpoint->orbital_storage;

#if GPU_ENABLED
    cublasStatus_t custat;
    cublasOperation_t cu_transT = CUBLAS_OP_T, cu_transN = CUBLAS_OP_N;
#endif


    if(typeid(RmgType) == typeid(double)) {

        int idx, ion, gion, stop, ip, sindex, index2, ione=1, itwo=2, P0_BASIS;
        int *pidx;
        int i, j, nh, inh;
        int incx = 1, alloc, step, count;
        double *weiptr, *mptr, *dnmI, coeffR, coeffI, coeff2R, coeff2I;
        double *nwork, *nwork2, *pR, *pI, *psintR, *qqq;
        ION *iptr;
        SPECIES *sp;
        double coeffMatR[2*MAX_NL], coeffMatI[2*MAX_NL], rzero = 0.0, rone=1.0;
        double *sintR_compack, *sintI_compack;
        int istate, proj_index;


        if(pct.num_tot_proj == 0)
        {
            for(i = 0; i < num_states * P0_BASIS; i++)
            {
                work[i] = 0.0;
                Bns[i] = 0.0;
            }
            for(int idx = 0;idx < num_states * P0_BASIS;idx++)
                work2[idx] = psiR[idx];
            return;
        }
            

        alloc =P0_BASIS;

        alloc = pct.num_tot_proj * num_states;
    #if GPU_ENABLED
        sintR_compack = ct.gpu_host_temp1;
    #else
        sintR_compack = new double[alloc];
    #endif
        nwork = new double[alloc];

        for(i = 0; i < num_states * pct.num_tot_proj; i++)
            sintR_compack[i] = 0.0;
        /*Base index for sintR and sintI */


        for(istate = 0; istate < num_states; istate++)
        {
            sindex = kpoint->kidx * pct.num_nonloc_ions * num_states * ct.max_nl + istate * ct.max_nl;
            for (ion = 0; ion < pct.num_nonloc_ions; ion++)
            {
                proj_index = ion * ct.max_nl;
                psintR = &sintR[ion * num_states * ct.max_nl + sindex];
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
        cublasSetVector( num_states * pct.num_tot_proj, sizeof( double ), sintR_compack, ione, ct.gpu_work3, ione );

        cublasSetVector( pct.num_tot_proj*pct.num_tot_proj, sizeof( double ), pct.M_dnm, ione, ct.gpu_work1, ione );
        cublasDgemm (ct.cublas_handle, cu_transN, cu_transN, pct.num_tot_proj, num_states, pct.num_tot_proj, 
                &rone, ct.gpu_work1,  pct.num_tot_proj, ct.gpu_work3, pct.num_tot_proj,
                &rzero,  ct.gpu_work2, pct.num_tot_proj);
        cublasDgemm (ct.cublas_handle, cu_transN, cu_transN, P0_BASIS, num_states, pct.num_tot_proj, 
                &rone, ct.gpu_Bweight,  P0_BASIS, ct.gpu_work2, pct.num_tot_proj,
                &rzero,  ct.gpu_temp, P0_BASIS);
        cublasGetVector( num_states * P0_BASIS, sizeof( double ), ct.gpu_temp, ione, work, ione);


        cublasSetVector( num_states * P0_BASIS, sizeof( double ), psiR, ione, ct.gpu_temp, ione);
        cublasSetVector( pct.num_tot_proj*pct.num_tot_proj, sizeof( double ), pct.M_qqq, ione, ct.gpu_work1, ione );
        cublasDgemm (ct.cublas_handle, cu_transN, cu_transN, pct.num_tot_proj, num_states, pct.num_tot_proj, 
                &rone, ct.gpu_work1,  pct.num_tot_proj, ct.gpu_work3, pct.num_tot_proj,
                &rzero,  ct.gpu_work2, pct.num_tot_proj);
        cublasDgemm (ct.cublas_handle, cu_transN, cu_transN, P0_BASIS, num_states, pct.num_tot_proj, 
                &rone, ct.gpu_weight,  P0_BASIS, ct.gpu_work2, pct.num_tot_proj,
                &rone,  ct.gpu_temp, P0_BASIS);
        cublasGetVector( num_states * P0_BASIS, sizeof( double ), ct.gpu_temp, ione, work2, ione);

        cublasDgemm (ct.cublas_handle, cu_transN, cu_transN, P0_BASIS, num_states, pct.num_tot_proj, 
                &rone, ct.gpu_Bweight,  P0_BASIS, ct.gpu_work2, pct.num_tot_proj,
                &rzero,  ct.gpu_temp, P0_BASIS);
        cublasGetVector( num_states * P0_BASIS, sizeof( double ), ct.gpu_temp, ione, Bns, ione);

    #else

        dgemm (transa, transa, &pct.num_tot_proj, &num_states, &pct.num_tot_proj, 
                &rone, pct.M_dnm,  &pct.num_tot_proj, sintR_compack, &pct.num_tot_proj,
                &rzero,  nwork, &pct.num_tot_proj);
        dgemm (transa, transa, &P0_BASIS, &num_states, &pct.num_tot_proj, 
                &rone, pct.Bweight,  &P0_BASIS, nwork, &pct.num_tot_proj,
                &rzero,  work, &P0_BASIS);


        for(int idx = 0;idx < num_states * P0_BASIS;idx++)
            work2[idx] = psiR[idx];

        dgemm (transa, transa, &pct.num_tot_proj, &num_states, &pct.num_tot_proj, 
                &rone, pct.M_qqq,  &pct.num_tot_proj, sintR_compack, &pct.num_tot_proj,
                &rzero,  nwork, &pct.num_tot_proj);
        dgemm (transa, transa, &P0_BASIS, &num_states, &pct.num_tot_proj, 
                &rone, pct.weight,  &P0_BASIS, nwork, &pct.num_tot_proj,
                &rone,  work2, &P0_BASIS);
        dgemm (transa, transa, &P0_BASIS, &num_states, &pct.num_tot_proj, 
                &rone, pct.Bweight,  &P0_BASIS, nwork, &pct.num_tot_proj,
                &rzero,  Bns, &P0_BASIS);
    #endif


        delete [] nwork;
    #if !GPU_ENABLED
        delete [] sintR_compack;
    #endif


    }
    else {


    }

}
