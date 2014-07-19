
#include "transition.h"
#include "const.h"
#include "BaseThread.h"
#include "TradeImages.h"
#include "RmgTimer.h"
#include "RmgThread.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "State.h"
#include "Kpoint.h"
#include "../Headers/prototypes.h"

#include "GlobalSums.h"
#include "blas.h"

void app_nls_single (Kpoint<std::complex<double>> *kptr, double * psiR, double * psiI, double * workR, double * workI, double *work2R, double *work2I, double *sintR, double *sintI, int state);


void AppNls(Kpoint<double> *kpoint, double *sintR, double *sintI)
{

    int num_states = kpoint->get_nstates();
    int P0_BASIS = Rmg_G->get_P0_BASIS(1);

    const char *transa = "n";

    int idx, ion, gion, sindex;
    int i, j, nh, inh;
    int alloc;
    double *dnmI;
    double *nwork, *psintR, *qqq;
    ION *iptr;
    SPECIES *sp;
    double rzero = 0.0, rone=1.0;
    double *sintR_compack;
    int istate, proj_index;

    double *psiR = kpoint->orbital_storage;
    double *nv = kpoint->nv;
    double *ns = kpoint->ns;
    double *Bns = kpoint->Bns;

    if(pct.num_tot_proj == 0)
    {
        for(i = 0; i < num_states * P0_BASIS; i++)
        {
            nv[i] = 0.0;
            Bns[i] = 0.0;
        }
        for(int idx = 0;idx < num_states * P0_BASIS;idx++)
            ns[idx] = psiR[idx];

        return;
    }
            
    alloc = P0_BASIS;

    alloc = pct.num_tot_proj * num_states;
    sintR_compack = new double[alloc];
    nwork = new double[alloc];

    for(i = 0; i < num_states * pct.num_tot_proj; i++)
            sintR_compack[i] = 0.0;


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


    dgemm (transa, transa, &pct.num_tot_proj, &num_states, &pct.num_tot_proj, 
            &rone, (double *)pct.M_dnm,  &pct.num_tot_proj, (double *)sintR_compack, &pct.num_tot_proj,
            &rzero,  (double *)nwork, &pct.num_tot_proj);
    dgemm (transa, transa, &P0_BASIS, &num_states, &pct.num_tot_proj, 
            &rone, (double *)pct.Bweight,  &P0_BASIS, (double *)nwork, &pct.num_tot_proj,
            &rzero,  (double *)nv, &P0_BASIS);


    for(int idx = 0;idx < num_states * P0_BASIS;idx++)
        ns[idx] = psiR[idx];

    dgemm (transa, transa, &pct.num_tot_proj, &num_states, &pct.num_tot_proj, 
            &rone, (double *)pct.M_qqq,  &pct.num_tot_proj, (double *)sintR_compack, &pct.num_tot_proj,
            &rzero,  (double *)nwork, &pct.num_tot_proj);

    if(!ct.norm_conserving_pp) {
        dgemm (transa, transa, &P0_BASIS, &num_states, &pct.num_tot_proj, 
            &rone, (double *)pct.weight,  &P0_BASIS, (double *)nwork, &pct.num_tot_proj,
            &rone,  (double *)ns, &P0_BASIS);
    }

    dgemm (transa, transa, &P0_BASIS, &num_states, &pct.num_tot_proj, 
            &rone, (double *)pct.Bweight,  &P0_BASIS, (double *)nwork, &pct.num_tot_proj,
            &rzero,  (double *)Bns, &P0_BASIS);

    delete [] nwork;
    delete [] sintR_compack;

}

void AppNls(Kpoint<std::complex<double>> *kpoint, double *sintR, double *sintI)
{

    for(int st = 0;st < kpoint->nstates;st++) {

        double *nvR = &pct.nv[2*st * kpoint->pbasis]; 
        double *nvI = nvR + kpoint->pbasis;
        double *nsR = &pct.ns[2*st * kpoint->pbasis]; 
        double *nsI = nsR + kpoint->pbasis;
        double *tmp_psiR = (double *)kpoint->Kstates[st].psi;
        double *tmp_psiI = tmp_psiR + kpoint->pbasis;

        pack_to_standard(tmp_psiR, 1, kpoint->pbasis);
        app_nls_single(kpoint, tmp_psiR, tmp_psiI, nvR, nvI, nsR, nsI, sintR, sintI, st);

        pack_to_complex(tmp_psiR, 1, kpoint->pbasis);
        pack_to_complex(nvR, 1, kpoint->pbasis);
        pack_to_complex(nsR, 1, kpoint->pbasis);
    }

}




void app_nls_single (Kpoint<std::complex<double>> *kptr, double * psiR, double * psiI, double *nvR, double *nvI, double *nsR, double *nsI, double *sintR, double *sintI, int state)
{

    int ion, gion, sindex;
    int i, j, nh, inh;
    int incx = 1, alloc;
    double *weiptr, *mptr, *dnmI, coeffR, coeffI, coeff2R, coeff2I;
    double *pR, *pI, *psintR, *qqq;
    ION *iptr;
    SPECIES *sp;

    double *psintI;


    alloc = kptr->pbasis;
    if (alloc < ct.max_nlpoints)
        alloc = ct.max_nlpoints;
    double *nworkR = new double[alloc];
    double *nworkI = new double[alloc];
    double *nwork2R = new double[alloc];
    double *nwork2I = new double[alloc];

    for(int idx = 0;idx < kptr->pbasis; idx++)
        nworkR[idx] = 0.0;
    for(int idx = 0;idx < kptr->pbasis; idx++)
        nworkI[idx] = 0.0;
    for(int idx = 0;idx < kptr->pbasis; idx++)
        nwork2R[idx] = 0.0;
    for(int idx = 0;idx < kptr->pbasis; idx++)
        nwork2I[idx] = 0.0;

    /*Base index for sintR and sintI */
    sindex = kptr->kidx * pct.num_nonloc_ions * ct.num_states * ct.max_nl + state * ct.max_nl;



    /* Zero out the work array */
    for (int idx = 0; idx <kptr->pbasis; idx++)
        nvR[idx] = 0.0;
    for (int idx = 0; idx <kptr->pbasis; idx++)
        nvI[idx] = 0.0;


    // Copy wavefunctions
    for(int idx = 0;idx < kptr->pbasis;idx++)
        nsR[idx] = psiR[idx];
    for(int idx = 0;idx < kptr->pbasis;idx++)
        nsI[idx] = psiI[idx];

    /* Loop over ions once again */
    for (ion = 0; ion < pct.num_nonloc_ions; ion++)
    {

        /*Actual index of the ion under consideration*/
        gion = pct.nonloc_ions_list[ion];

        /*This needs to be here, since nonlocal ions include those that have overlap due to either beta
         * or Q and here we only need those that overlap due to beta*/
        if (pct.idxptrlen[gion])
        {
            iptr = &ct.ions[gion];
            sp = &ct.sp[iptr->species];
       
            nh = sp->nh;

            psintR = &sintR[ion * ct.num_states * ct.max_nl + sindex];
            psintI = &sintI[ion * ct.num_states * ct.max_nl + sindex];

            //weiptr = &pct.weight[gion];
            weiptr = &pct.weight[ion * ct.max_nl * kptr->pbasis];

            dnmI = pct.dnmI[gion];
            qqq = pct.qqq[gion];

            pR = pct.phaseptr[gion];
            pR += 2 * kptr->kidx * kptr->pbasis;
            pI = pR + kptr->pbasis;

            for (i = 0; i < nh; i++)
            {
                mptr = weiptr + i * kptr->pbasis;
                coeffR = 0.0;
                coeffI = 0.0;
                coeff2R = 0.0;
                coeff2I = 0.0;
                inh = i * nh;
                for (j = 0; j < nh; j++)
                {
                    coeffR += dnmI[inh + j] * psintR[j];
                    coeff2R += qqq[inh + j] * psintR[j];
                    coeffI += dnmI[inh + j] * psintI[j];
                    coeff2I += qqq[inh + j] * psintI[j];
                }                   /* end for j */
                QMD_daxpy (kptr->pbasis, coeffR, mptr, incx, nworkR, incx);
                QMD_daxpy (kptr->pbasis, coeff2R, mptr, incx, nwork2R, incx);
                QMD_daxpy (kptr->pbasis, coeffI, mptr, incx, nworkI, incx);
                QMD_daxpy (kptr->pbasis, coeff2I, mptr, incx, nwork2I, incx);
            }                       /*end for i */

            /* Write back the results */
            for (int idx = 0; idx < kptr->pbasis; idx++)
            {
                nvR[idx] += (nworkR[idx] * pR[idx] + nworkI[idx] * pI[idx]);
                nsR[idx] += (nwork2R[idx] * pR[idx] + nwork2I[idx] * pI[idx]);
            }

            for (int idx = 0; idx < kptr->pbasis; idx++)
            {
                nvI[idx] += (-nworkR[idx] * pI[idx] + nworkI[idx] * pR[idx]);
                nsI[idx] += (-nwork2R[idx] * pI[idx] + nwork2I[idx] * pR[idx]);
            }                       /* end for */

            for (int idx = 0; idx < kptr->pbasis; idx++)
            {
                nworkI[idx] = 0.0;
                nwork2I[idx] = 0.0;
                nworkR[idx] = 0.0;
                nwork2R[idx] = 0.0;
            }

        }



    }                           /* end for */

    delete [] nwork2I;
    delete [] nwork2R;
    delete [] nworkI;
    delete [] nworkR;




}                               /* end app_nl */

