
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
#include "vhartree.h"
#include "packfuncs.h"
#include "blas.h"

template <typename OrbitalType> void AppNls(Kpoint<OrbitalType> *kpoint, double *sintR, double *sintI);

template void AppNls<double>(Kpoint<double> *, double *, double *);
template void AppNls<std::complex<double> >(Kpoint<std::complex<double>> *, double *, double *);

template <typename OrbitalType> void AppNls(Kpoint<OrbitalType> *kpoint, double *sintR, double *sintI)
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
    OrbitalType *sintR_compack, *sintI_compack;
    int istate, proj_index;

    OrbitalType *psiR = kpoint->orbital_storage;
    OrbitalType *work = kpoint->nv;
    OrbitalType *work2 = kpoint->ns;
    OrbitalType *Bns = kpoint->Bns;

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
            
    alloc = P0_BASIS;

    alloc = pct.num_tot_proj * num_states;
    sintR_compack = new OrbitalType[alloc];
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


    if(typeid(OrbitalType) == typeid(double)) {

        dgemm (transa, transa, &pct.num_tot_proj, &num_states, &pct.num_tot_proj, 
                &rone, (double *)pct.M_dnm,  &pct.num_tot_proj, (double *)sintR_compack, &pct.num_tot_proj,
                &rzero,  (double *)nwork, &pct.num_tot_proj);
        dgemm (transa, transa, &P0_BASIS, &num_states, &pct.num_tot_proj, 
                &rone, (double *)pct.Bweight,  &P0_BASIS, (double *)nwork, &pct.num_tot_proj,
                &rzero,  (double *)work, &P0_BASIS);


        for(int idx = 0;idx < num_states * P0_BASIS;idx++)
            work2[idx] = psiR[idx];

        dgemm (transa, transa, &pct.num_tot_proj, &num_states, &pct.num_tot_proj, 
                &rone, (double *)pct.M_qqq,  &pct.num_tot_proj, (double *)sintR_compack, &pct.num_tot_proj,
                &rzero,  (double *)nwork, &pct.num_tot_proj);

        if(!ct.norm_conserving_pp) {
            dgemm (transa, transa, &P0_BASIS, &num_states, &pct.num_tot_proj, 
                &rone, (double *)pct.weight,  &P0_BASIS, (double *)nwork, &pct.num_tot_proj,
                &rone,  (double *)work2, &P0_BASIS);
        }

        dgemm (transa, transa, &P0_BASIS, &num_states, &pct.num_tot_proj, 
                &rone, (double *)pct.Bweight,  &P0_BASIS, (double *)nwork, &pct.num_tot_proj,
                &rzero,  (double *)Bns, &P0_BASIS);

    }
    else if(typeid(OrbitalType) == typeid(std::complex<double>)) {


    }


        delete [] nwork;
        delete [] sintR_compack;


}
