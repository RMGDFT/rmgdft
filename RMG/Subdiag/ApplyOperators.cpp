#include <complex>
#include "FiniteDiff.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_alloc.h"
#include "rmg_error.h"
#include "rmgthreads.h"
#include "RmgTimer.h"
#include "RmgThread.h"
#include "GlobalSums.h"
#include "Subdiag.h"

#include "prototypes.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"

template void ApplyOperators<double>(Kpoint<double> *, int, double *, double *, double *);
template void ApplyOperators<std::complex<double> >(Kpoint<std::complex<double>> *, int, std::complex<double> *, std::complex<double> *, double *);

// Applies A and B operators to one wavefunction
template <typename KpointType>
void ApplyOperators (Kpoint<KpointType> *kptr, int istate, KpointType *a_psi, KpointType *b_psi, double *vtot)
{
    BaseGrid *G = kptr->G;
    TradeImages *T = kptr->T;
    Lattice *L = &Rmg_L;
    STATE *sp = &kptr->kstates[istate];

    double vel = L->get_omega() / (G->get_NX_GRID(1) * G->get_NY_GRID(1) * G->get_NZ_GRID(1));

    int dimx = G->get_PX0_GRID(1);
    int dimy = G->get_PY0_GRID(1);
    int dimz = G->get_PZ0_GRID(1);
    int P0_BASIS = dimx * dimy * dimz;

    /* Generate 2*V*psi and store it in a smoothing grid and store in sg_twovpsi */
    if((ct.potential_acceleration_constant_step > 0.0) || (ct.potential_acceleration_poisson_step > 0.0)) {
        if(ct.scf_steps > 0)
        {
            if(ct.kohn_sham_fd_order == APP_CI_FOURTH)
                T->trade_imagesx (sp->dvhxc, vtot, dimx, dimy, dimz, 1, FULL_TRADE);

            if(ct.kohn_sham_fd_order == APP_CI_SIXTH)
                T->trade_imagesx (sp->dvhxc, vtot, dimx, dimy, dimz, 2, FULL_TRADE);
        }
    }

    /* A operating on psi stored in work3 */
    AppCilrDriver(T, kptr->Kstates[istate].psi, a_psi, b_psi, vtot, dimx, dimy, dimz, 
                  G->get_hxgrid(1), G->get_hygrid(1), G->get_hzgrid(1), ct.kohn_sham_fd_order);


    // Add in non-local
    for(int idx = 0; idx < P0_BASIS; idx++) {

        a_psi[idx] += TWO * kptr->nv[istate * kptr->pbasis + idx];

    }


    for (int idx = 0; idx < P0_BASIS; idx++) {

        a_psi[idx] = 0.5 * vel * a_psi[idx];

    }

    for(int idx = 0; idx < P0_BASIS; idx++) {

        b_psi[idx] += kptr->Bns[istate * kptr->pbasis + idx];

    }


} // end ApplyOperators


