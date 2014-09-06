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
    KpointType *psi = kptr->Kstates[istate].psi;

    double vel = L->get_omega() / (G->get_NX_GRID(1) * G->get_NY_GRID(1) * G->get_NZ_GRID(1));
    double hxgrid = G->get_hxgrid(1);
    double hygrid = G->get_hygrid(1);
    double hzgrid = G->get_hzgrid(1);

    int dimx = G->get_PX0_GRID(1);
    int dimy = G->get_PY0_GRID(1);
    int dimz = G->get_PZ0_GRID(1);
    int pbasis = dimx * dimy * dimz;
    int sbasis = (dimx + 4) * (dimy + 4) * (dimz + 4);

    bool potential_acceleration = ((ct.potential_acceleration_constant_step > 0.0) || (ct.potential_acceleration_poisson_step > 0.0));
    potential_acceleration = potential_acceleration & (ct.scf_steps > 0);


    // Apply A operator to psi
    CPP_app_cil_driver (L, T, psi, a_psi, dimx, dimy, dimz, hxgrid, hygrid, hzgrid, ct.kohn_sham_fd_order);


    // Apply B operator to psi
    CPP_app_cir_driver (L, T, psi, b_psi, dimx, dimy, dimz, ct.kohn_sham_fd_order);


    // if complex orbitals apply gradient to orbital and compute dot products
    std::complex<double> *kdr = new std::complex<double>[pbasis]();

    if(typeid(KpointType) == typeid(std::complex<double>)) {

        KpointType *gx = new KpointType[pbasis];
        KpointType *gy = new KpointType[pbasis];
        KpointType *gz = new KpointType[pbasis];

        CPP_app_grad_driver (L, T, psi, gx, gy, gz, dimx, dimy, dimz, hxgrid, hygrid, hzgrid, ct.kohn_sham_fd_order);

        std::complex<double> I_t(0.0, 1.0);
        for(int idx = 0;idx < pbasis;idx++) {

            kdr[idx] = -I_t * (kptr->kvec[0] * (std::complex<double>)gx[idx] +
                               kptr->kvec[1] * (std::complex<double>)gy[idx] +
                               kptr->kvec[2] * (std::complex<double>)gz[idx]);
        }

        delete [] gz;
        delete [] gy;
        delete [] gx;

    }

    // Generate 2*V*psi
    KpointType *sg_twovpsi_t = new KpointType[sbasis];
    if(potential_acceleration) {
        CPP_genvpsi (psi, sg_twovpsi_t, sp->dvhxc, (void *)kdr, kptr->kmag, dimx, dimy, dimz);
    }
    else {
        CPP_genvpsi (psi, sg_twovpsi_t, vtot, (void *)kdr, kptr->kmag, dimx, dimy, dimz);
    }

    // B operating on 2*V*psi stored in work
    KpointType *work_t = new KpointType[sbasis];
    CPP_app_cir_driver (L, T, sg_twovpsi_t, work_t, dimx, dimy, dimz, ct.kohn_sham_fd_order);

    for(int idx = 0; idx < pbasis; idx++) {

        a_psi[idx] = work_t[idx] - a_psi[idx];

    }

    // Add in non-local which has already had B applied in AppNls
    for(int idx = 0; idx < pbasis; idx++) {

        a_psi[idx] += TWO * kptr->nv[istate * kptr->pbasis + idx];

    }

    for (int idx = 0; idx < pbasis; idx++) {

        a_psi[idx] = 0.5 * vel * a_psi[idx];

    }

    // Add in already applied Bns to b_psi
    for(int idx = 0; idx < pbasis; idx++) {

        b_psi[idx] += kptr->Bns[istate * kptr->pbasis + idx];

    }

    delete [] kdr;
    delete [] work_t;
    delete [] sg_twovpsi_t;

} // end ApplyOperators


