/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "transition.h"
#include "const.h"
#include "State.h"
#include "Kpoint.h"
#include "TradeImages.h"
#include "RmgTimer.h"
#include "RmgThread.h"
#include "GlobalSums.h"
#include "rmgthreads.h"
#include "vhartree.h"
#include "packfuncs.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "Kpoint.h"
#include "Subdiag.h"
#include "Functional.h"
#include "Solvers.h"
#include "../../RMG/Headers/prototypes.h"
#include "RmgParallelFft.h"

#include "blas.h"
#include "../Headers/prototypes_tddft.h"



template void RmgTddft<double> (double *, double *, double *,
          double *, double *, double *, double *, Kpoint<double> **);
template void RmgTddft<std::complex<double> > (double *, double *, double *,
          double *, double *, double *, double *, Kpoint<std::complex<double>> **);
template <typename OrbitalType> void RmgTddft (double * vxc, double * vh, double * vnuc, 
        double * rho, double * rho_oppo, double * rhocore, double * rhoc, Kpoint<OrbitalType> **Kptr)
{

    double t3;
    double *vtot, *vtot_psi, *new_rho;
    double time1;
    double t[3];                  /* SCF checks and average potential */
    double max_unocc_res = 0.0;

    int dimx = Rmg_G->get_PX0_GRID(Rmg_G->get_default_FG_RATIO());
    int dimy = Rmg_G->get_PY0_GRID(Rmg_G->get_default_FG_RATIO());
    int dimz = Rmg_G->get_PZ0_GRID(Rmg_G->get_default_FG_RATIO());
    int FP0_BASIS = dimx * dimy * dimz;

    int n2, numst, P0_BASIS,i, ione =1;
    int pre_steps;
    int Ieldyn = 1, iprint = 0;

    /* to hold the send data and receive data of eigenvalues */
    double *rho_tot=NULL;   

    time1 = my_crtc ();

    P0_BASIS =  Rmg_G->get_P0_BASIS(1);
    FP0_BASIS = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);

    numst = ct.num_states; 
    n2 = numst * numst;

    double *Hmatrix = new double[n2];
    double *Hmatrix_old = new double[n2];
    double *Smatrix = new double[n2];
    double *Akick = new double[n2];
    double *Pn0 = new double[2*n2];
    double *Pn1 = new double[2*n2];
    double *vh_old = new double[FP0_BASIS];
    double *vxc_old = new double[FP0_BASIS];
    double *vh_corr_old = new double[FP0_BASIS];
    double *vh_corr = new double[FP0_BASIS];



    // Loop over k-points
    if(ct.num_kpts != 1) 
    {
        rmg_printf(" \n  TDDFT does not support multiple k-points \n");
        fflush(NULL);
        exit(0);
    }




    /* allocate memory for eigenvalue send array and receive array */

    new_rho = new double[FP0_BASIS];
    vtot = new double[FP0_BASIS];
    vtot_psi = new double[P0_BASIS];

    for (int idx = 0; idx < FP0_BASIS; idx++) vtot[idx] = 0.0;
    init_efield(vtot);
    GetVtotPsi (vtot_psi, vtot, Rmg_G->default_FG_RATIO);
    HmatrixUpdate(Kptr[0], vtot_psi, (OrbitalType *)Akick);

    /* save old vhxc + vnuc */
    for (int idx = 0; idx < FP0_BASIS; idx++) {
        vtot[idx] = vxc[idx] + vh[idx] + vnuc[idx];
    }

    // Transfer vtot from the fine grid to the wavefunction grid
    GetVtotPsi (vtot_psi, vtot, Rmg_G->default_FG_RATIO);

    /*Generate the Dnm_I */
    get_ddd (vtot);

    Betaxpsi (Kptr[0]);
    HSmatrix (Kptr[0], vtot_psi, (OrbitalType *)Hmatrix, (OrbitalType *)Smatrix);

    dcopy(&n2, Hmatrix, &ione, Hmatrix_old, &ione);



    pre_steps = 0;
    if(pre_steps == 0)
        for(i = 0; i < n2; i++) Hmatrix[i] += Akick[i];

    for(i = 0; i < 2* n2; i++) Pn0[i] = 0.0;

    for(i = 0; i < ct.nel/2; i++) Pn0[i * numst + i] = 2.0;

    for(i = 0; i < numst; i++) Kptr[0]->Kstates[i].occupation[0] = Pn0[i*numst + i]; 

    for(ct.scf_steps = 0; ct.scf_steps < ct.max_scf_steps; ct.scf_steps++)
    {

        RmgTimer *RT2a = new RmgTimer("1-TOTAL: ELDYN");
        eldyn_(&numst, Smatrix, Hmatrix, Pn0, Pn1, &Ieldyn, &iprint);
        delete(RT2a);
        for(i = 0; i < numst; i++) Kptr[0]->Kstates[i].occupation[0] = Pn0[i*numst + i]; 
        GetNewRho(Kptr, rho);

        dcopy(&FP0_BASIS, vh_corr, &ione, vh_corr_old, &ione);
        dcopy(&FP0_BASIS, vh, &ione, vh_old, &ione);
        dcopy(&FP0_BASIS, vxc, &ione, vxc_old, &ione);

        get_vxc(rho, rho_oppo, rhocore, vxc);
        VhDriver(rho_tot, rhoc, vh, vh_corr);
        for (int idx = 0; idx < FP0_BASIS; idx++) {
            vtot[idx] = vxc[idx] + vh[idx] +vh_corr[idx]
                -vxc_old[idx] -vh_old[idx] - vh_corr[idx];
        }

        GetVtotPsi (vtot_psi, vtot, Rmg_G->default_FG_RATIO);
        HmatrixUpdate(Kptr[0], vtot_psi, (OrbitalType *)Hmatrix);

        for(i = 0; i < n2; i++) Hmatrix[i] += Hmatrix_old[i];
        dcopy(&n2, Hmatrix, &ione, Hmatrix_old, &ione);

    }


}
