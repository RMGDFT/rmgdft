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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#ifdef USE_NUMA
    #include <numa.h>
#endif
#include <sys/mman.h>
#include "transition.h"
#include "const.h"
#include "RmgTimer.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "Kpoint.h"
#include "Subdiag.h"
#include "Functional.h"
#include "GpuAlloc.h"
#include "ErrorFuncs.h"
#include "RmgException.h"
#include "Functional.h"
#include "Solvers.h"
#include "Atomic.h"
#include "RmgParallelFft.h"
#include "FiniteDiff.h"
#include "Scalapack.h"
#include "Elpa.h"
#include "GatherScatter.h"



template <typename OrbitalType> void Reinit (double * vh, double * rho, double * rho_oppo, double * rhocore, double * rhoc,
        double * vnuc, double * vxc,  Kpoint<OrbitalType> **Kptr);
template void Reinit<double>(double*, double*, double*, double*, double*, double*, double*, Kpoint<double>**);
template void Reinit<std::complex<double> >(double*, double*, double*, double*, double*, double*, double*, Kpoint<std::complex <double> >**);

template <typename OrbitalType> void Reinit (double * vh, double * rho, double * rho_oppo, double * rhocore, double * rhoc,
        double * vnuc, double * vxc,  Kpoint<OrbitalType> **Kptr)
{

    RmgTimer RT0("ReInit");
    int kpt;
    int species;

    // after cell relax, change atom's position from crystal to cartesian

    for(size_t ion= 0; ion < Atoms.size(); ion++)
    {
        Rmg_L.to_cartesian(Atoms[ion].xtal, Atoms[ion].crds);
    }

    // Initialize some commonly used plans for our parallel ffts
    FftInitPlans();

    // If ecutwfc was set then adjust filter factor
    if(ct.ecutwfc > 0.0)
    {
        double tpiba2 = 4.0 * PI * PI / (Rmg_L.celldm[0] * Rmg_L.celldm[0]);
        double ecut = coarse_pwaves->gmax * tpiba2;

        if(ecut < ct.ecutwfc)
        {
            rmg_printf("WARNING: The value of ecutwfc you have selected is to large for the specified grid.  %7.2f %7.2f\n", ct.ecutwfc, ecut);
        }
        else
        {
            coarse_pwaves->gcut = ct.ecutwfc / tpiba2;
            fine_pwaves->gcut = coarse_pwaves->gcut * Rmg_G->default_FG_RATIO * Rmg_G->default_FG_RATIO;
        }

    }

    ct.hmaxgrid = Rmg_L.get_xside() * Rmg_G->get_hxgrid(1);
    if (Rmg_L.get_yside() * Rmg_G->get_hygrid(1) > ct.hmaxgrid)
        ct.hmaxgrid = Rmg_L.get_yside() * Rmg_G->get_hygrid(1);
    if (Rmg_L.get_zside() * Rmg_G->get_hzgrid(1) > ct.hmaxgrid)
        ct.hmaxgrid = Rmg_L.get_zside() * Rmg_G->get_hzgrid(1);

    if(ct.ecutrho <= 0.0) ct.ecutrho = (2.0 *PI/ct.hmaxgrid) *(2.0 *PI/ct.hmaxgrid);
    ct.hmingrid = Rmg_L.get_xside() * Rmg_G->get_hxgrid(1);
    if (Rmg_L.get_yside() * Rmg_G->get_hygrid(1) < ct.hmingrid)
        ct.hmingrid = Rmg_L.get_yside() * Rmg_G->get_hygrid(1);
    if (Rmg_L.get_zside() * Rmg_G->get_hzgrid(1) < ct.hmingrid)
        ct.hmingrid = Rmg_L.get_zside() * Rmg_G->get_hzgrid(1);

    double density[3];
    double planar_anisotropy = GetPlanarAnisotropy(density);
    if (planar_anisotropy > 1.1)
    {
        if (pct.imgpe == 0)
        {
	    printf ("\n");
	    printf ("Coordinate planes\n");
	    printf ("  Planar Anisotropy: %5.3f\n", planar_anisotropy);
	    printf ("  A0-A1 density: %9.3f\n", density[0]);
	    printf ("  A0-A2 density: %9.3f\n", density[1]);
	    printf ("  A1-A2 density: %9.3f\n", density[2]);
            printf ("Planar Anisotropy is large: %f", planar_anisotropy);
        }
    }


//  rescale the kpoint with correct lattice vectors.
    for (int kpt = 0; kpt < ct.num_kpts; kpt++) {
        double v1, v2, v3;
        v1 = 0.0;
        v2 = 0.0;
        v3 = 0.0;

        for(int ir = 0; ir<3; ir++)
        {
            v1 = ct.kp[kpt].kpt[0] *Rmg_L.b0[0]
                + ct.kp[kpt].kpt[1] *Rmg_L.b1[0] 
                + ct.kp[kpt].kpt[2] *Rmg_L.b2[0];
            v2 = ct.kp[kpt].kpt[0] *Rmg_L.b0[1]
                + ct.kp[kpt].kpt[1] *Rmg_L.b1[1] 
                + ct.kp[kpt].kpt[2] *Rmg_L.b2[1];
            v3 = ct.kp[kpt].kpt[0] *Rmg_L.b0[2]
                + ct.kp[kpt].kpt[1] *Rmg_L.b1[2] 
                + ct.kp[kpt].kpt[2] *Rmg_L.b2[2];
        }

        ct.kp[kpt].kvec[0] = v1 * twoPI;
        ct.kp[kpt].kvec[1] = v2 * twoPI;
        ct.kp[kpt].kvec[2] = v3 * twoPI;
        ct.kp[kpt].kmag = (v1 * v1 + v2 * v2 + v3 * v3) * twoPI * twoPI;

    }

    for (int kpt = 0; kpt < ct.num_kpts_pe; kpt++) {
        int kpt1 = kpt + pct.kstart;
        Kptr[kpt]->kp.kvec[0] = ct.kp[kpt1].kvec[0];
        Kptr[kpt]->kp.kvec[1] = ct.kp[kpt1].kvec[1];
        Kptr[kpt]->kp.kvec[2] = ct.kp[kpt1].kvec[2];
    }



    /*Set max_nldim */
    ct.max_nldim = 0;
    for (species = 0; species < ct.num_species; species++)
    {

        /* Get species type */
        SPECIES *sp = &Species[species];
        ct.max_nldim = std::max(ct.max_nldim, sp->nldim);


        if(ct.localize_projectors)
        {
            // Grid object local to this MPI process
            delete sp->OG;
            sp->OG = new BaseGrid(sp->nldim, sp->nldim, sp->nldim, 1, 1, 1, 0, 1);
            BaseGrid *OG = (BaseGrid *)sp->OG;
            OG->set_rank(0, pct.my_comm);

            // Lattice object for the localized projector region. Global vectors need to be
            // scaled so that they correspond to the local region.
            Lattice *L = new Lattice();
            double a0[3], a1[3], a2[3], s1, celldm[6], omega;
            s1 = (double)sp->nldim / (double)Rmg_G->get_NX_GRID(1);
            a0[0] = s1*Rmg_L.get_a0(0);a0[1] = s1*Rmg_L.get_a0(1);a0[2] = s1*Rmg_L.get_a0(2);
            s1 = (double)sp->nldim / (double)Rmg_G->get_NY_GRID(1);
            a1[0] = s1*Rmg_L.get_a1(0);a1[1] = s1*Rmg_L.get_a1(1);a1[2] = s1*Rmg_L.get_a1(2);
            s1 = (double)sp->nldim / (double)Rmg_G->get_NZ_GRID(1);
            a2[0] = s1*Rmg_L.get_a2(0);a2[1] = s1*Rmg_L.get_a2(1);a2[2] = s1*Rmg_L.get_a2(2);

            L->set_ibrav_type(None);
            L->latgen(celldm, &omega, a0, a1, a2, true);

            delete sp->prj_pwave;
            sp->prj_pwave = new Pw(*OG, *L, 1, false);
            delete L;
        }
        else
        {
            delete sp->prj_pwave;
            sp->prj_pwave = new Pw(*Rmg_G, Rmg_L, 1, false);
        }

    }


    ct.max_nldim = 0;
    for(auto &sp : Species) ct.max_nldim = std::max(ct.max_nldim, sp.nldim);

    /*Do forward transform for each species and store results on the coarse grid */
    RmgTimer *RT1 = new RmgTimer("2-ReInit: weights");
    for(auto& sp : Species) sp.InitWeights (ct.localize_projectors);
    for(auto& sp : Species) sp.InitOrbitals (ct.atomic_orbital_type);
    delete(RT1);

    bool potential_acceleration = (ct.potential_acceleration_constant_step > 1.0e-10);

    // Zero out dvh array if potential acceleration is enabled
    if(potential_acceleration)
    {

        for (kpt = 0; kpt < ct.num_kpts_pe; kpt++)
        {
            int stop = Kptr[kpt]->ndvh * Kptr[kpt]->pbasis * pct.coalesce_factor;
            for(int i=0;i < stop;i++) Kptr[kpt]->dvh[i] = 0.0;
        }
    }

    // Reinitialize autotuning finite differencing
    std::unordered_map<LaplacianCoeff *, LaplacianCoeff *> T;
    for(auto& a:FiniteDiff::FdCoeffs)
    {
        T.insert_or_assign(a.second, a.second);
    }
    FiniteDiff::FdCoeffs.clear();
    for(auto& a:T) delete a.second;
    SetLaplacian();
    GetFdFactor();
    // Kpoint version
    //for (int kpt =0; kpt < ct.num_kpts_pe; kpt++) Kptr[kpt]->GetFdFactor();

    //OptimizeFdCoeff();

}                               /* end init */


