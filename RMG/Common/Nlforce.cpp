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
#include "const.h"
#include "RmgTimer.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "Kpoint.h"
#include "prototypes.h"


#include "TradeImages.h"
#include "RmgException.h"
#include "Lattice.h"
#include "FiniteDiff.h"
#include "transition.h"



/*Set this to 1 to write out true NL force and the part
 * that comes from eigenvalues*/
#define VERBOSE 0

template void Nlforce<double> (double *, Kpoint<double> **Kptr, double *force_nl);
template void Nlforce<std::complex<double> > (double * , Kpoint<std::complex<double>> **Kptr, double *force_nl);

template <typename OrbitalType> void Nlforce (double * veff, Kpoint<OrbitalType> **Kptr, double *force_nl)
{
    int ion, isp, index, gion, nion;
    int nh, size, size1;
    double *gamma, *par_gamma, *par_omega;
    SPECIES *sp;
    ION *iptr;
    int num_ions;
    fftw_plan p2;
    std::complex<double> *in, *out;
    double  *qforce;
    double *tmp_force_gamma, *tmp_force_omega;
    int fpt0;

    OrbitalType *gx, *gy, *gz;
    OrbitalType *psi, *psi_x, *psi_y, *psi_z;
    double hxxgrid, hyygrid, hzzgrid;
    double hxgrid, hygrid, hzgrid;

    int FPX0_GRID, FPY0_GRID, FPZ0_GRID, FP0_BASIS;
    int PX0_GRID, PY0_GRID, PZ0_GRID, P0_BASIS;
    int num_occupied;

    hxxgrid = get_hxxgrid();
    hyygrid = get_hyygrid();
    hzzgrid = get_hzzgrid();
    hxgrid = get_hxgrid();
    hygrid = get_hygrid();
    hzgrid = get_hzgrid();

    FPX0_GRID = get_FPX0_GRID();
    FPY0_GRID = get_FPY0_GRID();
    FPZ0_GRID = get_FPZ0_GRID();
    PX0_GRID = get_PX0_GRID();
    PY0_GRID = get_PY0_GRID();
    PZ0_GRID = get_PZ0_GRID();
    FP0_BASIS = FPX0_GRID * FPY0_GRID * FPZ0_GRID;



    gx = new OrbitalType[FP0_BASIS];
    gy = new OrbitalType[FP0_BASIS];
    gz = new OrbitalType[FP0_BASIS];
    RmgTimer *RT1;

    size1 = ct.num_kpts * ct.num_states * ct.num_ions * ct.max_nl;
    fpt0 = ct.fpt[0];


    num_ions = ct.num_ions;


    /*Array for q-force */
    qforce = new double[ 3 * num_ions];
    tmp_force_gamma = new double[ 3 * num_ions];
    tmp_force_omega = new double[ 3 * num_ions];

    for (isp = 0; isp < 3 * num_ions; isp++)
    {
        qforce[isp] = 0.0;
        tmp_force_gamma[isp] = 0.0;
        tmp_force_omega[isp] = 0.0;
    }


    /*max for nh * (nh + 1) / 2 */
    size = (ct.max_nl + 1) * ct.max_nl / 2;
    
    gamma = new double[ size];
    par_gamma = new double[ 6 * size];
    par_omega = par_gamma + 3 * size;


    RT1 = new RmgTimer("Force: non-local: betaxpsi");

    if (ct.force_derivate_type == WAVEFUNCTION_DERIVATIVE)
    {
        if(!ct.is_gamma)
        { printf("\n WARNING:  need more test for wavefunction derivative in force calculation for non-gamma point calculation");
            fflush(NULL);
            exit(0);
        } 
        for (int kpt = 0; kpt < ct.num_kpts; kpt++)
        {

            num_occupied = 0;
            for(int st = 0; st < ct.num_states; st++)
            {
                if(Kptr[kpt]->Kstates[st].occupation[0] < 1.0e-10) break;
                num_occupied++;
                psi = Kptr[kpt]->Kstates[st].psi;
                psi_x = Kptr[kpt]->Kstates[st + ct.num_states].psi;
                psi_y = Kptr[kpt]->Kstates[st + 2*ct.num_states].psi;
                psi_z = Kptr[kpt]->Kstates[st + 3*ct.num_states].psi;
                CPP_app_grad_driver (&Rmg_L, Rmg_T, psi, psi_x, psi_y, psi_z, PX0_GRID, PY0_GRID, PZ0_GRID, hxgrid, hygrid, hzgrid, 6);
            }


            Betaxpsi(Kptr[kpt], 1*Kptr[kpt]->nstates, num_occupied, Kptr[kpt]->sint_derx, Kptr[kpt]->nl_weight);
            Betaxpsi(Kptr[kpt], 2*Kptr[kpt]->nstates, num_occupied, Kptr[kpt]->sint_dery, Kptr[kpt]->nl_weight);
            Betaxpsi(Kptr[kpt], 3*Kptr[kpt]->nstates, num_occupied, Kptr[kpt]->sint_derz, Kptr[kpt]->nl_weight);

            for(int i = 0; i < pct.num_nonloc_ions * ct.num_states * ct.max_nl; i++)
            {
                Kptr[kpt]->sint_derx[i] *= -1.0;
                Kptr[kpt]->sint_dery[i] *= -1.0;
                Kptr[kpt]->sint_derz[i] *= -1.0;
            }

        }
    }
    else
    {

        for (int kpt = 0; kpt < ct.num_kpts; kpt++)
        {
            Betaxpsi(Kptr[kpt], 0, Kptr[kpt]->nstates, Kptr[kpt]->sint_derx, Kptr[kpt]->nl_weight_derx);
            Betaxpsi(Kptr[kpt], 0, Kptr[kpt]->nstates, Kptr[kpt]->sint_dery, Kptr[kpt]->nl_weight_dery);
            Betaxpsi(Kptr[kpt], 0, Kptr[kpt]->nstates, Kptr[kpt]->sint_derz, Kptr[kpt]->nl_weight_derz);
        }
    }
    delete RT1;


    RT1 = new RmgTimer("Force: non-local: veff grad");
    CPP_app_grad_driver (&Rmg_L, Rmg_T, veff, (double *)gx, (double *)gy, (double *)gz, FPX0_GRID, FPY0_GRID, FPZ0_GRID, hxxgrid, hyygrid, hzzgrid, 6);
    delete RT1;


    for (ion = 0; ion < pct.num_nonloc_ions; ion++)
    {
        /*Actual index of the ion under consideration*/
        gion = pct.nonloc_ions_list[ion];

        iptr = &ct.ions[gion];

        nh = ct.sp[iptr->species].nh;

        RT1 = new RmgTimer("Force: non-local: get gamma");
        GetGamma (gamma, ion, nh, Kptr);
        delete RT1;
        RT1 = new RmgTimer("Force: non-local: nlforce_par_Q");
        nlforce_par_Q ((double *)gx, (double *)gy, (double *)gz, gamma, gion, iptr, nh, &qforce[3 * gion]);
        delete RT1;

    }                           /*end for(ion=0; ion<ions_max; ion++) */

    delete [] gx;
    delete [] gy;
    delete [] gz;


    for(int i = 0; i < ct.num_ions * 3; i++) 
    {
        qforce[i] *= get_vel_f();
    }







    /*Loop over ions again */
    nion = -1;
    for (ion = 0; ion < pct.num_owned_ions; ion++)
    {
        /*Global index of owned ion*/
        gion = pct.owned_ions_list[ion];

        /* Figure out index of owned ion in nonloc_ions_list array, store it in nion*/
        do {

            nion++;
            if (nion >= pct.num_nonloc_ions)
            {
                printf("\n Could not find matching entry in pct.nonloc_ions_list for owned ion %d", gion);
                rmg_error_handler(__FILE__, __LINE__, "Could not find matching entry in pct.nonloc_ions_list for owned ion ");
            }

        } while (pct.nonloc_ions_list[nion] != gion);

        iptr = &ct.ions[gion];


        nh = ct.sp[iptr->species].nh;

        /*partial_gamma(ion,par_gamma,par_omega, iptr, nh, p1, p2); */
        RT1 = new RmgTimer("Force: non-local: partial gamma");
        PartialGamma (gion, par_gamma, par_omega, nion, nh, Kptr);
        delete RT1;
        RT1 = new RmgTimer("Force: non-local: nlforce_par_gamma");
        nlforce_par_gamma (par_gamma, gion, nh, &tmp_force_gamma[3*gion]);
        delete RT1;

        RT1 = new RmgTimer("Force: non-local: nlforce_par_omega");
        nlforce_par_omega (par_omega, gion, nh, &tmp_force_omega[3*gion]);
        delete RT1;

    }                           /*end for(ion=0; ion<num_ions; ion++) */


    
    for(int i = 0; i < ct.num_ions * 3; i++) 
    {
        force_nl[i] += qforce[i];
        force_nl[i] += tmp_force_gamma[i];
        force_nl[i] += tmp_force_omega[i];
    }


#if VERBOSE  
    output_force(qforce, "Non-local forces: QForce");
    output_force(tmp_force_gamma, "Non-local forces: der_gamma term");
    output_force(tmp_force_omega, "Non-local forces: der_omega term");
#endif


    delete[] par_gamma;
    delete[] gamma;
    delete[] tmp_force_gamma;
    delete[] tmp_force_omega;
    delete[] qforce;


}

