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
#include "GpuAlloc.h"


#include "TradeImages.h"
#include "RmgException.h"
#include "Lattice.h"
#include "FiniteDiff.h"
#include "transition.h"
#include "GlobalSums.h"



/*Set this to 1 to write out true NL force and the part
 * that comes from eigenvalues*/
#define VERBOSE 0

template void Nlforce<double> (double *, Kpoint<double> **Kptr, double *force_nl);
template void Nlforce<std::complex<double> > (double * , Kpoint<std::complex<double>> **Kptr, double *force_nl);

template <typename OrbitalType> void Nlforce (double * veff, Kpoint<OrbitalType> **Kptr, double *force_nl)
{
    int ion, isp, gion, nion;
    int nh, size;
    double *gamma, *par_gamma, *par_omega;
    ION *iptr;
    int num_ions;
    double  *qforce;
    double *tmp_force_gamma, *tmp_force_omega;

    OrbitalType *psi, *psi_x, *psi_y, *psi_z;
    double hxxgrid, hyygrid, hzzgrid;

    int FPX0_GRID, FPY0_GRID, FPZ0_GRID, FP0_BASIS;
    int PX0_GRID, PY0_GRID, PZ0_GRID, P0_BASIS;
    int num_occupied;
    std::complex<double> I_t(0.0, 1.0);

    hxxgrid = get_hxxgrid();
    hyygrid = get_hyygrid();
    hzzgrid = get_hzzgrid();

    FPX0_GRID = get_FPX0_GRID();
    FPY0_GRID = get_FPY0_GRID();
    FPZ0_GRID = get_FPZ0_GRID();
    PX0_GRID = get_PX0_GRID();
    PY0_GRID = get_PY0_GRID();
    PZ0_GRID = get_PZ0_GRID();
    FP0_BASIS = FPX0_GRID * FPY0_GRID * FPZ0_GRID;
    P0_BASIS = PX0_GRID * PY0_GRID * PZ0_GRID;


    RmgTimer *RT1;

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
    int max_nl2 = (ct.max_nl + 1) * ct.max_nl / 2;
    
    gamma = new double[ max_nl2];
//    par_gamma = new double[ 6 * max_nl2];
//    par_omega = par_gamma + 3 * max_nl2;
    double *par_gamma_allions = new double[3 * pct.num_owned_ions * max_nl2];
    double *par_omega_allions = new double[3 * pct.num_owned_ions * max_nl2];
    for(int i = 0; i < 3 * pct.num_owned_ions * max_nl2; i++)
    {
        par_gamma_allions[i] = 0.0;
        par_omega_allions[i] = 0.0;
    }

    
    size =  pct.num_nonloc_ions * ct.state_block_size * ct.max_nl; 
#if GPU_ENABLED
    OrbitalType *sint_der = (OrbitalType *)GpuMallocHost(3*size * sizeof(OrbitalType));
    OrbitalType *sint_derx = sint_der + 0 * size;
    OrbitalType *sint_dery = sint_der + 1 * size;
    OrbitalType *sint_derz = sint_der + 2 * size;
#else
    OrbitalType *sint_der = new OrbitalType[3*size];
    OrbitalType *sint_derx = sint_der + 0 * size;
    OrbitalType *sint_dery = sint_der + 1 * size;
    OrbitalType *sint_derz = sint_der + 2 * size;
#endif
    


//  determin the number of occupied states for all kpoints.

    num_occupied = 0;
    for (int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
    {
        for(int st = 0; st < ct.num_states; st++)
        {
            if(abs(Kptr[kpt]->Kstates[st].occupation[0]) < 1.0e-10) 
            {
                num_occupied = std::max(num_occupied, st);
                break;
            }
        }
    }


    int num_state_block = (num_occupied +ct.state_block_size -1) / ct.state_block_size;
    int *state_start = new int[num_state_block];
    int *state_end = new int[num_state_block];
    for(int ib = 0; ib < num_state_block; ib++)
    {
        state_start[ib] = ib * ct.state_block_size;
        state_end[ib] = (ib+1) * ct.state_block_size;
        if(state_end[ib] > num_occupied) state_end[ib] = num_occupied;
    }

    int pbasis = Kptr[0]->pbasis;


    if(ct.alloc_states < ct.num_states + 3 * ct.state_block_size)     
    {
       printf("\n allco_states %d should be larger than ct.num_states %d + 3* ct.state_block_size, %d", ct.alloc_states, ct.num_states,
ct.state_block_size);
       exit(0);
    }
    
    for (int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
    {

        for(int ib = 0; ib < num_state_block; ib++)
        {
            for(int st = state_start[ib]; st < state_end[ib]; st++)
            {
                psi = Kptr[kpt]->Kstates[st].psi;
                psi_x = Kptr[kpt]->Kstates[ct.num_states].psi + (st-state_start[ib]) * pbasis;
                psi_y = psi_x + ct.state_block_size*pbasis;
                psi_z = psi_x +2* ct.state_block_size*pbasis;
                ApplyGradient(psi, psi_x, psi_y, psi_z, ct.force_grad_order, "Coarse");

                if(!ct.is_gamma)
                {
                    std::complex<double> *psi_C, *psi_xC, *psi_yC, *psi_zC;
                    psi_C = (std::complex<double> *) psi;
                    psi_xC = (std::complex<double> *) psi_x;
                    psi_yC = (std::complex<double> *) psi_y;
                    psi_zC = (std::complex<double> *) psi_z;
                    for(int i = 0; i < P0_BASIS; i++) 
                    {
                        psi_xC[i] += I_t *  Kptr[kpt]->kvec[0] * psi_C[i];
                        psi_yC[i] += I_t *  Kptr[kpt]->kvec[1] * psi_C[i];
                        psi_zC[i] += I_t *  Kptr[kpt]->kvec[2] * psi_C[i];
                    }
                }


            }


            int num_state_thisblock = state_end[ib] - state_start[ib];


            RT1 = new RmgTimer("2-Force: non-local-betaxpsi");
            Betaxpsi(Kptr[kpt], ct.num_states,                       num_state_thisblock, sint_derx, Kptr[kpt]->nl_weight);
            Betaxpsi(Kptr[kpt], ct.num_states+ ct.state_block_size,  num_state_thisblock, sint_dery, Kptr[kpt]->nl_weight);
            Betaxpsi(Kptr[kpt], ct.num_states+2*ct.state_block_size, num_state_thisblock, sint_derz, Kptr[kpt]->nl_weight);

            for(int i = 0; i < pct.num_nonloc_ions * num_state_thisblock * ct.max_nl; i++)
            {
                sint_derx[i] *= -1.0;
                sint_dery[i] *= -1.0;
                sint_derz[i] *= -1.0;
            }
            delete RT1;

            RT1 = new RmgTimer("2-Force: non-local-partial gamma");
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
                par_gamma = &par_gamma_allions[ion * 3 * max_nl2];
                par_omega = &par_omega_allions[ion * 3 * max_nl2];

                PartialGamma (gion, par_gamma, par_omega, nion, nh, Kptr, state_start[ib], state_end[ib], sint_derx, sint_dery, sint_derz);

            }
            delete RT1;

        }

    }

    delete [] state_end;
    delete [] state_start;

    GlobalSums(par_gamma_allions, 3*pct.num_owned_ions*max_nl2, pct.kpsub_comm);
    GlobalSums(par_omega_allions, 3*pct.num_owned_ions*max_nl2, pct.kpsub_comm);

    RT1 = new RmgTimer("2-Force: non-local-veff grad");
    OrbitalType *gx = new OrbitalType[FP0_BASIS];
    OrbitalType *gy = new OrbitalType[FP0_BASIS];
    OrbitalType *gz = new OrbitalType[FP0_BASIS];

    CPP_app_grad_driver (&Rmg_L, Rmg_T, veff, (double *)gx, (double *)gy, (double *)gz, FPX0_GRID, FPY0_GRID, FPZ0_GRID, hxxgrid, hyygrid, hzzgrid, ct.kohn_sham_fd_order);
    delete RT1;


    for (ion = 0; ion < pct.num_nonloc_ions; ion++)
    {
        /*Actual index of the ion under consideration*/
        gion = pct.nonloc_ions_list[ion];

        iptr = &ct.ions[gion];

        nh = ct.sp[iptr->species].nh;

        RT1 = new RmgTimer("2-Force: non-local-get gamma");
        GetGamma (gamma, ion, nh, Kptr);
        delete RT1;
        RT1 = new RmgTimer("2-Force: non-local-nlforce_par_Q");
        nlforce_par_Q ((double *)gx, (double *)gy, (double *)gz, gamma, gion, iptr, nh, &qforce[3 * gion]);
        delete RT1;

    }                           /*end for(ion=0; ion<ions_max; ion++) */

    delete [] gz;
    delete [] gy;
    delete [] gx;


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
        //PartialGamma (gion, par_gamma, par_omega, nion, nh, Kptr);
        par_gamma = &par_gamma_allions[ion * 3 * max_nl2];
        par_omega = &par_omega_allions[ion * 3 * max_nl2];
        RT1 = new RmgTimer("2-Force: non-local-nlforce_par_gamma");
        nlforce_par_gamma (par_gamma, gion, nh, &tmp_force_gamma[3*gion]);
        delete RT1;

        RT1 = new RmgTimer("2-Force: non-local-nlforce_par_omega");
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


    //    delete[] par_gamma;
    delete[] par_gamma_allions;
    delete[] par_omega_allions;
#if GPU_ENABLED
    GpuFreeHost(sint_der);
#else
    delete[] sint_der;
#endif

    delete[] gamma;
    delete[] tmp_force_omega;
    delete[] tmp_force_gamma;
    delete[] qforce;


}


