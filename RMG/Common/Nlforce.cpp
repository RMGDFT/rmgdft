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
#include "GpuAlloc.h"
#include "blas.h"

//#include "TradeImages.h"
#include "RmgException.h"
#include "Lattice.h"
#include "FiniteDiff.h"
#include "transition.h"
#include "GlobalSums.h"
#include "prototypes_rmg.h"


void nlforce_par_Q (double *gx, double *gy, double *gz, double * gamma, int ion, ION * iptr, int nh, double * forces);


/*Set this to 1 to write out true NL force and the part
 * that comes from eigenvalues*/

template void Nlforce<double> (double *, Kpoint<double> **Kptr, double *force_nl);
template void Nlforce<std::complex<double> > (double * , Kpoint<std::complex<double>> **Kptr, double *force_nl);

template <typename OrbitalType> void Nlforce (double * veff, Kpoint<OrbitalType> **Kptr, double *force_nl)
{
    int ion, isp, gion, nion;
    int nh, size;
    OrbitalType *gamma, *par_gamma, *par_omega;
    ION *iptr;
    int num_ions;
    double  *qforce;
    double *tmp_force_gamma, *tmp_force_omega;

    OrbitalType *psi, *psi_x, *psi_y, *psi_z;

    int num_occupied;
    std::complex<double> I_t(0.0, 1.0);

    int FP0_BASIS = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);

    int num_nonloc_ions = Kptr[0]->BetaProjector->get_num_nonloc_ions();
    int num_owned_ions = Kptr[0]->BetaProjector->get_num_owned_ions();
    int *nonloc_ions_list = Kptr[0]->BetaProjector->get_nonloc_ions_list();
    int *owned_ions_list = Kptr[0]->BetaProjector->get_owned_ions_list();


    RmgTimer *RT1;

    num_ions = Atoms.size();


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
    int max_nl2 = ct.max_nl * ct.max_nl;
    
    gamma = new OrbitalType[ max_nl2];
//    par_gamma = new double[ 6 * max_nl2];
//    par_omega = par_gamma + 3 * max_nl2;
    int factor = ct.noncoll_factor * ct.noncoll_factor;
    OrbitalType *par_gamma_allions = new OrbitalType[3 * num_owned_ions * max_nl2 * factor];
    OrbitalType *par_omega_allions = new OrbitalType[3 * num_owned_ions * max_nl2 * factor];
    for(int i = 0; i < 3 * num_owned_ions * max_nl2 * factor; i++)
    {
        par_gamma_allions[i] = 0.0;
        par_omega_allions[i] = 0.0;
    }

    
    size =  num_nonloc_ions * ct.state_block_size * ct.max_nl * ct.noncoll_factor; 
    size += 1;
    OrbitalType *sint_der = (OrbitalType *)GpuMallocManaged(3*size * sizeof(OrbitalType));
    OrbitalType *sint_derx = sint_der + 0 * size;
    OrbitalType *sint_dery = sint_der + 1 * size;
    OrbitalType *sint_derz = sint_der + 2 * size;


//  determine the number of occupied states for all kpoints.

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
    int pbasis_noncoll = pbasis * ct.noncoll_factor;


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
                psi_x = Kptr[kpt]->Kstates[ct.num_states].psi + (st-state_start[ib]) * pbasis_noncoll;
                psi_y = psi_x + ct.state_block_size*pbasis_noncoll;
                psi_z = psi_x +2* ct.state_block_size*pbasis_noncoll;
                ApplyGradient(psi, psi_x, psi_y, psi_z, ct.force_grad_order, "Coarse");
                if(ct.noncoll)
                    ApplyGradient(psi+pbasis, psi_x+pbasis, psi_y+pbasis, psi_z+pbasis, ct.force_grad_order, "Coarse");

                if(!ct.is_gamma)
                {
                    std::complex<double> *psi_C, *psi_xC, *psi_yC, *psi_zC;
                    psi_C = (std::complex<double> *) psi;
                    psi_xC = (std::complex<double> *) psi_x;
                    psi_yC = (std::complex<double> *) psi_y;
                    psi_zC = (std::complex<double> *) psi_z;
                    for(int i = 0; i < pbasis_noncoll; i++) 
                    {
                        psi_xC[i] += I_t *  Kptr[kpt]->kp.kvec[0] * psi_C[i];
                        psi_yC[i] += I_t *  Kptr[kpt]->kp.kvec[1] * psi_C[i];
                        psi_zC[i] += I_t *  Kptr[kpt]->kp.kvec[2] * psi_C[i];
                    }
                }


            }


            int num_state_thisblock = state_end[ib] - state_start[ib];


            RT1 = new RmgTimer("2-Force: non-local-betaxpsi");
            int st_start = ct.num_states * ct.noncoll_factor;
            int st_thisblock = num_state_thisblock * ct.noncoll_factor;
            int st_block = ct.state_block_size * ct.noncoll_factor;
            Betaxpsi(Kptr[kpt], st_start,              st_thisblock, sint_derx);
            Betaxpsi(Kptr[kpt], st_start +   st_block, st_thisblock, sint_dery);
            Betaxpsi(Kptr[kpt], st_start + 2*st_block, st_thisblock, sint_derz);

            for(int i = 0; i < num_nonloc_ions * st_thisblock * ct.max_nl; i++)
            {
                sint_derx[i] *= -1.0;
                sint_dery[i] *= -1.0;
                sint_derz[i] *= -1.0;
            }
            delete RT1;

            RT1 = new RmgTimer("2-Force: non-local-partial gamma");
            nion = -1;
            for (ion = 0; ion < num_owned_ions; ion++)
            {
                /*Global index of owned ion*/
                gion = owned_ions_list[ion];

                /* Figure out index of owned ion in nonloc_ions_list array, store it in nion*/
                do {

                    nion++;
                    if (nion >= num_nonloc_ions)
                    {
                        printf("\n Could not find matching entry in nonloc_ions_list for owned ion %d", gion);
                        rmg_error_handler(__FILE__, __LINE__, "Could not find matching entry in nonloc_ions_list for owned ion ");
                    }

                } while (nonloc_ions_list[nion] != gion);

                iptr = &Atoms[gion];


                nh = Species[iptr->species].nh;

                /*partial_gamma(ion,par_gamma,par_omega, iptr, nh, p1, p2); */
                par_gamma = &par_gamma_allions[ion * 3 * max_nl2 * factor];
                par_omega = &par_omega_allions[ion * 3 * max_nl2 * factor];

                PartialGamma (kpt, gion, par_gamma, par_omega, nion, nh, Kptr, state_start[ib], state_end[ib], sint_derx, sint_dery, sint_derz);

            }
            delete RT1;

        }

    }

    delete [] state_end;
    delete [] state_start;

    factor *= sizeof(OrbitalType)/sizeof(double);
    GlobalSums(par_gamma_allions, 3*num_owned_ions*max_nl2 * factor, pct.kpsub_comm);
    GlobalSums(par_omega_allions, 3*num_owned_ions*max_nl2 * factor, pct.kpsub_comm);

    RT1 = new RmgTimer("2-Force: non-local-veff grad");
    OrbitalType *gx = new OrbitalType[FP0_BASIS];
    OrbitalType *gy = new OrbitalType[FP0_BASIS];
    OrbitalType *gz = new OrbitalType[FP0_BASIS];

    ApplyGradient (veff, (double *)gx, (double *)gy, (double *)gz, ct.force_grad_order, "Fine");

    delete RT1;


    for (ion = 0; ion < num_nonloc_ions; ion++)
    {
        /*Actual index of the ion under consideration*/
        gion = nonloc_ions_list[ion];

        iptr = &Atoms[gion];

        nh = Species[iptr->species].nh;

        RT1 = new RmgTimer("2-Force: non-local-get gamma");
        GetGamma ((double *)gamma, ion, nh, Kptr);
        delete RT1;
        RT1 = new RmgTimer("2-Force: non-local-nlforce_par_Q");
        nlforce_par_Q ((double *)gx, (double *)gy, (double *)gz, (double *)gamma, gion, iptr, nh, &qforce[3 * gion]);
        delete RT1;

    }                           /*end for(ion=0; ion<ions_max; ion++) */

    delete [] gz;
    delete [] gy;
    delete [] gx;


    double vel_f = get_vel_f();
    for(int i = 0; i < num_ions * 3; i++) qforce[i] *= vel_f;



    /*Loop over ions again */
    nion = -1;
    for (ion = 0; ion < num_owned_ions; ion++)
    {
        /*Global index of owned ion*/
        gion = owned_ions_list[ion];

        /* Figure out index of owned ion in nonloc_ions_list array, store it in nion*/
        do {

            nion++;
            if (nion >= num_nonloc_ions)
            {
                printf("\n Could not find matching entry in nonloc_ions_list for owned ion %d", gion);
                rmg_error_handler(__FILE__, __LINE__, "Could not find matching entry in nonloc_ions_list for owned ion ");
            }

        } while (nonloc_ions_list[nion] != gion);

        iptr = &Atoms[gion];


        nh = Species[iptr->species].nh;

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



    for(int i = 0; i < num_ions * 3; i++) 
    {
        force_nl[i] += qforce[i];
        force_nl[i] += tmp_force_gamma[i];
        force_nl[i] += tmp_force_omega[i];
    }


    if(ct.verbose)
    {
        output_force(qforce, "Non-local forces: QForce");
        output_force(tmp_force_gamma, "Non-local forces: der_gamma term");
        output_force(tmp_force_omega, "Non-local forces: der_omega term");
    }


    //    delete[] par_gamma;
    delete[] par_gamma_allions;
    delete[] par_omega_allions;
    GpuFreeManaged(sint_der);

    delete[] gamma;
    delete[] tmp_force_omega;
    delete[] tmp_force_gamma;
    delete[] qforce;


}


void nlforce_par_Q (double *gx, double *gy, double *gz, double * gamma, int ion, ION * iptr, int nh, double * forces)
{
    int idx2, n, m, count, icount, size;
    int *pidx;

    count = Atoms[ion].Qindex.size();
    pidx = Atoms[ion].Qindex.data();

    double *Qnm;

    Qnm = Atoms[ion].augfunc.data();

    if (count)
    {
        size = (nh * (nh + 1)) / 2;

        idx2 = 0;
        for (n = 0; n < nh; n++)
        {
            for (m = n; m < nh; m++)
            {
                if (m != n) gamma[idx2] *= 2.0;
                idx2++;
            }
        }

        double one = 1.0, zero = 0.0, *tmp_arr;
        int ione = 1;
        tmp_arr = new double[count];

        dgemm("N", "N", &count, &ione, &size, &one, Qnm, &count, gamma, &size, &zero, tmp_arr, &count); 

        for (icount = 0; icount < count; icount++)
        {
            forces[0] -= gx[pidx[icount]] * tmp_arr[icount];
            forces[1] -= gy[pidx[icount]] * tmp_arr[icount];
            forces[2] -= gz[pidx[icount]] * tmp_arr[icount];
        }
        delete [] tmp_arr;
    }


}

template void nlforce_par_gamma<double> (double * par_gamma, int ion, int nh, double *force);
template void nlforce_par_gamma<std::complex<double>> (std::complex<double> * par_gamma, int ion, int nh, double *force);
template <typename T> void nlforce_par_gamma (T * par_gamma, int ion, int nh, double *force)
{
    int idx1, size, n, m;
    if(ct.noncoll) 
    {
        printf("WARNING noncoll forces are not correct");
    }

    T *gamma_x, *gamma_y, *gamma_z;
    double  *dnmI;

    size = nh * nh * ct.noncoll_factor * ct.noncoll_factor;

    gamma_x = par_gamma;
    gamma_y = gamma_x + size;
    gamma_z = gamma_y + size;

    dnmI = Atoms[ion].dnmI;

    for (n = 0; n < nh; n++)
    {
        for (m = 0; m < nh; m++)
        {
            idx1 = n * nh + m;
            force[0] += dnmI[idx1] * std::real(gamma_x[idx1]);
            force[1] += dnmI[idx1] * std::real(gamma_y[idx1]);
            force[2] += dnmI[idx1] * std::real(gamma_z[idx1]);
        }
    }

}


template void nlforce_par_omega<std::complex<double>> (std::complex<double> *par_omega, int ion, int nh, double *force);
template void nlforce_par_omega<double> (double *par_omega, int ion, int nh, double *force);
template <typename OrbitalType> void nlforce_par_omega (OrbitalType *par_omega, int ion, int nh, double *force)
{
    int idx1, size, n, m;

    OrbitalType *omega_x, *omega_y, *omega_z;

    size = nh * nh * ct.noncoll_factor * ct.noncoll_factor;

    omega_x = par_omega;
    omega_y = omega_x + size;
    omega_z = omega_y + size;

    if(ct.noncoll)
    {
        printf("WARNING noncoll forces are not correct");
    }

    double *qqq;
    qqq = Atoms[ion].qqq;

    for (n = 0; n < nh; n++)
    {
        for (m = 0; m < nh; m++)
        {
            idx1 = n * nh + m;
            force[0] -= qqq[idx1] * std::real(omega_x[idx1]);
            force[1] -= qqq[idx1] * std::real(omega_y[idx1]);
            force[2] -= qqq[idx1] * std::real(omega_z[idx1]);

        }
    }

}

