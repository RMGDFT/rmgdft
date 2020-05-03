/*
 *
 * Copyright 2019 The RMG Project Developers. See the COPYRIGHT file 
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
#include <complex>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <libgen.h>

#include "const.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "Kpoint.h"
#include "RmgSumAll.h"
#include "RmgGemm.h"
#include "GpuAlloc.h"
#include "FiniteDiff.h"
#include "LdaU.h"
#include "transition.h"



template LdaU<double>::LdaU(Kpoint<double> &);
template LdaU<std::complex<double>>::LdaU(Kpoint<std::complex<double>> &);
template LdaU<double>::~LdaU(void);
template LdaU<std::complex<double>>::~LdaU(void);

template void LdaU<double>::calc_ns_occ(double *, int, int);
template void LdaU<std::complex<double>>::calc_ns_occ(std::complex<double> *, int, int);
template void LdaU<double>::app_vhubbard(double *, double *, int, int);
template void LdaU<std::complex<double>>::app_vhubbard(std::complex<double> *, std::complex<double> *, int, int);
template void LdaU<double>::write_ldaU(void);
template void LdaU<std::complex<double>>::write_ldaU(void);


template <class KpointType> LdaU<KpointType>::LdaU(Kpoint<KpointType> &kp) : K(kp)
{
    this->ldaU_m = (2*ct.max_ldaU_l+1);

    this->Hubbard_J.resize(boost::extents[ct.num_species][3]);
    this->ns_occ.resize(boost::extents[ct.nspin][Atoms.size()][this->ldaU_m][this->ldaU_m]);
}

// Computes the LDA+U occupation matrix. If sint_compack_in is not NULL then it uses that array instead of
// allocating it's own
template <class KpointType> void LdaU<KpointType>::calc_ns_occ(KpointType *sint, int first_state, int num_states)
{
    int num_tot_proj = K.OrbitalProjector->get_num_tot_proj();
    int num_nonloc_ions = K.OrbitalProjector->get_num_nonloc_ions();
    int *nonloc_ions_list = K.OrbitalProjector->get_nonloc_ions_list();
    int pstride = K.OrbitalProjector->get_pstride();

    size_t alloc = (size_t)num_tot_proj * (size_t)ct.max_states * ct.noncoll_factor;
    KpointType *sint_compack = new KpointType[alloc]();

    if(first_state != 0) 
    {
        printf("\n first_state in calc_ns_occ must be 0 but it is  %d", first_state);
        rmg_error_handler(__FILE__, __LINE__, "wrong first_state");
    }
    boost::multi_array_ref<KpointType, 4> nsint{sint_compack, boost::extents[K.nstates][ct.noncoll_factor][Atoms.size()][pstride]};

    // Repack the sint array
    for(int istate = 0; istate < num_states; istate++)
    {
        size_t sindex = (istate + first_state) * num_nonloc_ions * pstride * ct.noncoll_factor;
        for(int ispin = 0; ispin < ct.noncoll_factor; ispin++)
        { 
            for (int ion = 0; ion < num_nonloc_ions; ion++)
            {
                int proj_index = ion * pstride;
                KpointType *psint = &sint[proj_index + sindex + ispin * num_nonloc_ions* pstride];
                int gion = nonloc_ions_list[ion];
                for (int i = 0; i < pstride; i++)
                {
                    //sint_compack[istate * num_tot_proj + proj_index + i] = psint[i];
                    nsint[istate][ispin][gion][i] = psint[i];

                }
            }
        }
    }


    for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
    {
        for(int is1 = 0; is1 < ct.noncoll_factor; is1++)
        {
            for(int is2 = 0; is2 < ct.noncoll_factor; is2++)
            {

                for(int i=0;i < this->ldaU_m;i++)
                {
                    for(int j=0;j < this->ldaU_m;j++)
                    {
                        std::complex<double> occ(0.0, 0.0); 
                        for(int st=0;st < K.nstates;st++)
                        {
                            occ = occ + K.Kstates[st].occupation[0] * nsint[st][is1][ion][i] * std::conj(nsint[st][is2][ion][j]);
                        }
                        ns_occ[is1 * ct.noncoll_factor + is2][ion][i][j] = occ * K.kp.kweight;
                    }
                }
            }
        }
    }

    // Need to sum over k-points and symmetrize here then may need to reimpose hermiticity

    // Get occupation matrix from opposite spin case
    if(ct.nspin == 2)
    {
        MPI_Status status;
        int len = 2 * Atoms.size() * pstride * pstride;
        double *sendbuf = (double *)ns_occ.data();
        double *recvbuf = sendbuf + len;
        MPI_Sendrecv(sendbuf, len, MPI_DOUBLE, (pct.spinpe+1)%2, pct.gridpe,
                recvbuf, len, MPI_DOUBLE, (pct.spinpe+1)%2, pct.gridpe, pct.spin_comm, &status);

    }
    delete [] sint_compack;
}


// Writes out LDA+U information
template <class KpointType> void LdaU<KpointType>::write_ldaU(void)
{
    if(pct.imgpe!=0) return;

    for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
    {
        if(Species[Atoms[ion].species].num_ldaU_orbitals)
        {
            fprintf(ct.logfile, "  ion %lu  LDA+U occupation matrix_real\n", ion);
            for(int is1 = 0; is1 < ct.noncoll_factor; is1++)
            {
                for(int i=0;i < ldaU_m;i++)
                {
                    for(int is2 = 0; is2 < ct.noncoll_factor; is2++)
                    {
                        for(int j=0;j < ldaU_m;j++)
                        {
                            fprintf(ct.logfile, "%7.4f ", std::real(ns_occ[is1*ct.noncoll_factor + is2][ion][i][j]));
                        }
                    }
                    fprintf(ct.logfile, "\n");
                }
            }
        }
    }
}


// Applies the hubbard potential to orbitals V_hub|psi>
template <class KpointType> void LdaU<KpointType>::app_vhubbard(KpointType *v_hub_x_psi, KpointType *sint, int first_state, int num_states)
{
    KpointType ONE_t(1.0), zero_t(0.0);

    int num_tot_proj = K.OrbitalProjector->get_num_tot_proj();
    int num_nonloc_ions = K.OrbitalProjector->get_num_nonloc_ions();
    int pstride = K.OrbitalProjector->get_pstride();

    // allocate memory for sint_compack;
    size_t alloc = (size_t)num_tot_proj * (size_t)ct.max_states * ct.noncoll_factor;
    KpointType *sint_compack = new KpointType[alloc]();
    KpointType *nwork = new KpointType[alloc];

    // and for the diagonal part of ns_occ
    size_t alloc1 = (size_t)pstride * (size_t)Atoms.size() * ct.noncoll_factor;
    KpointType *lambda = new KpointType[alloc1 * alloc1](); 
    std::complex<double> *lambda_C = (std::complex<double> *)lambda;
    boost::multi_array_ref<KpointType, 6> nlambda{lambda,
        boost::extents[ct.noncoll_factor][Atoms.size()][pstride][ct.noncoll_factor][Atoms.size()][pstride]};
    boost::multi_array_ref<std::complex<double>, 6> nlambda_C{lambda_C,
        boost::extents[ct.noncoll_factor][Atoms.size()][pstride][ct.noncoll_factor][Atoms.size()][pstride]};

    // Repack the sint array
    boost::multi_array_ref<KpointType, 4> nsint{sint_compack, boost::extents[K.nstates][ct.noncoll_factor][Atoms.size()][pstride]};
    int *nonloc_ions_list = K.OrbitalProjector->get_nonloc_ions_list();

    for(int istate = 0; istate < num_states; istate++)
    {
        size_t sindex = (istate + first_state) * num_nonloc_ions * pstride * ct.noncoll_factor;
        for(int ispin = 0; ispin < ct.noncoll_factor; ispin++)
        { 
            for (int ion = 0; ion < num_nonloc_ions; ion++)
            {
                int proj_index = ion * pstride;
                KpointType *psint = &sint[proj_index + sindex + ispin * num_nonloc_ions* pstride];
                int gion = nonloc_ions_list[ion];
                for (int i = 0; i < pstride; i++)
                {
                    //sint_compack[istate * num_tot_proj + proj_index + i] = psint[i];
                    nsint[istate][ispin][gion][i] = psint[i];

                }
            }
        }
    }


    // Put the diagonal part of ns_occ into a separate array
    this->Ehub = 0.0;
    this->Ecorrect = 0.0;

    for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
    {
        SPECIES &AtomType = Species[Atoms[ion].species];
        double Ueff = AtomType.Hubbard_U / 2.0;       // FIXME: Have to deal with more complicated cases later

        //Tr(ns_occ * ns_occ)
        for(int is1 = 0; is1 < ct.noncoll_factor; is1++)
        {
            for(int is2 = 0; is2 < ct.noncoll_factor; is2++)
            {

                int ispin = is1 * ct.noncoll_factor + is2;
                for(int i=0;i < pstride;i++)
                    for(int j=0;j < pstride;j++)
                    {
                        this->Ehub -= Ueff * std::real( ns_occ[ispin][ion][i][j] * ns_occ[ispin][ion][j][i]);
                    }

                if(is1 == is2)
                {
                    for(int i=0;i < pstride;i++)
                    {
                        this->Ehub += Ueff * std::real(ns_occ[ispin][ion][i][i]);
                        this->Ecorrect += Ueff * std::real(ns_occ[ispin][ion][i][i] * ns_occ[ispin][ion][i][i]);
                    }
                }
            }
        }
    }

    //MPI_Allreduce(MPI_IN_PLACE, &this->Ehub, 1, MPI_DOUBLE, MPI_SUM, pct.spin_comm);
    //MPI_Allreduce(MPI_IN_PLACE, &this->Ecorrect, 1, MPI_DOUBLE, MPI_SUM, pct.spin_comm);

    if((ct.scf_steps > 0) || (ct.runflag == RESTART))
    {
        for(int is1 = 0; is1 < ct.noncoll_factor; is1++)
            for(int is2 = 0; is2 < ct.noncoll_factor; is2++)
            {
                int ispin = is1 * ct.noncoll_factor + is2;
                for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
                {
                    SPECIES &AtomType = Species[Atoms[ion].species];
                    double Ueff = AtomType.Hubbard_U / 2.0;       // FIXME: Have to deal with more complicated cases later

                    for(int i=0;i < pstride;i++)
                    {
                        for(int j=0;j < pstride;j++)
                        {
                            if(ct.is_gamma)
                            {
                                if(i==j && is1 == is2)
                                    nlambda[is1][ion][i][is2][ion][j] = std::real(Ueff * (1.0 - 2.0*ns_occ[ispin][ion][i][j]));
                                else
                                    nlambda[is1][ion][i][is2][ion][j] = -std::real(Ueff * 2.0*ns_occ[ispin][ion][i][j]);

                            }
                            else
                            {
                                if(i==j && is1 == is2)
                                    nlambda_C[is1][ion][i][is2][ion][j] = Ueff * (1.0 - 2.0*ns_occ[ispin][ion][i][j]);
                                else
                                    nlambda_C[is1][ion][i][is2][ion][j] = -Ueff * 2.0*ns_occ[ispin][ion][i][j];
                            }
                        }
                    }
                }
            }


        char *transa = "n";
        int num_tot_proj_nc = num_tot_proj * ct.noncoll_factor;
        RmgGemm (transa, transa, num_tot_proj_nc, num_states, num_tot_proj_nc,
                ONE_t, lambda, num_tot_proj_nc, sint_compack, num_tot_proj_nc,
                zero_t, nwork, num_tot_proj_nc);

        int num_states_nc = num_states * ct.noncoll_factor;
        RmgGemm (transa, transa, K.pbasis, num_states_nc, num_tot_proj,
                ONE_t, K.orbital_weight, K.pbasis, nwork, num_tot_proj,
                ONE_t, v_hub_x_psi, K.pbasis);

    }

    delete [] lambda;
    delete [] nwork;
    delete [] sint_compack;
}


// Destructor
template <class KpointType> LdaU<KpointType>::~LdaU(void)
{
}

// calculate LDA+u force for this kpoint
template void LdaU<double>::calc_force(double *sint, double *force_ldau);
template void LdaU<std::complex<double>>::calc_force(std::complex<double> *sint, double *force_ldau);
template <class KpointType> void LdaU<KpointType>::calc_force(KpointType *sint, double *force_ldau)
{
    KpointType *psi, *psi_x, *psi_y, *psi_z;

    std::complex<double> I_t(0.0, 1.0);

    int num_tot_proj = K.OrbitalProjector->get_num_tot_proj();
    int num_nonloc_ions = K.OrbitalProjector->get_num_nonloc_ions();
    int num_owned_ions = K.OrbitalProjector->get_num_owned_ions();

    int *nonloc_ions_list = K.OrbitalProjector->get_nonloc_ions_list();
    int *owned_ions_list = K.OrbitalProjector->get_owned_ions_list();

    int pstride = K.OrbitalProjector->get_pstride();

    size_t size =  (size_t)num_tot_proj * ct.state_block_size * sizeof(KpointType);
    KpointType *sint_der = (KpointType *)GpuMallocManaged(3*size * sizeof(KpointType));
    KpointType *sint_derx = sint_der + 0 * size;
    KpointType *sint_dery = sint_der + 1 * size;
    KpointType *sint_derz = sint_der + 2 * size;

    for(int idx = 0; idx < (int)Atoms.size() * 3; idx++) force_ldau[idx] = 0.0;
    //  determine the number of occupied states for all kpoints.

    int num_occupied = 0;
    for(int st = 0; st < ct.num_states; st++)
    {
        if(abs(K.Kstates[st].occupation[0]) < 1.0e-10) 
        {
            num_occupied = std::max(num_occupied, st);
            break;
        }
    }


    int num_state_block = (num_occupied +ct.state_block_size -1) / ct.state_block_size;
    int *state_start = new int[num_state_block];
    int *state_end = new int[num_state_block];
    std::complex<double> *par_occ_x = new std::complex<double>[Atoms.size() * pstride * pstride]();
    std::complex<double> *par_occ_y = new std::complex<double>[Atoms.size() * pstride * pstride]();
    std::complex<double> *par_occ_z = new std::complex<double>[Atoms.size() * pstride * pstride]();

    for(int ib = 0; ib < num_state_block; ib++)
    {
        state_start[ib] = ib * ct.state_block_size;
        state_end[ib] = (ib+1) * ct.state_block_size;
        if(state_end[ib] > num_occupied) state_end[ib] = num_occupied;
    }

    int pbasis = K.pbasis;
    int pbasis_noncoll = pbasis * ct.noncoll_factor;


    if(ct.alloc_states < ct.num_states + 3 * ct.state_block_size)     
    {
        printf("\n allco_states %d should be larger than ct.num_states %d + 3* ct.state_block_size, %d", ct.alloc_states, ct.num_states,
                ct.state_block_size);
        exit(0);
    }

    for(int ib = 0; ib < num_state_block; ib++)
    {
        for(int st = state_start[ib]; st < state_end[ib]; st++)
        {
            psi = K.Kstates[st].psi;
            psi_x = K.Kstates[ct.num_states].psi + (st-state_start[ib]) * pbasis_noncoll;
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
                    psi_xC[i] += I_t *  K.kp.kvec[0] * psi_C[i];
                    psi_yC[i] += I_t *  K.kp.kvec[1] * psi_C[i];
                    psi_zC[i] += I_t *  K.kp.kvec[2] * psi_C[i];
                }
            }


        }


        int num_state_thisblock = state_end[ib] - state_start[ib];


        RmgTimer *RT1 = new RmgTimer("2-Force: LDAU: non-local-betaxpsi");
        int st_start = ct.num_states * ct.noncoll_factor;
        int st_thisblock = num_state_thisblock * ct.noncoll_factor;
        int st_block = ct.state_block_size * ct.noncoll_factor;

        K.OrbitalProjector->project(&K, sint_derx, st_start,              st_thisblock, K.orbital_weight);
        K.OrbitalProjector->project(&K, sint_dery, st_start +   st_block, st_thisblock, K.orbital_weight);
        K.OrbitalProjector->project(&K, sint_derz, st_start + 2*st_block, st_thisblock, K.orbital_weight);

        //        for(int i = 0; i < num_nonloc_ions * st_thisblock * pstride; i++)
        //        {
        //            sint_derx[i] *= -1.0;
        //            sint_dery[i] *= -1.0;
        //            sint_derz[i] *= -1.0;
        //        }
        delete RT1;

        RT1 = new RmgTimer("2-Force: LDAU: non-local-partial");
        int nion = -1;
        for (int ion = 0; ion < num_owned_ions; ion++)
        {
            /*Global index of owned ion*/
            int gion = owned_ions_list[ion];

            /* Figure out index of owned ion in nonloc_ions_list array, store it in nion*/
            do {

                nion++;
                if (nion >= num_nonloc_ions)
                {
                    printf("\n Could not find matching entry in nonloc_ions_list for owned ion %d", gion);
                    rmg_error_handler(__FILE__, __LINE__, "Could not find matching entry in nonloc_ions_list for owned ion ");
                }

            } while (nonloc_ions_list[nion] != gion);

            for(int st = state_start[ib]; st < state_end[ib]; st++)
            {
                double occ_st = K.Kstates[st].occupation[0];
                int sindex = st * num_nonloc_ions * pstride + nion * pstride;
                int der_index = (st - state_start[ib]) * num_nonloc_ions * pstride + nion * pstride;

                for (int m1 = 0; m1 < pstride; m1++)
                {
                    for (int m2 = 0; m2 < pstride; m2++)
                    {
                        par_occ_x[ion * pstride * pstride + m1 * pstride + m2] += 
                            occ_st*(sint_derx[der_index + m1] * std::conj(sint[sindex + m2]) + 
                                    sint[sindex + m1] * MyConj(sint_derx[der_index + m2])); 

                        par_occ_y[ion * pstride * pstride + m1 * pstride + m2] += 
                            occ_st*(sint_dery[der_index + m1] * std::conj(sint[sindex + m2]) + 
                                    sint[sindex + m1] * MyConj(sint_dery[der_index + m2])); 
                        par_occ_z[ion * pstride * pstride + m1 * pstride + m2] += 
                            occ_st*(sint_derz[der_index + m1] * std::conj(sint[sindex + m2]) + 
                                    sint[sindex + m1] * MyConj(sint_derz[der_index + m2])); 
                    }
                }

            }
        }
        delete RT1;

    }

    std::complex<double> sum_x, sum_y, sum_z;
    for (int ion = 0; ion < num_owned_ions; ion++)
    {
        int gion = owned_ions_list[ion];
        sum_x = 0.0;
        sum_y = 0.0;
        sum_z = 0.0;
        for(int m1 = 0; m1 < pstride; m1++)
        {
            sum_x += par_occ_x[ion *pstride *pstride + m1 * pstride + m1];
            sum_y += par_occ_y[ion *pstride *pstride + m1 * pstride + m1];
            sum_z += par_occ_z[ion *pstride *pstride + m1 * pstride + m1];
            for(int m2 = 0; m2 <pstride; m2++)
            {
                sum_x -= 2.0 * ns_occ[0][gion][m2][m1] * par_occ_x[ion *pstride *pstride + m1 * pstride + m2];
                sum_y -= 2.0 * ns_occ[0][gion][m2][m1] * par_occ_y[ion *pstride *pstride + m1 * pstride + m2];
                sum_z -= 2.0 * ns_occ[0][gion][m2][m1] * par_occ_z[ion *pstride *pstride + m1 * pstride + m2];
            }

        }

        SPECIES &AtomType = Species[Atoms[gion].species];
        double Ueff = AtomType.Hubbard_U / 2.0; 

        force_ldau[gion * 3 + 0] = -Ueff * std::real(sum_x);
        force_ldau[gion * 3 + 1] = -Ueff * std::real(sum_y);
        force_ldau[gion * 3 + 2] = -Ueff * std::real(sum_z);

        if(ct.verbose && pct.gridpe == 0)
        {
            printf("  ion %d LDA+U occupation matrix\n", ion);
            for(int i=0;i < ldaU_m;i++)
            {
                for(int j=0;j < ldaU_m;j++)
                {
                    printf("%7.4f   ", std::abs(ns_occ[0][ion][i][j]));
                }
                printf("\n");
            }
            printf("  ion %d LDA+U occupation_x matrix\n", ion);
            for(int i=0;i < ldaU_m;i++)
            {
                for(int j=0;j < ldaU_m;j++)
                {
                    printf("%7.4f   ", std::abs(par_occ_x[ion *pstride *pstride + i * pstride + j]));
                }
                printf("\n");
            }
        }

    }


    delete [] state_end;
    delete [] state_start;
    delete [] par_occ_x;
    delete [] par_occ_y;
    delete [] par_occ_z;
    GpuFreeManaged(sint_der);

}
