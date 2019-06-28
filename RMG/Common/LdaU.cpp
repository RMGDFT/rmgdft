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
#include "LdaU.h"



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
    this->ldaU_m = 2*ct.max_ldaU_l+1;

    this->Hubbard_J.resize(boost::extents[ct.num_species][3]);
    this->ns_occ.resize(boost::extents[ct.spin_flag+1][ct.num_ions][2*ct.max_ldaU_l+1][2*ct.max_ldaU_l+1]);
}

// Computes the LDA+U occupation matrix. If sint_compack_in is not NULL then it uses that array instead of
// allocating it's own
template <class KpointType> void LdaU<KpointType>::calc_ns_occ(KpointType *sint, int first_state, int num_states)
{
    int num_tot_proj = K.OrbitalProjector->get_num_tot_proj();
    int num_nonloc_ions = K.OrbitalProjector->get_num_nonloc_ions();
    int *nonloc_ions_list = K.OrbitalProjector->get_nonloc_ions_list();
    int pstride = K.OrbitalProjector->get_pstride();

    size_t alloc = (size_t)num_tot_proj * (size_t)ct.max_states;
    KpointType *sint_compack = new KpointType[alloc]();

    boost::multi_array_ref<KpointType, 3> nsint{sint_compack, boost::extents[K.nstates][ct.num_ions][pstride]};

    // Repack the sint array
    for(int istate = 0; istate < num_states; istate++)
    {
        int sindex = (istate + first_state) * num_nonloc_ions * pstride;
        for (int ion = 0; ion < num_nonloc_ions; ion++)
        {
            int proj_index = ion * pstride;
            KpointType *psint = &sint[proj_index + sindex];
            int gion = nonloc_ions_list[ion];
            for (int i = 0; i < pstride; i++)
            {
                //sint_compack[istate * num_tot_proj + proj_index + i] = psint[i];
                nsint[istate][gion][i] = psint[i];
            }
        }
    }



    for(int ion=0;ion < ct.num_ions;ion++)
    {
        for(int i=0;i < this->ldaU_m;i++)
        {
            for(int j=0;j < this->ldaU_m;j++)
            {
                std::complex<double> occ(0.0, 0.0); 
                for(int st=0;st < K.nstates;st++)
                {
                    occ = occ + K.Kstates[st].occupation[0] * nsint[st][ion][i] * nsint[st][ion][j];
                }
                ns_occ[0][ion][i][j] = occ * K.kweight;
            }
        }
    }

    // Impose hermiticity
    for(int ion=0;ion < ct.num_ions;ion++)
    {
        for(int i=0;i < this->ldaU_m;i++)
        {
            for(int j=i+1;j < this->ldaU_m;j++)
            {
                ns_occ[0][ion][i][j] = ns_occ[0][ion][j][i];
            }
        }
    }

    // Need to sum over k-points and symmetrize here then may need to reimpose hermiticity

    // Get occupation matrix from opposite spin case
    if(ct.spin_flag)
    {
        MPI_Status status;
        int len = 2 * ct.num_ions * pstride * pstride;
        double *sendbuf = (double *)ns_occ.data();
        double *recvbuf = sendbuf + len;
        MPI_Sendrecv(sendbuf, len, MPI_DOUBLE, (pct.spinpe+1)%2, pct.gridpe,
            recvbuf, len, MPI_DOUBLE, (pct.spinpe+1)%2, pct.gridpe, pct.spin_comm, &status);

    }
}


// Writes out LDA+U information
template <class KpointType> void LdaU<KpointType>::write_ldaU(void)
{
    if(pct.imgpe!=0) return;

    for(int ion=0;ion < ct.num_ions;ion++)
    {
        ION *iptr = &ct.ions[ion];
        SPECIES *sp = &ct.sp[iptr->species];
        if(sp->num_ldaU_orbitals)
        {
            for(int ispin=0;ispin < ct.spin_flag+1;ispin++)
            {
                fprintf(ct.logfile, "  ion %d spin %d LDA+U occupation matrix\n", ion, ispin);
                for(int i=0;i < ldaU_m;i++)
                {
                    for(int j=0;j < ldaU_m;j++)
                    {
                        fprintf(ct.logfile, "%7.4f   ", std::abs(ns_occ[ispin][ion][i][j]));
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
    KpointType ZERO_t(0.0);
    KpointType ONE_t(1.0);

    int num_tot_proj = K.OrbitalProjector->get_num_tot_proj();
    int num_nonloc_ions = K.OrbitalProjector->get_num_nonloc_ions();
    int pstride = K.OrbitalProjector->get_pstride();

    // allocate memory for sint_compack;
    size_t alloc = (size_t)num_tot_proj * (size_t)ct.max_states;
    KpointType *sint_compack = new KpointType[alloc]();
    KpointType *nwork = new KpointType[alloc];

    // and for the diagonal part of ns_occ
    size_t alloc1 = (size_t)pstride * (size_t)ct.num_ions;
    KpointType *lambda = new KpointType[alloc1](); 
    std::complex<double> *lambda_C = (std::complex<double> *)lambda;
    boost::multi_array_ref<KpointType, 2> nlambda{lambda, boost::extents[ct.num_ions][pstride]};
    boost::multi_array_ref<std::complex<double>, 2> nlambda_C{lambda_C, boost::extents[ct.num_ions][pstride]};

    // Repack the sint array
    boost::multi_array_ref<KpointType, 3> nsint{sint_compack, boost::extents[K.nstates][ct.num_ions][pstride]};
    int *nonloc_ions_list = K.OrbitalProjector->get_nonloc_ions_list();
    for(int istate = 0; istate < num_states; istate++)
    {
        int sindex = (istate + first_state) * num_nonloc_ions * pstride;
        for (int ion = 0; ion < num_nonloc_ions; ion++)
        {
            int proj_index = ion * pstride;
            KpointType *psint = &sint[proj_index + sindex];
            int gion = nonloc_ions_list[ion];
            for (int i = 0; i < pstride; i++)
            {
                //sint_compack[istate * num_tot_proj + proj_index + i] = psint[i];
                nsint[istate][gion][i] = psint[i];
            }
        }
    }


    // Put the diagonal part of ns_occ into a separate array
    this->Ehub = 0.0;
    this->Ecorrect = 0.0;
    for(int ion=0;ion < ct.num_ions;ion++)
    {
        SPECIES *sp = &ct.sp[ct.ions[ion].species];
        double Ueff = sp->Hubbard_U / 2.0;       // FIXME: Have to deal with more complicated cases later
        
        for(int i=0;i < pstride;i++)
        {
            this->Ehub += Ueff * std::real(ns_occ[0][ion][i][i] * (1.0 - ns_occ[0][ion][i][i]));
            this->Ecorrect += Ueff * std::real(ns_occ[0][ion][i][i] * ns_occ[0][ion][i][i]);
            if(ct.is_gamma)
            {
                nlambda[ion][i] = std::real(Ueff * (1.0 - 2.0*ns_occ[0][ion][i][i]));
            }
            else
            {
                nlambda_C[ion][i] = Ueff * (1.0 - 2.0*ns_occ[0][ion][i][i]);
            }
        }
    }

    //MPI_Allreduce(MPI_IN_PLACE, &this->Ehub, 1, MPI_DOUBLE, MPI_SUM, pct.spin_comm);
    //MPI_Allreduce(MPI_IN_PLACE, &this->Ecorrect, 1, MPI_DOUBLE, MPI_SUM, pct.spin_comm);

#if 0
    // Get the eigenvectors of the occupation matrix
    if(ct.ldaU_mode == LDA_PLUS_U_SIMPLE)
    {
        double *work = new double[ct.num_ions];
        for(int ion=0;ion < ct.num_ions;ion++)
        {
       
                     
        }
        delete [] work;
    }

    for(int jj = 0;jj < K.nstates;jj++) {
        for(int ii = 0;ii < num_tot_proj;ii++) {
            nwork[jj*num_tot_proj + ii] = M_dnm[ii] * sint_compack[jj*num_tot_proj + ii];
        }
    }
#endif

    if((ct.scf_steps > 0) || (ct.runflag == RESTART))
    {
        char *transa = "n";
        for(int jj = 0;jj < num_states;jj++) {
            for(int ii = 0;ii < num_tot_proj;ii++) {
                nwork[jj*num_tot_proj + ii] = lambda[ii] * sint_compack[jj*num_tot_proj + ii];
            }
        }

        RmgGemm (transa, transa, K.pbasis, num_states, num_tot_proj,
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
