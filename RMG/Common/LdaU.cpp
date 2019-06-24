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

template void LdaU<double>::calc_ns_occ(double *);
template void LdaU<std::complex<double>>::calc_ns_occ(std::complex<double> *);
template void LdaU<double>::write_ldaU(void);
template void LdaU<std::complex<double>>::write_ldaU(void);


template <class KpointType> LdaU<KpointType>::LdaU(Kpoint<KpointType> &kp) : K(kp)
{
    this->ldaU_m = 2*ct.max_ldaU_l+1;

    this->Hubbard_J.resize(boost::extents[ct.num_species][3]);
    this->ns_occ.resize(boost::extents[ct.spin_flag+1][ct.num_ions][2*ct.max_ldaU_l+1][2*ct.max_ldaU_l+1]);
}

// Computes the LDA+U occupation matrix
template <class KpointType> void LdaU<KpointType>::calc_ns_occ(KpointType *sint)
{
    int num_tot_proj = K.OrbitalProjector->get_num_tot_proj();
    int num_nonloc_ions = K.OrbitalProjector->get_num_nonloc_ions();
    int *nonloc_ions_list = K.OrbitalProjector->get_nonloc_ions_list();
    int pstride = K.OrbitalProjector->get_pstride();

    size_t alloc = (size_t)num_tot_proj * (size_t)ct.max_states;
    KpointType *sint_compack = new KpointType[alloc]();

    boost::multi_array_ref<KpointType, 3> nsint{sint_compack, boost::extents[K.nstates][ct.num_ions][pstride]};

    // Repack the sint array
    for(int istate = 0; istate < K.nstates; istate++)
    {
        int sindex = istate * num_nonloc_ions * pstride;
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
                double occ = 0.0; 
                for(int st=0;st < K.nstates;st++)
                {
                    occ = occ + K.Kstates[st].occupation[0] * std::real(nsint[st][ion][i] * nsint[st][ion][j]);
                }
                ns_occ[0][ion][i][j] = occ * K.kweight;
            }
        }
    }

    // Get opposite spins
    if(ct.spin_flag)
    {
        MPI_Status status;
        int len = ct.num_ions * pstride * pstride;
        double *sendbuf = ns_occ.data();
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
                        fprintf(ct.logfile, "%7.4f   ", ns_occ[ispin][ion][i][j]);
                    }
                    fprintf(ct.logfile, "\n");
                }
            }
        }
    }
}

// Destructor
template <class KpointType> LdaU<KpointType>::~LdaU(void)
{
}
