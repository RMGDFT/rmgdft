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



template LdaU<double>::LdaU(Kpoint<double> *, int, int);
template LdaU<std::complex<double>>::LdaU(Kpoint<std::complex<double>> *, int, int);
template LdaU<double>::~LdaU(void);
template LdaU<std::complex<double>>::~LdaU(void);

template void LdaU<double>::calc_ns_occ(double *);
template void LdaU<std::complex<double>>::calc_ns_occ(std::complex<double> *);



template <class KpointType> LdaU<KpointType>::LdaU(Kpoint<KpointType> *kptr, int num_ions, int max_ldaU_l) : 
                        ns_occ(boost::extents[num_ions][2*max_ldaU_l+1][2*max_ldaU_l+1])

{
    this->ldaU_m = 2*max_ldaU_l+1;
    this->K = kptr;
}

// Computes the LDA+U occupation matrix
template <class KpointType> void LdaU<KpointType>::calc_ns_occ(KpointType *sint)
{
    int nstates = this->K->nstates;
    int num_tot_proj = this->K->OrbitalProjector->get_num_tot_proj();
    int num_nonloc_ions = this->K->OrbitalProjector->get_num_nonloc_ions();
    int *nonloc_ions_list = this->K->OrbitalProjector->get_nonloc_ions_list();
    int pstride = this->K->OrbitalProjector->get_pstride();

    size_t alloc = (size_t)num_tot_proj * (size_t)ct.max_states;
    KpointType *sint_compack = new KpointType[alloc]();

    boost::multi_array_ref<KpointType, 3> nsint{sint_compack, boost::extents[nstates][ct.num_ions][pstride]};

    // Repack the sint array
    for(int istate = 0; istate < nstates; istate++)
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
                for(int st=0;st < nstates;st++)
                {
                    occ = occ + this->K->Kstates[st].occupation[0] * std::real(nsint[st][ion][i] * nsint[st][ion][j]);
                }
                this->ns_occ[ion][i][j] = occ * this->K->kweight;
if(pct.gridpe==0)
{
printf("PPPP  %d  %d  %d  %14.8e\n",ion,i,j,this->ns_occ[ion][i][j]);
}
            }
            }
    }

    delete [] sint_compack;
}

// Destructor
template <class KpointType> LdaU<KpointType>::~LdaU(void)
{
}
