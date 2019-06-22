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



template LdaU<double>::LdaU(int, int, int);
template LdaU<std::complex<double>>::LdaU(int, int, int);
template LdaU<double>::~LdaU(void);
template LdaU<std::complex<double>>::~LdaU(void);

template void LdaU<double>::calc_ns_occ(void);
template void LdaU<std::complex<double>>::calc_ns_occ(void);



template <class KpointType> LdaU<KpointType>::LdaU(int num_ions, int nspin, int max_ldaU_l) : 
                        ns_occ(boost::extents[num_ions][nspin][2*max_ldaU_l+1][2*max_ldaU_l+1])

{
    this->ldaU_m = 2*max_ldaU_l+1;
}

// Computes the LDA+U occupation matrix
template <class KpointType> void LdaU<KpointType>::calc_ns_occ(void)
{
    int nspin = ct.spin_flag + 1;

    // Zero it out
    for(int ion=0;ion < ct.num_ions;ion++)
    {
        for(int ispin=0;ispin < nspin;ispin++)
        {
            for(int i=0;i < this->ldaU_m;i++)
            {
                for(int j=0;j < this->ldaU_m;i++)
                {
                    this->ns_occ[ion][ispin][i][j] = 0.0;
                }
            }
        }
    }

    // Now generate the occupation matrix

}

// Destructor
template <class KpointType> LdaU<KpointType>::~LdaU(void)
{
}
