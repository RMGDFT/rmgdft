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
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "AtomicInterpolate.h"
#include "transition.h"
#include "Atomic.h"
#include "RmgParallelFft.h"
#include "RmgException.h"

// This is used to initialize vnuc with option ct.localize_localpp = "false"
// on the high density grid. Each object can have a sum representation
// where the sum is over ions and a local representation which separates
// the contributions of individual ions and is used to calculate forces. 
// The object types are.
//

void InitDelocalizedLocalpp(double *vlocpp_r)
{


    // this term is included in rhoc and vh term ct.ES and EigSums
    double tpiba = 2.0 * PI / Rmg_L.celldm[0];
    double tpiba2 = tpiba * tpiba;
    double gsquare, k[3];
    SPECIES *sp;
    ION *iptr1;
    double kr, t1, gval, rc, rc2;
    std::complex<double> strfac;

    std::complex<double> *vlocpp_g = new std::complex<double>[fine_pwaves->pbasis];


    for(int ig=0;ig < fine_pwaves->pbasis;ig++) vlocpp_g[ig] = 0.0;

    // for each type of atoms, initilize localpp in g space
    for (int isp = 0; isp < ct.num_species; isp++)
    {
        /* Get species type */
        sp = &ct.sp[isp];
        double Zj = sp->zvalence;

        rc = sp->rc;
        rc2 = rc * rc;

        for(int ig=0;ig < fine_pwaves->pbasis;ig++) 
        {
            if(fine_pwaves->gmags[ig] < fine_pwaves->gcut)
            {
                gsquare = fine_pwaves->gmags[ig] * tpiba2;
                gval = sqrt(gsquare);
                k[0] = fine_pwaves->g[ig].a[0] * tpiba;
                k[1] = fine_pwaves->g[ig].a[1] * tpiba;
                k[2] = fine_pwaves->g[ig].a[2] * tpiba;

                // calculating structure factors for each type of atoms
                strfac = 0.0;
                for (int i = 0; i < ct.num_ions; i++)
                {

                    iptr1 = &ct.ions[i];
                    if(iptr1->species == isp) 
                    {
                        kr = iptr1->crds[0] * k[0] + iptr1->crds[1] * k[1] + iptr1->crds[2] * k[2];
                        strfac +=  std::exp(std::complex<double>(0.0, -kr));
                    }
                }


                t1 = AtomicInterpolateInline_Ggrid(sp->localpp_g, gval);
                vlocpp_g[ig] += strfac * t1 ;
                //if(gsquare > 1.0e-6) vlocpp_g[ig] += strfac * Zj * std::exp(-rc2 * gsquare/4.0)/gsquare;
            }
        }
    }

    PfftInverse(vlocpp_g, vlocpp_g, *fine_pwaves);
    for(int ig = 0; ig < fine_pwaves->pbasis; ig++) vlocpp_r[ig] = std::real(vlocpp_g[ig])/Rmg_L.get_omega();


    delete []vlocpp_g;

}

