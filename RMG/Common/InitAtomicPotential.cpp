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

void InitAtomicPotential(double *sumobject)
{


    // this term is included in rhoc and vh term ct.ES and EigSums
    double tpiba = 2.0 * PI / Rmg_L.celldm[0];
    double tpiba2 = tpiba * tpiba;
    double gsquare, k[3];
    SPECIES *sp;
    ION *iptr1;
    double kr, t1, gval, rc, rc2;
    std::complex<double> strfac;
    int P0_BASIS = coarse_pwaves->pbasis;
    int N = ct.prolong_order;
N=12;
    Prolong Pr(2, N, ct.cmix, *Rmg_T,  Rmg_L, *Rmg_G);
    std::vector<std::array<double, 3>> shifts(216);
    double w[216];
    std::fill(w, w+216, 0.0);
    std::fill(sumobject, sumobject+P0_BASIS, 0.0);

    double h = 0.35, norm=0.0;
    int index=0;
    for(int ix=0;ix<N/2;ix++)
    {
        for(int iy=0;iy<N/2;iy++)
        {
            for(int iz=0;iz<N/2;iz++)
            {
                double dist = sqrt(ix*ix + iy*iy + iz*iz);
                int idist = ix*ix + iy*iy + iz*iz;
                if((std::abs(ix) != 0 && std::abs(iy) == 0 && std::abs(iz) == 0) ||
                   (std::abs(ix) == 0 && std::abs(iy) != 0 && std::abs(iz) == 0) ||
                   (std::abs(ix) == 0 && std::abs(iy) == 0 && std::abs(iz) != 0))
                {
                    shifts[index][0] = -h * (double)ix;
                    shifts[index][1] = -h * (double)iy;
                    shifts[index][2] = -h * (double)iz;
                    shifts[index+1][0] = h * (double)ix;
                    shifts[index+1][1] = h * (double)iy;
                    shifts[index+1][2] = h * (double)iz;

                    double factorx = Pr.a[1][ix];
                    double factory = Pr.a[1][iy];
                    double factorz = Pr.a[1][iz];
                    w[index] = factorx*factory*factorz;
                    norm += w[index];

                    factorx = Pr.a[1][N-ix-1];
                    factory = Pr.a[1][N-iy-1];
                    factorz = Pr.a[1][N-iz-1];
                    w[index+1] = factorx*factory*factorz;
                    norm += w[index+1];
                    index+=2;
                }
            }
        }
    }
    w[index] = 1.0;
    norm += 1.0;
    shifts[index][0] = 0.0;
    shifts[index][1] = 0.0;
    shifts[index][2] = 0.0;
    index++;

    std::complex<double> *temp_g = new std::complex<double>[coarse_pwaves->pbasis]();

   // Responsibility for freeing lobject lies in the calling routine
    size_t alloc = (size_t)ct.num_ions * (size_t)P0_BASIS + 128;

    double gcut = sqrt(coarse_pwaves->gcut*tpiba2);

for(int sidx = 0;sidx < index;sidx++)
{
        for(int idx = 0;idx < P0_BASIS;idx++) temp_g[idx] = 0.0;
 
        // for each type of atoms, initilize localpp in g space
        for (int isp = 0; isp < ct.num_species; isp++)
        {
            /* Get species type */
            sp = &Species[isp];

            for(size_t ig=0;ig < coarse_pwaves->pbasis;ig++) 
            {
                if(!coarse_pwaves->gmask[ig]) continue;
                {
                    gsquare = coarse_pwaves->gmags[ig] * tpiba2;
                    gval = sqrt(gsquare);
                    if(gval >= gcut) continue;
                    k[0] = coarse_pwaves->g[ig].a[0] * tpiba;
                    k[1] = coarse_pwaves->g[ig].a[1] * tpiba;
                    k[2] = coarse_pwaves->g[ig].a[2] * tpiba;

                    t1 = AtomicInterpolateInline_Ggrid(sp->localpp_g, gval);

                    // calculating structure factors for each type of atoms
                    strfac = 0.0;
                    for (int i = 0; i < ct.num_ions; i++)
                    {
                        iptr1 = &Atoms[i];
                        if(iptr1->species == isp) 
                        {
                            /*Calculate the phase factor for delocalized case */
                            std::array<double, 3> crds;
                            crds[0] = iptr1->crds[0];
                            crds[1] = iptr1->crds[1];
                            crds[2] = iptr1->crds[2];
                            crds[0] += shifts[sidx][0];
                            crds[1] += shifts[sidx][1];
                            crds[2] += shifts[sidx][2];
                            kr = crds[0] * k[0] + crds[1] * k[1] + crds[2] * k[2];
                            strfac +=  std::exp(std::complex<double>(0.0, -kr));
                        }
                    }
                    temp_g[ig] += strfac * t1 ;

                }
            }
        }


        double omega = Rmg_L.get_omega();
        coarse_pwaves->FftInverse(temp_g, temp_g);
        for(size_t ig = 0; ig < coarse_pwaves->pbasis; ig++) sumobject[ig] += w[sidx]*std::real(temp_g[ig])/omega;
}

    delete []temp_g;

}  

