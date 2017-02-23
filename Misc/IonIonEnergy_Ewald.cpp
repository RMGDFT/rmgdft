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


#include <complex>
#include "const.h"
#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "State.h"
#include "Kpoint.h"
#include "Functional.h"
#include "GlobalSums.h"
#include "transition.h"
#include "RmgSumAll.h"
#include <boost/math/special_functions/erf.hpp>


/* Evaluate total ion-ion energy by ewald method 
   written by Wenchang Lu NCSU*/
static std::complex<double> structure_factor(double *k);

double IonIonEnergy_Ewald ()
{
    double r, x, y, z;
    ION *iptr1, *iptr2;
    int i, j;
    
    double total_ii;
    double t1;

    double sigma = 0.0;
    if(!ct.localize_localpp) sigma = 3.0; 
    for (i = 0; i<ct.num_species; i++) 
        sigma = std::max(sigma, ct.sp[i].rc);


    //  the ion-ion term in real space
    //  check how many unit cell to make erfc(nL/sqrt(2.0)*sigma) < tol ==1.0e-8
    //  erfc(4.06) = 9.37e-9 so we set the largest r/(sqrt(2)*sigma) = 4.06


    double rcutoff = 4.06;
    int num_cell_x = int(rcutoff * sqrt(2.0) * sigma/Rmg_L.get_xside()) + 1;
    int num_cell_y = int(rcutoff * sqrt(2.0) * sigma/Rmg_L.get_yside()) + 1;
    int num_cell_z = int(rcutoff * sqrt(2.0) * sigma/Rmg_L.get_zside()) + 1;


    double ii_real_space = 0.0;

    for (i = pct.gridpe; i < ct.num_ions; i+=pct.grid_npes)
    {

        iptr1 = &ct.ions[i];
        double Zi = ct.sp[iptr1->species].zvalence;

        for (j = 0; j < ct.num_ions; j++)
        {

            iptr2 = &ct.ions[j];
            double Zj = ct.sp[iptr2->species].zvalence;
            t1 = sqrt (sigma);
            if(ct.localize_localpp)
                t1 = sqrt (ct.sp[iptr1->species].rc * ct.sp[iptr1->species].rc +
                        ct.sp[iptr2->species].rc * ct.sp[iptr2->species].rc);
            
            for(int ix = -num_cell_x; ix<= num_cell_x; ix++)
            for(int iy = -num_cell_y; iy<= num_cell_y; iy++)
            for(int iz = -num_cell_z; iz<= num_cell_z; iz++)
            {
                x = iptr1->crds[0] - iptr2->crds[0] + ix * Rmg_L.a0[0] + iy * Rmg_L.a1[0] + iz * Rmg_L.a2[0];
                y = iptr1->crds[1] - iptr2->crds[1] + ix * Rmg_L.a0[1] + iy * Rmg_L.a1[1] + iz * Rmg_L.a2[1];
                z = iptr1->crds[2] - iptr2->crds[2] + ix * Rmg_L.a0[2] + iy * Rmg_L.a1[2] + iz * Rmg_L.a2[2];
                r = sqrt(x*x + y*y + z*z);

                // r= 0 means two atoms are the same one.
                if(ct.localize_localpp)
                {
                    if(r > 1.0e-5) ii_real_space += Zi * Zj/r * boost::math::erfc(r/t1);
                }
                else
                {
                    if(r > 1.0e-5) ii_real_space += Zi * Zj * boost::math::erfc(t1*r) / r;
                }

            }
        }

    }

//   reciprocal space term
    double ii_kspace = 0.0;

    // this term is included in rhoc and vh term ct.ES and EigSums
    // so it is not necessary to include it when using localized localpp
    // but when using delocalized rhoc does not exist so we need it
    if(!ct.localize_localpp)
    {
        if(pct.gridpe == 0) ii_kspace = -ct.nel*ct.nel / sigma / 4.0;
        double tpiba = 2.0 * PI / Rmg_L.celldm[0];
        double tpiba2 = tpiba * tpiba;
        double gsquare, k[3];
        std::complex<double> S;

        for(int ig=0;ig < fine_pwaves->pbasis;ig++)
        {
            if(fine_pwaves->gmags[ig] > 1.0e-6)
            {
                gsquare = fine_pwaves->gmags[ig] * tpiba2;
                k[0] = fine_pwaves->g[ig].a[0] * tpiba;
                k[1] = fine_pwaves->g[ig].a[1] * tpiba;
                k[2] = fine_pwaves->g[ig].a[2] * tpiba;

                S = structure_factor(k);

                ii_kspace += std::norm(S) * exp(-gsquare/sigma/4.0)/ fine_pwaves->gmags[ig] / tpiba2;
            }
         }
         ii_kspace = 4.0*PI/Rmg_L.omega * ii_kspace;
    }


    // term self 
    double ii_self = 0.0;

    if(ct.localize_localpp)
    {
        for (i = pct.gridpe; i < ct.num_ions; i+=pct.grid_npes)
            ii_self -= ct.sp[ct.ions[i].species].zvalence *
                ct.sp[ct.ions[i].species].zvalence / (sqrt (2.0 * PI) * ct.sp[ct.ions[i].species].rc);
    }
    else
    {
        for (i = pct.gridpe; i < ct.num_ions; i+=pct.grid_npes)
            ii_self -= ct.sp[ct.ions[i].species].zvalence *
                ct.sp[ct.ions[i].species].zvalence * sqrt (4.0 * sigma / PI);
        ii_self *= 0.5;
    }

    total_ii =  0.5 * ii_real_space + 0.5 * ii_kspace + ii_self;
    MPI_Allreduce(MPI_IN_PLACE, &total_ii, 1, MPI_DOUBLE, MPI_SUM, pct.grid_comm);

    return total_ii;

}

static std::complex<double> structure_factor(double *k)
{
    ION *iptr1;
    double kr;
    std::complex<double> S = 0.0;

    for (int i = 0; i < ct.num_ions; i++)
    {

        iptr1 = &ct.ions[i];
        double Zi = ct.sp[iptr1->species].zvalence;
        kr = iptr1->crds[0] * k[0] + iptr1->crds[1] * k[1] + iptr1->crds[2] * k[2];
        S +=  Zi * std::exp(std::complex<double>(0.0, kr));
    }

    return S;
}




