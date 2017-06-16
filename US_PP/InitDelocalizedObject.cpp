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

void InitDelocalizedObject(double *sumobject, double * &ionobject, int object_type, bool compute_lobject)
{


    // this term is included in rhoc and vh term ct.ES and EigSums
    double tpiba = 2.0 * PI / Rmg_L.celldm[0];
    double tpiba2 = tpiba * tpiba;
    double gsquare, k[3];
    SPECIES *sp;
    ION *iptr1;
    double kr, t1, gval, rc, rc2;
    std::complex<double> strfac;
    size_t ig_storage;
    int FP0_BASIS = fine_pwaves->pbasis;

    ig_storage = fine_pwaves->pbasis;
    if(compute_lobject) ig_storage = ct.num_ions * fine_pwaves->pbasis;

    std::complex<double> *temp_g = new std::complex<double>[ig_storage];

   // Responsibility for freeing lobject lies in the calling routine
    size_t alloc = (size_t)ct.num_ions * (size_t)FP0_BASIS + 128;
    if(compute_lobject) ionobject = new double[alloc]();

    // fully delocalized objects
    pct.num_loc_ions = ct.num_ions;
    for(int ion = 0;ion < ct.num_ions;ion++)  pct.loc_ions_list[ion] = ion;

    for(unsigned int ig=0;ig < ig_storage; ig++) temp_g[ig] = 0.0;

    // for each type of atoms, initilize localpp in g space
    for (int isp = 0; isp < ct.num_species; isp++)
    {
        /* Get species type */
        sp = &ct.sp[isp];
        double Zv = sp->zvalence;

        rc = sp->rc;
        rc2 = rc * rc;
        for(int ig=0;ig < fine_pwaves->pbasis;ig++) 
        {
            if(fine_pwaves->gmags[ig] < ct.filter_factor*ct.filter_factor*fine_pwaves->gcut)
            {
                gsquare = fine_pwaves->gmags[ig] * tpiba2;
                gval = sqrt(gsquare);
                k[0] = fine_pwaves->g[ig].a[0] * tpiba;
                k[1] = fine_pwaves->g[ig].a[1] * tpiba;
                k[2] = fine_pwaves->g[ig].a[2] * tpiba;

                switch(object_type) 
                {
                    case ATOMIC_LOCAL_PP:
                        t1 = AtomicInterpolateInline_Ggrid(sp->localpp_g, gval);
                        break;

                    case ATOMIC_RHO:
                        t1 = AtomicInterpolateInline_Ggrid(sp->arho_g, gval);
                        break;

                    case ATOMIC_RHOCOMP:
                        t1= Zv * exp (-rc2 * gsquare/4.0);
                        break;

                    case ATOMIC_RHOCORE: 
                        if(sp->nlccflag) 
                        {
                            t1 = AtomicInterpolateInline_Ggrid(sp->rhocore_g, gval);
                        }
                        else
                        {
                            t1 = 0.0;
                        }
                        break;

                    default:
                        throw RmgFatalException() << "Undefined local object type" << 
                            " in " << __FILE__ << " at line " << __LINE__ << "\n";

                }

                if(compute_lobject) 
                {
                    for (int ion = 0; ion < ct.num_ions; ion++)
                    {

                        iptr1 = &ct.ions[ion];
                        if(iptr1->species == isp) 
                        {
                            kr = iptr1->crds[0] * k[0] + iptr1->crds[1] * k[1] + iptr1->crds[2] * k[2];
                            strfac =  std::exp(std::complex<double>(0.0, -kr));
                            temp_g[(size_t)ion * (size_t)FP0_BASIS + ig] += strfac * t1;
                        }
                    }
                }
                else
                { 
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


                    temp_g[ig] += strfac * t1 ;
                    //if(gsquare > 1.0e-6) temp_g[ig] += strfac * Zj * std::exp(-rc2 * gsquare/4.0)/gsquare;
                }
            }
        }
    }


    if(compute_lobject) 
    {
        for (int ion = 0; ion < ct.num_ions; ion++)
        {
            std::complex<double> *ptr_g = &temp_g[(size_t)ion * (size_t)FP0_BASIS];
            PfftInverse(ptr_g, ptr_g, *fine_pwaves);
            for(int ig = 0; ig < fine_pwaves->pbasis; ig++) 
                ionobject[(size_t)ion * (size_t)FP0_BASIS + ig] = std::real(ptr_g[ig])/Rmg_L.get_omega();
        }
    }
    else
    {
        PfftInverse(temp_g, temp_g, *fine_pwaves);
        for(int ig = 0; ig < fine_pwaves->pbasis; ig++) sumobject[ig] = std::real(temp_g[ig])/Rmg_L.get_omega();
    }


    delete []temp_g;


    if(compute_lobject) return;

    // Renormalize atomic rho
    if(object_type == ATOMIC_RHO) {
        double t2 = 0.0;
        for (int idx = 0; idx < FP0_BASIS; idx++) t2 += sumobject[idx];
        t2 = get_vel_f() *  real_sum_all (t2, pct.grid_comm);
        double t1 = ct.nel / t2;
        double difference = fabs(t1 - 1.0);
        if ((ct.verbose == 1) || (difference > 0.05))
        {
            if (pct.imgpe == 0)
                printf ("\n LCAO initialization: Normalization constant for initial atomic charge is %f\n", t1);
        }

        for(int idx = 0;idx < FP0_BASIS;idx++) sumobject[idx] *= t1;

    }

    // Core charges may have a small negative component because of filtering
    if(object_type == ATOMIC_RHOCORE) 
    {
        for (int idx = 0; idx < FP0_BASIS; idx++) 
        {
            if(sumobject[idx] < 0.0) sumobject[idx] = 0.0;
        }
    }

    if(object_type == ATOMIC_RHOCOMP)
    {

        /* Check compensating charges */
        ct.crho = 0.0;
        for (int idx = 0; idx < FP0_BASIS; idx++) ct.crho += sumobject[idx];

        ct.crho = ct.crho * get_vel_f();
        ct.crho = real_sum_all (ct.crho, pct.grid_comm);  /* sum over pct.grid_comm  */

        if (ct.verbose)
        {
            if (pct.imgpe==0)
                printf("\nCompensating charge is %.8e\n", ct.crho);
        }
    }

    if(ct.runflag == 5 || ct.runflag == 6 || ct.forceflag == TDDFT) return;

    if(object_type == ATOMIC_LOCAL_PP) init_efield (sumobject);

}  

