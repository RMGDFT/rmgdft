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

#include <math.h>
#include <float.h>
#include <complex>
#include <iostream>
#include <fstream>
#include <fcntl.h>
#include <sys/stat.h>

#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "RmgException.h"
#include "RmgSumAll.h"
#include "transition.h"


#include "RmgParallelFft.h"

void VhcorrPeriodicPart(double *vh_x, double *vh_y, double *vh_z, double alpha, double *r0)
{

    double gsquare, gx, g_r0;
    std::complex<double> phase_r0;
    ptrdiff_t densgrid[3];
    pfft_plan plan_forward, plan_back;
    densgrid[0] = fine_pwaves->global_dimx;
    densgrid[1] = fine_pwaves->global_dimy;
    densgrid[2] = fine_pwaves->global_dimz;

    double g2cut = (sqrt(fine_pwaves->gmax))*(sqrt(fine_pwaves->gmax));
    int global_basis = fine_pwaves->global_basis;
    int pbasis = fine_pwaves->pbasis;

    std::complex<double> *crho = new std::complex<double>[pbasis];
    plan_forward = pfft_plan_dft_3d(densgrid,
            (double (*)[2])crho,
            (double (*)[2])crho,
            pct.pfft_comm,
            PFFT_FORWARD,
            PFFT_TRANSPOSED_NONE|PFFT_ESTIMATE);

    plan_back = pfft_plan_dft_3d(densgrid,
            (pfft_complex *)crho,
            (pfft_complex *)crho,
            pct.pfft_comm,
            PFFT_BACKWARD,
            PFFT_TRANSPOSED_NONE|PFFT_ESTIMATE);


    //  for(int ig=0;ig < pbasis;ig++) {
    //      if(pwaves.gmags[ig] > g2cut) {
    //          crho[ig] = std::complex<double>(0.0, 0.0);
    //      }
    //  }


    double tpiba = 2.0 * PI / Rmg_L.celldm[0];
    double tpiba2 = tpiba * tpiba;

    for(int ig=0;ig < pbasis;ig++) {
        if(fine_pwaves->gmags[ig] > 1.0e-6) 
        {   
            gsquare = fine_pwaves->gmags[ig] * tpiba2;
            g_r0  = fine_pwaves->g[ig].a[0] *r0[0];
            g_r0 += fine_pwaves->g[ig].a[1] *r0[1];
            g_r0 += fine_pwaves->g[ig].a[2] *r0[2];
            g_r0 *= tpiba;
            phase_r0 = exp(std::complex<double>(0.0, -g_r0));

            gx  = fine_pwaves->g[ig].a[0] * tpiba;
            crho[ig] = exp(-alpha * alpha * gsquare/4.0)/gsquare * gx *phase_r0;
        }
        else
        {
            crho[ig] = 0.0;
        }
        
    }

    pfft_execute_dft(plan_back, (double (*)[2])crho, (double (*)[2])crho);
    for(int i = 0;i < pbasis;i++) vh_x[i] = std::imag(crho[i]) * 4 * PI/Rmg_L.get_omega();


    for(int ig=0;ig < pbasis;ig++) {
        if(fine_pwaves->gmags[ig] > 1.0e-6) 
        {   
            gsquare = fine_pwaves->gmags[ig] * tpiba2;
            g_r0  = fine_pwaves->g[ig].a[0] *r0[0];
            g_r0 += fine_pwaves->g[ig].a[1] *r0[1];
            g_r0 += fine_pwaves->g[ig].a[2] *r0[2];
            g_r0 *= tpiba;
            phase_r0 = exp(std::complex<double>(0.0, -g_r0));

            gx  = fine_pwaves->g[ig].a[1] * tpiba;
            crho[ig] = exp(-alpha * alpha * gsquare/4.0)/gsquare * gx *phase_r0;
        }
        else
        {
            crho[ig] = 0.0;
        }
        
    }

    pfft_execute_dft(plan_back, (double (*)[2])crho, (double (*)[2])crho);
    for(int i = 0;i < pbasis;i++) vh_y[i] = std::imag(crho[i]) * 4 * PI/Rmg_L.get_omega();

    for(int ig=0;ig < pbasis;ig++) {
        if(fine_pwaves->gmags[ig] > 1.0e-6) 
        {   
            gsquare = fine_pwaves->gmags[ig] * tpiba2;
            g_r0  = fine_pwaves->g[ig].a[0] *r0[0];
            g_r0 += fine_pwaves->g[ig].a[1] *r0[1];
            g_r0 += fine_pwaves->g[ig].a[2] *r0[2];
            g_r0 *= tpiba;
            phase_r0 = exp(std::complex<double>(0.0, -g_r0));

            //printf("\n  iggg %d %e %e %e %e %e", ig, fine_pwaves->g[ig].a[2],g_r0, g_r0/PI, phase_r0)

            gx  = fine_pwaves->g[ig].a[2] * tpiba;
            crho[ig] = exp(-alpha * alpha * gsquare/4.0)/gsquare * gx *phase_r0;
        }
        else
        {
            crho[ig] = 0.0;
        }
        
    }

    pfft_execute_dft(plan_back, (double (*)[2])crho, (double (*)[2])crho);
    for(int i = 0;i < pbasis;i++) vh_z[i] = std::imag(crho[i]) * 4 * PI/Rmg_L.get_omega();




    pfft_destroy_plan(plan_back);
    pfft_destroy_plan(plan_forward);
    delete [] crho;

}

