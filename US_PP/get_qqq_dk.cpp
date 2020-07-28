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
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "main.h"
#include "transition.h"
#include "GlobalSums.h"

void get_qqq_dk (double dk_xtal[3], std::complex<double> *qqq_dk, std::complex<double> *qqq_dk_so)
{
    int idx, i, j, ion;
    int nh, ncount, icount;
    std::complex<double> *qqq, sum;
    ION *iptr;
    SPECIES *sp;

    for(size_t idx = 0; idx < Atoms.size() * ct.max_nl * ct.max_nl; idx++) qqq_dk[idx] = 0.0;
    int px_offset = Rmg_G->get_PX_OFFSET(ct.FG_RATIO);
    int py_offset = Rmg_G->get_PY_OFFSET(ct.FG_RATIO);
    int pz_offset = Rmg_G->get_PZ_OFFSET(ct.FG_RATIO);
    int px0_grid = Rmg_G->get_PX0_GRID(ct.FG_RATIO);
    int py0_grid = Rmg_G->get_PY0_GRID(ct.FG_RATIO);
    int pz0_grid = Rmg_G->get_PZ0_GRID(ct.FG_RATIO);
    double hxgrid = Rmg_G->get_hxgrid(ct.FG_RATIO);
    double hygrid = Rmg_G->get_hygrid(ct.FG_RATIO);
    double hzgrid = Rmg_G->get_hzgrid(ct.FG_RATIO);
    int nx_grid = Rmg_G->get_NX_GRID(ct.FG_RATIO);
    int ny_grid = Rmg_G->get_NY_GRID(ct.FG_RATIO);
    int nz_grid = Rmg_G->get_NZ_GRID(ct.FG_RATIO);
    std::complex<double> *phase_dk = new std::complex<double>[px0_grid * py0_grid * pz0_grid];
    for(int ix = 0; ix < px0_grid; ix++)
    {
        for(int iy = 0; iy < py0_grid; iy++)
        {
            for(int iz = 0; iz < pz0_grid; iz++)
            {
                int ixx = ix + px_offset;
       //      if(ixx > nx_grid/2) ixx -= nx_grid;
               int iyy = iy + py_offset;
        //       if(iyy > ny_grid/2) iyy -= ny_grid;
                int izz = iz + pz_offset;
         //       if(izz > nz_grid/2) izz -= nz_grid;
                double kr = hxgrid * ixx * dk_xtal[0] + hygrid * iyy * dk_xtal[1] + hzgrid * izz * dk_xtal[2];
                phase_dk[ix * py0_grid * pz0_grid + iy * pz0_grid + iz] = std::exp(std::complex<double>(0.0,   kr * twoPI));

            //if(ix == 0) printf("\n ccc %d %d  %f %f %f", iy, iz, kr, phase_dk[iy*pz0_grid+iz]);
            }
        }
    }

    for (ion = 0; ion < ct.num_ions; ion++)
    {
        iptr = &Atoms[ion];
        sp = &Species[iptr->species];

        nh = sp->nh;
        qqq = &qqq_dk[ion * ct.max_nl * ct.max_nl];
        ncount = Atoms[ion].Qindex.size();
        double ktao = dk_xtal[0] * Atoms[ion].xtal[0] +  dk_xtal[1] * Atoms[ion].xtal[1] + dk_xtal[2] * Atoms[ion].xtal[2];
        ktao = 0.0;
        std::complex<double> phase_ion = std::exp(std::complex<double>(0.0,  -ktao * twoPI));

        idx = 0;
        if(!ct.norm_conserving_pp) {

            for (i = 0; i < nh; i++)
            {
                for (j = i; j < nh; j++)
                {
                    sum = 0.0;
                    if (ncount)
                    {
                        for (icount = 0; icount < ncount; icount++)
                        {
                            sum += (double)Atoms[ion].augfunc[icount + idx * ncount] * phase_dk[Atoms[ion].Qindex[icount]] * phase_ion;
                            int idx = Atoms[ion].Qindex[icount];
                            int ix = idx/pz0_grid/py0_grid;
                            int iy = (idx/pz0_grid)%py0_grid;
                            int iz = idx%pz0_grid;
                         // if(i == 0 && j == 0 && ion == 0 && std::abs(Atoms[ion].augfunc[icount]) > 0.01)
                         //       printf("\n aaa %d %d %d %f  %f %f", ix, iy, iz, Atoms[ion].augfunc[icount], phase_dk[idx] * phase_ion);
                            
                        }
                    }
                    GlobalSums (&sum, 1, pct.grid_comm);
                    sum = sum * get_vel_f();
                    if (pct.gridpe == 0  && 0)
                        printf ( "\n i=%d j=%d q_cal=%15.8f  %15.8f q_rel=%15.8f",
                                i, j, std::real(sum), std::imag(sum), sp->qqq[i][j]);

                    qqq[i * nh + j] = sum;
                    if (i != j)
                        qqq[j * nh + i] = qqq[i * nh + j];

                    idx++;
                }                   /*end for j */
            }                       /*end for i */
        } // end if for norm conserving

    }                           /*end for ion */

    if(!ct.noncoll) return;

    //implement eq. 18 in Corso and Conte, PRB 71 115106(2005)
    for(size_t idx = 0; idx < 4* Atoms.size() * ct.max_nl * ct.max_nl; idx++) qqq_dk_so[idx] = 0.0;
    std::complex<double> *qqq_so;

    for (ion = 0; ion < ct.num_ions; ion++)
    {
        iptr = &Atoms[ion];
        sp = &Species[iptr->species];

        nh = sp->nh;
        qqq = &qqq_dk[ion * ct.max_nl * ct.max_nl];
        qqq_so = &qqq_dk_so[ion * ct.max_nl * ct.max_nl * 4];

        if(sp->is_spinorb)
        {
            for(int ih =0; ih < nh; ih++)
                for(int jh =0; jh < nh; jh++)
                {
                    for(int is1 = 0; is1 < 2; is1++)
                        for(int is2 = 0; is2 < 2; is2++)
                        {
                            for(int m1 = 0; m1 < nh; m1++)
                                for(int m2 = 0; m2 < nh; m2++)
                                    for(int is0 = 0; is0 < 2; is0++)
                                    {
                                        qqq_so[ih*nh + jh + (is1*2 + is2) * nh * nh] += 
                                            qqq[m1*nh+m2] * sp->fcoef_so[ih][m1][is1*2+is0] *
                                            sp->fcoef_so[m2][jh][is0*2+is2];  
                                    }
                        }
                }
        }
        else
        {
            for(int idx = 0; idx < nh*nh; idx++)
            {
                qqq_so[idx + 0 * nh*nh] = qqq[idx];
                qqq_so[idx + 3 * nh*nh] = qqq[idx];
            }
        }
    }

}
