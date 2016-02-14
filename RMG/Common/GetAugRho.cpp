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
#include <string.h>
#include "transition.h"
#include "const.h"
#include "RmgTimer.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "Kpoint.h"
#include <complex>
#include "../Headers/prototypes.h"


void GetAugRho(double *augrho)
{

    int max_product = (ct.max_nl + 1) * ct.max_nl / 2;
    double *qtpr;

    int FP0_BASIS = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);

    double *product = new double[max_product];

    // Zero out augrho
    for(int idx = 0;idx < FP0_BASIS;idx++) augrho[idx] = 0.0;


    if(!ct.norm_conserving_pp) {

        for (int ion = 0; ion < pct.num_nonloc_ions; ion++)
        {
            int gion = pct.nonloc_ions_list[ion];
            
            if (pct.Qidxptrlen[gion])
            {
                
                ION *iptr = &ct.ions[gion];
           
                int nh = ct.sp[iptr->species].nh;
                
                int *ivec = pct.Qindex[gion];
                int ncount = pct.Qidxptrlen[gion];
                double *qnmI = pct.augfunc[gion];

                for (int i=0; i < max_product; i++)
                    product[i] = 0.0;


                int idx = 0;
                for (int i = 0; i < nh; i++)
                {
                    for (int j = i; j < nh; j++)
                    {
                        qtpr = qnmI + idx * ncount;
                        for (int icount = 0; icount < ncount; icount++)
                        {
                            augrho[ivec[icount]] += qtpr[icount] * product[idx];
                        }           /*end for icount */
                        idx++;
                    }               /*end for j */
                }                   /*end for i */


            }                       /*end if */

        }                           /*end for ion */

    }

    // Do I need to do this?
    symmetrize_rho (augrho);

    delete [] product;
}
