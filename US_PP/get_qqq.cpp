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

void get_qqq ()
{
    int idx, i, j, ion;
    int nh, ncount, icount;
    double *qqq, sum;
    ION *iptr;
    SPECIES *sp;
    FILE *ftpr=NULL;
    char filename[MAX_PATH];
    const bool SET=true;
//if(ct.norm_conserving_pp)return;
    if (pct.gridpe == 0 && verify_boolean ("write_pseudopotential_plots", &SET))
    {
	snprintf (filename, MAX_PATH, "q.txt");
        my_fopen (ftpr, filename, "w+");
    }


    for (ion = 0; ion < ct.num_ions; ion++)
    {
        iptr = &Atoms[ion];
        sp = &Species[iptr->species];

        nh = sp->nh;
        ncount = Atoms[ion].Qindex.size();

        if (iptr->qqq == NULL) {
            iptr->qqq = new double[nh * nh];
            qqq = iptr->qqq;
            for(idx = 0;idx < nh * nh;idx++) qqq[idx] = 0.0;
        }
        qqq = iptr->qqq;

        if (pct.gridpe == 0 && verify_boolean ("write_pseudopotential_plots", &SET))
            fprintf (ftpr, "%% for ion %d :\n", ion);

        idx = 0;
        if(!ct.norm_conserving_pp) {

            for (i = 0; i < nh; i++)
            {
                for (j = i; j < nh; j++)
                {
                    sum = 0.0;
                    if (ncount)
                    {
#pragma omp parallel for reduction(+:sum)
                        for (int icount = 0; icount < ncount; icount++)
                        {
                            sum += GetAugcharge(i, j, icount, ct.cg_coeff.data(), iptr);
                        }
                    }
                    sum = real_sum_all (sum, pct.grid_comm);
                    sum = sum * get_vel_f();
                    if (fabs (sum) < 1.0e-8)
                        sum = 0.0;
                    if (pct.gridpe == 0 && verify_boolean ("write_pseudopotential_plots", &SET))
                        fprintf (ftpr, "i=%d j=%d q_cal=%15.8f q_rel=%15.8f\n",
                                 i, j, sum, sp->qqq[i][j]);

                    qqq[i * nh + j] = sum;
                    if (i != j)
                        qqq[j * nh + i] = qqq[i * nh + j];

                    idx++;
                }                   /*end for j */
            }                       /*end for i */
        } // end if for norm conserving

        if (pct.gridpe == 0 && verify_boolean ("write_pseudopotential_plots", &SET))
            fprintf (ftpr, "\n");
    }                           /*end for ion */
    if (pct.gridpe == 0 && verify_boolean ("write_pseudopotential_plots", &SET))
        fclose (ftpr);

    if(!ct.noncoll) return;

  //implement eq. 18 in Corso and Conte, PRB 71 115106(2005)

    for (ion = 0; ion < ct.num_ions; ion++)
    {
        iptr = &Atoms[ion];
        sp = &Species[iptr->species];

        nh = sp->nh;
        qqq = iptr->qqq;
        if(iptr->qqq_so == NULL) 
            iptr->qqq_so = new std::complex<double>[nh*nh*4];

        for(int idx = 0; idx < nh*nh*4; idx++) iptr->qqq_so[idx] = 0.0;

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
                                        iptr->qqq_so[ih*nh + jh + (is1*2 + is2) * nh * nh] += 
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
                iptr->qqq_so[idx + 0 * nh*nh] = qqq[idx];
                iptr->qqq_so[idx + 3 * nh*nh] = qqq[idx];
            }
        }
    }

}
