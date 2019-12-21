#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "main.h"
#include "transition.h"

void get_ddd (double * veff, double *vxc)
{
    int idx, i, j, ion;
    int nh, ncount, icount;
    double *dnmI, *sum;
    int *ivec, sum_dim, sum_idx;
    ION *iptr;
    SPECIES *sp;

    int FP0_BASIS =  Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);

    /*Count the number of elements in sum array */
    sum_dim = 0;
    for (ion = 0; ion < ct.num_ions; ion++)
    {
        iptr = &Atoms[ion];
        sp = &Species[iptr->species];
        
        nh = sp->nh;

        /*Number of elements is sum of 1+2+3+...+nh */
        sum_dim += nh * (nh + 1) / 2;
    }

    sum = new double[sum_dim * ct.noncoll_factor * ct.noncoll_factor]();


    sum_idx = 0;

    for (ion = 0; ion < ct.num_ions; ion++)
    {
        iptr = &Atoms[ion];
        sp = &Species[iptr->species];

        ivec = Atoms[ion].Qindex.data();
        nh = sp->nh;
        ncount = Atoms[ion].Qindex.size();

        if (iptr->dnmI == NULL)
            iptr->dnmI = new double[nh * nh * ct.noncoll_factor * ct.noncoll_factor];

        idx = 0;
        for (i = 0; i < nh; i++)
        {
            for (j = i; j < nh; j++)
            {
                if (ncount)
                {
                    for (icount = 0; icount < ncount; icount++)
                    {
                        sum[sum_idx] += Atoms[ion].augfunc[icount + idx * ncount] * veff[ivec[icount]];
                        if(ct.noncoll)
                        {
                            sum[1*sum_dim + sum_idx] += Atoms[ion].augfunc[icount + idx * ncount] 
                                * vxc[1*FP0_BASIS + ivec[icount]];
                            sum[2*sum_dim + sum_idx] += Atoms[ion].augfunc[icount + idx * ncount] 
                                * vxc[2*FP0_BASIS + ivec[icount]];
                            sum[3*sum_dim + sum_idx] += Atoms[ion].augfunc[icount + idx * ncount] 
                                * vxc[3*FP0_BASIS + ivec[icount]];
                        }
                    }
                }               /*end if (ncount) */

                idx++;
                sum_idx++;
            }                   /*end for (j = i; j < nh; j++) */
        }                       /*end for (i = 0; i < nh; i++) */

    }                           /*end for (ion = 0; ion < ct.num_ions; ion++) */



    if (sum_idx != sum_dim)
        error_handler ("Problem with sum index");

    int num_sums = sum_dim * ct.noncoll_factor * ct.noncoll_factor;
    global_sums (sum, &num_sums, pct.grid_comm);

    for(int idx = 0; idx < num_sums; idx++) sum[idx] *= get_vel_f();

    sum_idx = 0;

    for (ion = 0; ion < ct.num_ions; ion++)
    {
        iptr = &Atoms[ion];
        sp = &Species[iptr->species];

        nh = sp->nh;

        dnmI = iptr->dnmI;

        for (i = 0; i < nh; i++)
        {
            for (j = i; j < nh; j++)
            {

                /*if (fabs (sum[sum_idx]) < 1.0e-10)
                  sum[sum_idx] = 0.0; */

                dnmI[i * nh + j] = sp->ddd0[i][j] + sum[sum_idx];
                if (i != j)
                    dnmI[j * nh + i] = dnmI[i * nh + j];

                if(ct.noncoll)
                {
                    dnmI[i * nh + j+ 1*nh*nh] = sum[1*sum_dim + sum_idx];
                    if (i != j)
                        dnmI[j * nh + i + 1*nh*nh] = dnmI[i * nh + j + 1* nh*nh];

                    dnmI[i * nh + j+ 2*nh*nh] = sum[2*sum_dim + sum_idx];
                    if (i != j)
                        dnmI[j * nh + i + 2*nh*nh] = dnmI[i * nh + j + 2* nh*nh];

                    dnmI[i * nh + j+ 3*nh*nh] = sum[3*sum_dim + sum_idx];
                    if (i != j)
                        dnmI[j * nh + i + 3*nh*nh] = dnmI[i * nh + j + 3* nh*nh];
                }

                sum_idx++;

            }                   /*end for j */
        }                       /*end for i */
    }                           /*end for ion */

    delete [] sum;
}
