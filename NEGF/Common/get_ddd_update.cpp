#include "negf_prototypes.h"
/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include "transition.h"


void get_ddd_update (double * veff)
{
    int idx, i, j, ion;
    int nh, ncount, icount;
    float *qnmI;
    double *dnmI, *sum;
    int *ivec, sum_dim, sum_idx;
    ION *iptr;
    SPECIES *sp;


    if(ct.noncoll) 
    {
        rmg_printf("\n need to change for noncollinear \n");
        fflush(NULL);
        exit(0);
    }

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

    my_calloc (sum, sum_dim, double);


    sum_idx = 0;

    for (ion = 0; ion < ct.num_ions; ion++)
    {
        iptr = &Atoms[ion];
        sp = &Species[iptr->species];

        ivec = Atoms[ion].Qindex.data();
        nh = sp->nh;
        ncount = Atoms[ion].Qindex.size();

        if (iptr->dnmI == NULL)
            iptr->dnmI = new double[nh * nh];

        idx = 0;
        for (i = 0; i < nh; i++)
        {
            for (j = i; j < nh; j++)
            {
                if (ncount)
                {
                    //qnmI = Atoms[ion].augfunc.data() + idx * ncount;
                    for (icount = 0; icount < ncount; icount++)
                    {
                        //sum[sum_idx] += qnmI[icount] * veff[ivec[icount]];
                        double Qr = GetAugcharge(i, j, icount, ct.cg_coeff.data(), iptr);
                        sum[sum_idx] += Qr * veff[ivec[icount]];
                    }
                }               /*end if (ncount) */

                sum[sum_idx] *= get_vel_f();

                idx++;
                sum_idx++;
            }                   /*end for (j = i; j < nh; j++) */
        }                       /*end for (i = 0; i < nh; i++) */

    }                           /*end for (ion = 0; ion < ct.num_ions; ion++) */



    if (sum_idx != sum_dim)
        rmg_error_handler (__FILE__, __LINE__, "Problem with sum index");

    global_sums (sum, &sum_dim, pct.grid_comm);

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

                dnmI[i * nh + j] = sum[sum_idx];
                if (i != j)
                    dnmI[j * nh + i] = dnmI[i * nh + j];

                sum_idx++;

            }                   /*end for j */
        }                       /*end for i */
    }                           /*end for ion */

    my_free (sum);
}
