#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "main.h"
#include "transition.h"

static std::complex<double> DnmTransform(int ih, int jh, int is1, int is2, double *Ia, SPECIES &sp);
void get_ddd (double * veff, double *vxc)
{
    int idx, ion;
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
        for (int i = 0; i < nh; i++)
        {
            for (int j = i; j < nh; j++)
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

        for (int i = 0; i < nh; i++)
        {
            for (int j = i; j < nh; j++)
            {

                /*if (fabs (sum[sum_idx]) < 1.0e-10)
                  sum[sum_idx] = 0.0; */

                dnmI[i * nh + j] = sp->ddd0[i][j] + sum[sum_idx];
                if (i != j)
                    dnmI[j * nh + i] = dnmI[i * nh + j];

                sum_idx++;

            }                   /*end for j */
        }                       /*end for i */
    }                           /*end for ion */

// implement eq 19 Corso and Conte PRB 75 115106(2005)

    if(!ct.noncoll) 
    {
        delete [] sum;
        return;
    }

    sum_idx = 0;
    double *Ia = new double[ct.max_nl * ct.max_nl * 4];
    for (ion = 0; ion < ct.num_ions; ion++)
    {
        iptr = &Atoms[ion];
        sp = &Species[iptr->species];

        nh = sp->nh;

        if(iptr->dnmI_so == NULL)
            iptr->dnmI_so = new std::complex<double>[nh*nh*4];

        for (int ih = 0; ih < nh; ih++)
            for (int jh = ih; jh < nh; jh++)
            {
                for(int is = 0; is < 4; is++)
                {
                    Ia[ih * nh + jh + is *nh*nh] = sum[is * sum_dim + sum_idx];
                    Ia[jh * nh + ih + is *nh*nh] = sum[is * sum_dim + sum_idx];
                }
                sum_idx++;
            }


        for (int ih = 0; ih < nh; ih++)
        {
            for (int jh = 0; jh < nh; jh++)
            {
                for(int is1 = 0; is1 <2; is1++)
                    for(int is2 = 0; is2 <2; is2++)
                    {
                        iptr->dnmI_so[ih *nh + jh + (is1*2+is2) * nh * nh] = DnmTransform(ih,jh,is1,is2,Ia,*sp);
                    }

                iptr->dnmI_so[ih *nh + jh + 0 * nh * nh] += sp->ddd0_so[ih][jh][0];
                iptr->dnmI_so[ih *nh + jh + 1 * nh * nh] += sp->ddd0_so[ih][jh][1];
                iptr->dnmI_so[ih *nh + jh + 2 * nh * nh] += sp->ddd0_so[ih][jh][2];
                iptr->dnmI_so[ih *nh + jh + 3 * nh * nh] += sp->ddd0_so[ih][jh][3];

            }                   /*end for j */
        }                       /*end for i */
    }                           /*end for ion */

    delete [] Ia;
    delete [] sum;
}

static std::complex<double> DnmTransform(int ih, int jh, int is1, int is2, double *Ia, SPECIES &sp)
{
    std::complex<double> value;
    value = 0.0;
    int nh = sp.nh;
    if(sp.is_spinorb)
    {
        for(int m1 = 0; m1 < nh; m1++)
            for(int m2 = 0; m2 < nh; m2++)
            {
                value += Ia[m1*nh+m2 + 0*nh*nh] *
                    (sp.fcoef_so[ih][m1][is1 *2 +0] * sp.fcoef_so[m2][jh][0*2 + is2] +
                     sp.fcoef_so[ih][m1][is1 *2 +1] * sp.fcoef_so[m2][jh][1*2 + is2]);
                value += Ia[m1*nh+m2 + 1*nh*nh] *
                    (sp.fcoef_so[ih][m1][is1 *2 +1] * sp.fcoef_so[m2][jh][0*2 + is2] +
                     sp.fcoef_so[ih][m1][is1 *2 +0] * sp.fcoef_so[m2][jh][1*2 + is2]);
                value += Ia[m1*nh+m2 + 2*nh*nh] * std::complex<double>(0.0, -1.0) * 
                    (sp.fcoef_so[ih][m1][is1 *2 +1] * sp.fcoef_so[m2][jh][0*2 + is2] -
                     sp.fcoef_so[ih][m1][is1 *2 +0] * sp.fcoef_so[m2][jh][1*2 + is2]);
                value += Ia[m1*nh+m2 + 3*nh*nh] *
                    (sp.fcoef_so[ih][m1][is1 *2 +0] * sp.fcoef_so[m2][jh][0*2 + is2] -
                     sp.fcoef_so[ih][m1][is1 *2 +1] * sp.fcoef_so[m2][jh][1*2 + is2]);
            }
    }
    else
    {
        if(is1 == 0 && is2 == 0)
            value = Ia[ih * nh + jh + 0*nh*nh] + Ia[ih *nh + jh + 3*nh*nh];
        if(is1 == 0 && is2 == 1)
            value = std::complex<double>(Ia[ih * nh + jh + 1*nh*nh], -Ia[ih *nh + jh + 2*nh*nh]);
        if(is1 == 1 && is2 == 0)
            value = std::complex<double>(Ia[ih * nh + jh + 1*nh*nh],  Ia[ih *nh + jh + 2*nh*nh]);
        if(is1 == 1 && is2 == 1)
            value = Ia[ih * nh + jh + 0*nh*nh] - Ia[ih *nh + jh + 3*nh*nh];
    }

    return value;
}
