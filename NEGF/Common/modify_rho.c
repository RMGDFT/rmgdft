/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include "twoParts.h"


void modify_rho (double * rho, double * rho_old)
{
    int idx, ione = 1;
    double t2;
    register double tcharge;
    int i, j, k;
    double total_charge, tcharge_fixed, t_fixed;
    int xoff, yoff, zoff;
    int test;
    int x1, x2, y1, y2, z1, z2;
    int idx1, idx2, item;
    double *array_global;
    int fpbasis;
    
    fpbasis = get_FP0_BASIS();

    if (ct.runflag == 111)
    {

        idx1 = get_FNX_GRID() * get_FPY0_GRID() * get_FPZ0_GRID();
        my_malloc_init( array_global, idx1, double );

        item = get_FNX_GRID() / 3;
        if (item * 3 - get_FNX_GRID() != 0)
            error_handler ("run flag=111");
        distribute_to_X_soft (rho, array_global);

        for (i = 0; i < get_FNX_GRID(); i++)
            for (j = 0; j < get_FPY0_GRID() * get_FPZ0_GRID(); j++)
            {
                idx = i * get_FPY0_GRID() * get_FPZ0_GRID() + j;
                idx1 = (i + item) * get_FPY0_GRID() * get_FPZ0_GRID() + j;
                idx2 = (i - item) * get_FPY0_GRID() * get_FPZ0_GRID() + j;

                if (i <= item)
                    array_global[idx] = array_global[idx1];
                if (i >= 2 * item)
                    array_global[idx] = array_global[idx2];
            }

        X_to_distribute_soft (array_global, rho);

        idx1 = get_FNX_GRID() * get_FPY0_GRID() * get_FPZ0_GRID();
        my_free(array_global);

    }

    /*  Fixed  charge density in some region -----Qingzhong     05/17/03/  */
    if (chargeDensityCompass.type == 2)
    {

        /*  normalize the charge density  */

        tcharge = 0.0;
        for (idx = 0; idx < get_FP0_BASIS(); idx++)
            tcharge += rho[idx];
        ct.tcharge = real_sum_all (tcharge, pct.grid_comm);
        for (idx = 0; idx < get_FP0_BASIS(); idx++)
            rho[idx] = rho_old[idx];


    }

    if (chargeDensityCompass.type == 0)
    {

        /*  normalize the charge density  */

        tcharge = 0.0;
        for (idx = 0; idx < get_FP0_BASIS(); idx++)
            tcharge += rho[idx];
        ct.tcharge = real_sum_all (tcharge, pct.grid_comm) * get_vel_f();
        if (pct.gridpe == 0)
            printf ("total charge %10.4f = %10.4f + %10.4f\n",
                    ct.tcharge, ct.tcharge - ct.nel, ct.nel);

        t2 = ct.nel / ct.tcharge;

        dscal (&fpbasis, &t2, &rho[0], &ione);
    }

    if (chargeDensityCompass.type == 1)
    {

        x1 = chargeDensityCompass.box1.x1;
        x2 = chargeDensityCompass.box1.x2;
        y1 = chargeDensityCompass.box1.y1;
        y2 = chargeDensityCompass.box1.y2;
        z1 = chargeDensityCompass.box1.z1;
        z2 = chargeDensityCompass.box1.z2;

        xoff = get_FPX_OFFSET();
        yoff = get_FPY_OFFSET();
        zoff = get_FPZ_OFFSET();

        total_charge = 0.0;
        tcharge_fixed = 0.0;
        for (i = 0; i < get_FPX0_GRID(); i++)
        {
            for (j = 0; j < get_FPY0_GRID(); j++)
            {
                for (k = 0; k < get_FPZ0_GRID(); k++)
                {
                    idx = i * get_FPY0_GRID() * get_FPZ0_GRID() + j * get_FPZ0_GRID() + k;

                    test = (((i + xoff) < x1) || ((i + xoff) >= x2) ||
                            ((j + yoff) < y1) || ((j + yoff) >= y2) ||
                            ((k + zoff) < z1) || ((k + zoff) >= z2));
                    if (!test)
                    {
                        total_charge += rho[i * get_FPY0_GRID() * get_FPZ0_GRID() + j * get_FPZ0_GRID() + k];

                    }
                    else
                    {
                        rho[idx] = rho_old[idx];
                        tcharge_fixed += rho_old[i * get_FPY0_GRID() * get_FPZ0_GRID() + j * get_FPZ0_GRID() + k];
                    }
                }
            }
        }

        t2 = real_sum_all (total_charge, pct.grid_comm) * get_vel_f();
        t_fixed = real_sum_all (tcharge_fixed, pct.grid_comm) * get_vel_f();

        if (pct.gridpe == 0)
            printf ("total charge %10.4f + %10.4f = %10.4f = %10.4f + %10.4f\n",
                    t2, t_fixed, t2 + t_fixed, t2 + t_fixed - ct.nel, ct.nel);

        /*
           t2 = (ct.nel - t_fixed) / (t2 * get_vel_f());
           t2 = 1.0 / get_vel_f();

           for (i = 0; i < get_FPX0_GRID(); i++)
           {
           for (j = 0; j < get_FPY0_GRID(); j++)
           {
           for (k = 0; k < get_FPZ0_GRID(); k++)
           {
           idx = i * get_FPY0_GRID() * get_FPZ0_GRID() + j * get_FPZ0_GRID() + k;
           test = (((i + xoff) < chargeDensityCompass.box1.x1)
           || ((i + xoff) >= chargeDensityCompass.box1.x2)
           || ((j + yoff) < chargeDensityCompass.box1.y1)
           || ((j + yoff) >= chargeDensityCompass.box1.y2)
           || ((k + zoff) < chargeDensityCompass.box1.z1)
           || ((k + zoff) >= chargeDensityCompass.box1.z2));
           if (!test)
           {
           rho[i * get_FPY0_GRID() * get_FPZ0_GRID() + j * get_FPZ0_GRID() + k] *= t2;

           }
           }
           }
           }
         */

    }
    MPI_Barrier(pct.img_comm);


}
