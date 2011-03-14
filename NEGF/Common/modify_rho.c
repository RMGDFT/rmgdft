/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "md.h"


void modify_rho (REAL * rho, REAL * rho_old)
{
    int idx, ione = 1;
    REAL t2;
    register double tcharge;
    int i, j, k;
    REAL total_charge, tcharge_fixed, t_fixed;
    int pex, pey, pez, xoff, yoff, zoff;
    int test;
    int x1, x2, y1, y2, z1, z2;
    int idx1, idx2, item;
    REAL *array_global;

    if (ct.runflag == 111)
    {

        idx1 = FNX_GRID * FPY0_GRID * FPZ0_GRID;
        my_malloc_init( array_global, idx1, REAL );

        item = FNX_GRID / 3;
        if (item * 3 - FNX_GRID != 0)
            error_handler ("run flag=111");
        distribute_to_X_soft (rho, array_global);

        for (i = 0; i < FNX_GRID; i++)
            for (j = 0; j < FPY0_GRID * FPZ0_GRID; j++)
            {
                idx = i * FPY0_GRID * FPZ0_GRID + j;
                idx1 = (i + item) * FPY0_GRID * FPZ0_GRID + j;
                idx2 = (i - item) * FPY0_GRID * FPZ0_GRID + j;

                if (i <= item)
                    array_global[idx] = array_global[idx1];
                if (i >= 2 * item)
                    array_global[idx] = array_global[idx2];
            }

        X_to_distribute_soft (array_global, rho);

        idx1 = FNX_GRID * FPY0_GRID * FPZ0_GRID;
        my_free(array_global);

    }

    /*  Fixed  charge density in some region -----Qingzhong     05/17/03/  */
    if (chargeDensityCompass.type == 2)
    {

        /*  normalize the charge density  */

        tcharge = 0.0;
        for (idx = 0; idx < FP0_BASIS; idx++)
            tcharge += rho[idx];
        ct.tcharge = real_sum_all (tcharge);
        for (idx = 0; idx < FP0_BASIS; idx++)
            rho[idx] = rho_old[idx];


    }

    if (chargeDensityCompass.type == 0)
    {

        /*  normalize the charge density  */

        tcharge = 0.0;
        for (idx = 0; idx < FP0_BASIS; idx++)
            tcharge += rho[idx];
        ct.tcharge = real_sum_all (tcharge) * ct.vel_f;
        if (pct.gridpe == 0)
            printf ("total charge %10.4f = %10.4f + %10.4f\n",
                    ct.tcharge, ct.tcharge - ct.nel, ct.nel);

        t2 = ct.nel / ct.tcharge;

        sscal (&FP0_BASIS, &t2, &rho[0], &ione);
    }

    if (chargeDensityCompass.type == 1)
    {

        x1 = chargeDensityCompass.box1.x1;
        x2 = chargeDensityCompass.box1.x2;
        y1 = chargeDensityCompass.box1.y1;
        y2 = chargeDensityCompass.box1.y2;
        z1 = chargeDensityCompass.box1.z1;
        z2 = chargeDensityCompass.box1.z2;

        pe2xyz (pct.gridpe, &pex, &pey, &pez);
        xoff = pex * FPX0_GRID;
        yoff = pey * FPY0_GRID;
        zoff = pez * FPZ0_GRID;

        total_charge = 0.0;
        tcharge_fixed = 0.0;
        for (i = 0; i < FPX0_GRID; i++)
        {
            for (j = 0; j < FPY0_GRID; j++)
            {
                for (k = 0; k < FPZ0_GRID; k++)
                {
                    idx = i * FPY0_GRID * FPZ0_GRID + j * FPZ0_GRID + k;

                    test = (((i + xoff) < x1) || ((i + xoff) >= x2) ||
                            ((j + yoff) < y1) || ((j + yoff) >= y2) ||
                            ((k + zoff) < z1) || ((k + zoff) >= z2));
                    if (!test)
                    {
                        total_charge += rho[i * FPY0_GRID * FPZ0_GRID + j * FPZ0_GRID + k];

                    }
                    else
                    {
                        rho[idx] = rho_old[idx];
                        tcharge_fixed += rho_old[i * FPY0_GRID * FPZ0_GRID + j * FPZ0_GRID + k];
                    }
                }
            }
        }

        t2 = real_sum_all (total_charge) * ct.vel_f;
        t_fixed = real_sum_all (tcharge_fixed) * ct.vel_f;

        if (pct.gridpe == 0)
            printf ("total charge %10.4f + %10.4f = %10.4f = %10.4f + %10.4f\n",
                    t2, t_fixed, t2 + t_fixed, t2 + t_fixed - ct.nel, ct.nel);

        /*
           t2 = (ct.nel - t_fixed) / (t2 * ct.vel_f);
           t2 = 1.0 / ct.vel_f;

           for (i = 0; i < FPX0_GRID; i++)
           {
           for (j = 0; j < FPY0_GRID; j++)
           {
           for (k = 0; k < FPZ0_GRID; k++)
           {
           idx = i * FPY0_GRID * FPZ0_GRID + j * FPZ0_GRID + k;
           test = (((i + xoff) < chargeDensityCompass.box1.x1)
           || ((i + xoff) >= chargeDensityCompass.box1.x2)
           || ((j + yoff) < chargeDensityCompass.box1.y1)
           || ((j + yoff) >= chargeDensityCompass.box1.y2)
           || ((k + zoff) < chargeDensityCompass.box1.z1)
           || ((k + zoff) >= chargeDensityCompass.box1.z2));
           if (!test)
           {
           rho[i * FPY0_GRID * FPZ0_GRID + j * FPZ0_GRID + k] *= t2;

           }
           }
           }
           }
         */

    }
    my_barrier ();


}
