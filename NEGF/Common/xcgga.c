/************************** SVN Revision Information **************************
 **    $Id: xcgga.c 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 
/****f* QMD-MGDFT/xcgga.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5
 * COPYRIGHT
 *   Copyright (C) 1995  Emil Briggs
 *   Copyright (C) 1998  Emil Briggs, Charles Brabec, Mark Wensell, 
 *                       Dan Sullivan, Chris Rapcewicz, Jerzy Bernholc
 *   Copyright (C) 2001  Emil Briggs, Wenchang Lu,
 *                       Marco Buongiorno Nardelli,Charles Brabec, 
 *                       Mark Wensell,Dan Sullivan, Chris Rapcewicz,
 *                       Jerzy Bernholc
 * FUNCTION
 *   void xcgga(P0_GRID *rho, P0_GRID *vxc, P0_GRID *exc, int mode)
 *   Functions to generate the exchange-correlation potential and
 *   energy density using the generalized gradient approximation.
 * INPUTS
 *   rho: total valence charge density
 *   mode: type of GGA 
 *         0 = LDA-PZ81
 *         1 = GGA-BLYP
 *         2 = GGA-XB-CP
 *         3 = GGA-XP-CP
 *         4 = GGA-PBE
 * OUTPUT
 *   vxc: exchange correlation potential
 *   exc: exchange correlation energy  
 * PARENTS
 *   get_te.c get_vxc.c
 * CHILDREN
 *   pack_ptos.c trade_images.c app_grad.c app6_del2.c    
 * SOURCE
 */


#include "md.h"
#include <float.h>
#include <math.h>


#define    crs    1.91915829267751281
#define    SMALL  1.e-08


void xcgga (REAL * rho1, REAL * vxc, REAL * exc, int mode)
{
    int ix, iy, iz, idx;
    REAL *gx, *gy, *gz, *agx, *agy, *agz, *agg, *d2rho, *nrho;
    REAL d, s, u, v, kf, us, uu, t, vv, ww;
    REAL pisq3, ex, vx, ec;
    REAL zet, rs, g, h, sk, gks2;
    REAL vcup, vcdn;
    REAL dvcup, dvcdn;
    REAL ecrs, eczet, alfc;
    REAL cpot, cen, xen;
    REAL dhalf, d1half[3], d2half, vcm0, fac;
    int ndim, lgga, lpot;


    pisq3 = 3.0 * PI * PI;


    /* Grab some memory */
    my_malloc_init( gx, 8 * FP0_BASIS, REAL );
    gy = gx + FP0_BASIS;
    gz = gy + FP0_BASIS;
    agx = gz + FP0_BASIS;
    agy = agx + FP0_BASIS;
    agz = agy + FP0_BASIS;
    agg = agz + FP0_BASIS;
    d2rho = agg + FP0_BASIS;


    /* Load the density into the smoothing grid */
    nrho = rho1;

    /* Generate the gradient of the density */
    app_grad (rho1, gx, gy, gz, FPX0_GRID, FPY0_GRID, FPZ0_GRID);


    /* Get the Laplacian of the density */
    app6_del2 (rho1, d2rho, FPX0_GRID, FPY0_GRID, FPZ0_GRID, ct.hxxgrid, ct.hyygrid, ct.hzzgrid);


    /* Absolute value of grad(rho) */
    for (idx = 0; idx < FP0_BASIS; idx++)
    {

        agg[idx] = sqrt (gx[idx] * gx[idx] + gy[idx] * gy[idx] + gz[idx] * gz[idx]);

    }                           /* end for */


    /* Get its gradient */
    app_grad (agg, agx, agy, agz, FPX0_GRID, FPY0_GRID, FPZ0_GRID);



    /* Now get the potential */
    ndim = 3;
    for (idx = 0; idx < FP0_BASIS; idx++)
    {

        if ((d = nrho[idx]) < SMALL)
        {
            d = SMALL;
            fac = exp (50. * (nrho[idx] / d - 1.0));
        }
        else
        {
            fac = 1.0;
        }
        kf = pow (pisq3 * d, 0.333333333333333);

        s = agg[idx] / (2.0 * kf * d);
        us = gx[idx] * agx[idx] + gy[idx] * agy[idx] + gz[idx] * agz[idx];
        /* 
           us = agg[idx] * d2rho[idx];
         */

        u = us / (d * d * 8.0 * kf * kf * kf);

        v = d2rho[idx] / (4.0 * d * kf * kf);

        if (mode == GGA_BLYP || mode == GGA_XB_CP)
        {

            /* Exchange potential from becke */
            xbecke (&d, &s, &u, &v, &ex, &vx);

        }
        else if (mode == GGA_PBE)
        {

            /* exchange potential from Perdew, Burke, Ernzerhof */
            lpot = 1;
            lgga = 1;
            exchpbe (&d, &s, &u, &v, &lpot, &lgga, &ex, &vx);
        }
        else if (mode == GGA_XP_CP)
        {

            /* Exchange potential from Perdew */
            exch (&d, &s, &u, &v, &ex, &vx);

        }                       /* end if */


        if (mode == GGA_BLYP)
        {
            if ((d = nrho[idx]) < 1.e-15)
            {
                cen = 0.0;
                cpot = 0.0;
            }
            else
            {
                if (d < SMALL)
                    d = SMALL;
                dhalf = d / 2.0;
                d1half[0] = gx[idx] / 2.0;
                d1half[1] = gy[idx] / 2.0;
                d1half[2] = gz[idx] / 2.0;
                d2half = d2rho[idx] / 2.0;
                corlyp(&dhalf, &dhalf, d1half, d1half, &d2half, &d2half, &cen, &cpot, &vcm0,&ndim);
            }
            vxc[idx] = fac * (cpot + vx);
            exc[idx] = fac * (cen + ex);
        }
        else if (mode == GGA_PBE)
        {
            zet = 0.0;         /* Spin up = spin down */
            rs = crs / kf;
            sk = 2.0 * sqrt (kf / PI);

            g = 1.0;

            gks2 = 2.0 * sk * g;


            t = agg[idx] / (d * gks2);
            uu = us / (d * d * gks2 * gks2 * gks2);
            vv = d2rho[idx] / (d * gks2 * gks2);
            ww = 0.0;          /* Non-spin polarized case */

            lpot = 1;
            lgga = 1;
            corpbe (&rs, &zet, &t, &uu, &vv, &ww, &lgga, &lpot, &ec, &vcup, &vcdn,
                    &h, &dvcup, &dvcdn);

            cpot = vcup + dvcup;
            cen = ec + h;

            vxc[idx] = cpot + vx;
            exc[idx] = cen + ex;

        }
        else
        {
            /* LSD contribution to correlation */
            zet = 0.0;         /* Spin up = spin down */
            rs = crs / kf;
            corlsd (&rs, &zet, &ec, &vcup, &vcdn, &ecrs, &eczet, &alfc);


            sk = 2.0 * sqrt (kf / PI);

#if 0
            /* commented out for speed since we are doing spin-unpolarized calculations */
            g = pow (1.0 + zet, 0.6666666666666666);
            g += pow (1.0 - zet, 0.6666666666666666);
            g = g / 2.0;
#endif
            g = 1.0;

            gks2 = 2.0 * sk * g;


            t = agg[idx] / (d * gks2);
            uu = us / (d * d * gks2 * gks2 * gks2);
            vv = d2rho[idx] / (d * gks2 * gks2);
            ww = 0.0;          /* Non-spin polarized case */

            corgga (&rs, &zet, &t, &uu, &vv, &ww, &h, &dvcup, &dvcdn, &kf, &sk,
                    &g, &ec, &ecrs, &eczet);

            cpot = vcup + dvcup;
            cen = ec + h;

            vxc[idx] = cpot + vx;
            exc[idx] = cen + ex;

        }

    }                           /* end for */


    /* Release our memory */
    my_free(gx);

}                               /* end xcgga */
