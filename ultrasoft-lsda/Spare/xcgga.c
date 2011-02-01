/************************** SVN Revision Information **************************
 **    $Id$    **
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


#include "main.h"
#include <float.h>
#include <math.h>


#define    crs    1.91915829267751281
#define    SMALL  1.e-8


void xcgga (REAL * rho, REAL * vxc, REAL * exc, int mode)
{
#if 1
    int ix, iy, iz, idx;
    FP0_GRID *gx, *gy, *gz, *agx, *agy, *agz, *agg, *d2rho;
    REAL d, s, u, v, kf, us, uu, t, vv, ww;
    REAL pisq3, ex, vx, ec;
    REAL zet, rs, g, h, sk, gks2;
    REAL vcup, vcdn;
    REAL dvcup, dvcdn;
    REAL ecrs, eczet, alfc;
    REAL cpot, cen, xen;
    REAL dhalf, d1half[3], d2half, vcm0, fac;
    int ndim = 3, lgga, lpot;

    pisq3 = THREE * PI * PI;


    /* Grab some memory */
    my_malloc (gx, 1, FP0_GRID);
    my_malloc (gy, 1, FP0_GRID);
    my_malloc (gz, 1, FP0_GRID);
    my_malloc (agx, 1, FP0_GRID);
    my_malloc (agy, 1, FP0_GRID);
    my_malloc (agz, 1, FP0_GRID);
    my_malloc (agg, 1, FP0_GRID);
    my_malloc (d2rho, 1, FP0_GRID);






    /* Generate the gradient of the density */
    app_gradf (rho, gx, gy, gz);


    /* Get the Laplacian of the density */
    app6_del2f (rho, d2rho);


    /* Absolute value of grad(rho) */
    for (idx = 0; idx < FP0_BASIS; idx++)
    {

        agg->s2[idx] = sqrt (gx->s2[idx] * gx->s2[idx] +
                             gy->s2[idx] * gy->s2[idx] + gz->s2[idx] * gz->s2[idx]);

    }                           /* end for */




    /* Get its gradient */
    app_gradf (agg->s2, agx, agy, agz);



    /* Now get the potential */
    for (idx = 0; idx < FP0_BASIS; idx++)
    {

        if ((d=rho[idx]) < SMALL && ct.scf_steps < 10)
        {
            d = SMALL;
            fac = exp (50. * (rho[idx] / d - 1.0));
        }
        else
        {
            fac = 1.0;
        }
        kf = pow (pisq3 * d, 0.333333333333333);

        s = agg->s2[idx] / (TWO * kf * d);
        us = gx->s2[idx] * agx->s2[idx] + gy->s2[idx] * agy->s2[idx] + gz->s2[idx] * agz->s2[idx];
        /* 
           us = agg->s2[idx] * d2rho->s2[idx];
         */

        u = us / (d * d * EIGHT * kf * kf * kf);

        v = d2rho->s2[idx] / (FOUR * d * kf * kf);

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
            if ((d = rho[idx]) < 1.e-15)
            {
                cen = 0.0;
                cpot = 0.0;
            }
            else
            {
                if (d < SMALL)
                    d = SMALL;
                dhalf = d / 2.0;
                d1half[0] = gx->s2[idx] / 2.0;
                d1half[1] = gy->s2[idx] / 2.0;
                d1half[2] = gz->s2[idx] / 2.0;
                d2half = d2rho->s2[idx] / 2.0;
                corlyp (&dhalf, &dhalf, d1half, d1half, &d2half, &d2half, &cen, &cpot, &vcm0,
                        &ndim);
            }
            vxc[idx] = fac * (cpot + vx);
            exc[idx] = fac * (cen + ex);
        }
        else if (mode == GGA_PBE)
        {
            zet = ZERO;         /* Spin up = spin down */
            rs = crs / kf;
            sk = TWO * sqrt (kf / PI);

            g = ONE;

            gks2 = TWO * sk * g;


            t = agg->s2[idx] / (d * gks2);
            uu = us / (d * d * gks2 * gks2 * gks2);
            vv = d2rho->s2[idx] / (d * gks2 * gks2);
            ww = ZERO;          /* Non-spin polarized case */

            lpot = 1;
            lgga = 1;
            corpbe (&rs, &zet, &t, &uu, &vv, &ww, &lgga, &lpot, &ec, &vcup, &vcdn,
                    &h, &dvcup, &dvcdn);

            cpot = vcup + dvcup;
            cen = ec + h;

            vxc[idx] = cpot + vx;
            exc[idx] = cen + ex;
            
	    if (d < 1.e-12)
	    { 
		    vxc[idx] = mu_pz(d);
		    exc[idx] = e_pz(d);

	    }

        }
        else
        {
            /* LSD contribution to correlation */
            zet = ZERO;         /* Spin up = spin down */
            rs = crs / kf;
            corlsd (&rs, &zet, &ec, &vcup, &vcdn, &ecrs, &eczet, &alfc);


            sk = TWO * sqrt (kf / PI);

#if 0
            /* commented out for speed since we are doing spin-unpolarized calculations */
            g = pow (1.0 + zet, 0.6666666666666666);
            g += pow (1.0 - zet, 0.6666666666666666);
            g = g / TWO;
#endif
            g = ONE;

            gks2 = TWO * sk * g;


            t = agg->s2[idx] / (d * gks2);
            uu = us / (d * d * gks2 * gks2 * gks2);
            vv = d2rho->s2[idx] / (d * gks2 * gks2);
            ww = ZERO;          /* Non-spin polarized case */

            corgga (&rs, &zet, &t, &uu, &vv, &ww, &h, &dvcup, &dvcdn, &kf, &sk,
                    &g, &ec, &ecrs, &eczet);

            cpot = vcup + dvcup;
            cen = ec + h;

            vxc[idx] = cpot + vx;
            exc[idx] = cen + ex;

        }

    }                           /* end for */
/*
    printf ("\nprint out 50_50 vxc\n");
    print_density_z_direction(50,50,vxc,FPX0_GRID,FPY0_GRID,FPZ0_GRID,20);
*/

    /* Release our memory */
    my_free (d2rho);
    my_free (agg);
    my_free (agz);
    my_free (agy);
    my_free (agx);
    my_free (gz);
    my_free (gy);
    my_free (gx);
#endif

}                               /* end xcgga */
