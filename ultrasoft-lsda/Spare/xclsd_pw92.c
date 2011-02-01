/************************** SVN Revision Information **************************
 **    $Id: xcgga.c 1142 2010-09-13 21:30:11Z btan $    **
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
/*#define    SMALL  1.e-8*/
/*#define    SMALL  1.e-15 */
#define    SMALL  1.e-15


/* for spin polarized system, return xcpotential for up and down*/
void xclsd_pw92(REAL * rho_up, REAL * rho_dn, REAL * vxc_up, REAL * exc)
{
#if 1
    int ix, iy, iz, idx;

    REAL d, kf, ecrs, eczet, alfc;
   
    REAL pisq3, ex, vx, ec, fac;
    REAL zet, rs, g, h, sk, gks2;
    REAL vc_up, vc_dn;
    REAL cen, xen;
    REAL cpot_up, cpot_dn;
    
    REAL vx_up, vx_dn;
    REAL ex_up, ex_dn;

    REAL rho2, vxc_dn[FP0_BASIS], charge;
   
    pisq3 = THREE * PI * PI;

    /* Now get the potential */
    for (idx = 0; idx < FP0_BASIS; idx++)
    {
	d = rho_up[idx] + rho_dn[idx];    
        
	//if (d < SMALL && ct.scf_steps < 50)
        if ( (rho_up[idx] < SMALL) || (rho_dn[idx] < SMALL))
	{   
            vxc_up[idx]=0.0;
	    vxc_dn[idx]=0.0;
	    exc[idx]=0.0;
            continue;
        }
        
	kf = pow (pisq3 * d, 0.333333333333333);

        /* Use Ex[up,dn] = 0.5 * (Ex[2*up] + Ex[2*dn]) */
         
        /* do spin up part first */ 
        rho2 = TWO * rho_up[idx];
        exchlsd (&rho2, &ex_up, &vx_up);
         
        /* Now do the spin down part */   
        rho2 = TWO * rho_dn[idx];
        exchlsd (&rho2, &ex_dn, &vx_dn);
         
        /* get the exchage energy */   
        ex = (ex_up * rho_up[idx] + ex_dn * rho_dn[idx]) / d;            /*(2.0 * d)*/         
        

        /* relative polarization for spin polarized */
        zet = (rho_up[idx]-rho_dn[idx])/d;
        rs = crs / kf;
        sk = TWO * sqrt (kf / PI);

        /* for spin polarized calculation */
        g = pow (1.0 + zet, 0.6666666666666666);
        g += pow (1.0 - zet, 0.6666666666666666);
        g = g / TWO;
            
        gks2 = TWO * sk * g;

        corlsd (&rs, &zet, &ec, &vc_up, &vc_dn, &ecrs, &eczet, &alfc);


        /* corelation potential for spin up and down respectively*/
        cpot_up = vc_up ;
        cpot_dn = vc_dn ;           

        cen = ec ;

        /* get exchange-correlation potential for spin up and down*/
        vxc_up[idx] = cpot_up + vx_up;
        vxc_dn[idx] = cpot_dn + vx_dn;

        exc[idx] = cen + ex;
	//dprintf("grid %d: rho %.4e and exchange energy %.4e and correlation energy %.4e\n", idx, d, ex, cen);

    }
  	    /* end for */


    /* check whether the summation over rho gives the total electron*/
    charge = 0.0 ;
    for (idx = 0; idx < FP0_BASIS; idx++)
	    charge += rho_up[idx];

    charge = real_sum_all(charge);
    charge = charge * ct.vel_f;


    idx = minimum(vxc_up);
    dprintf("minimum  grid %d:  charge %.4e, rho %.4e and exchange-correlation potential %.4e\n", idx, charge,  rho_up[idx]+rho_dn[idx], *(vxc_up+idx));
 
    idx = maximum(vxc_up);
    dprintf("maximum  grid %d:  charge %.4e, rho %.4e and exchange-correlation potential %.4e\n", idx, charge,  rho_up[idx]+rho_dn[idx], *(vxc_up+idx));

#endif
}                               /* end xclsd_pw92 */
