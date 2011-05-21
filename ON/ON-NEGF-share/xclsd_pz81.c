/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*

                            xclsd_pz81.c


    Functions to generate the exchange-correlation potential and
    energy density using the Ceperley - Alder exchange-correlation
    potential and energy as parameterized by Perdew and Zunger, 
    Phys. Rev. B23, 5048 (1981)

    This is the local spin-density version emulating Karsten Friehube's
    LSD XC routines from the Scheffler code. MGW 10/8/96.
*/



/* spin unpolarized constants */
#define		A_U		 0.0311
#define		B_U		-0.0480
#define		gamma_U		-0.1423
#define		beta1_U		 1.0529
#define		beta2_U		 0.3334
#define		C_U		 0.0020
#define		D_U		-0.0116

/* spin polarized constants */
#define		A_P		 0.01555
#define		B_P		-0.0269
#define		gamma_P		-0.0843
#define		beta1_P		 1.3981
#define		beta2_P		 0.2611
#define		C_P		 0.0007
#define		D_P		-0.0048

/* other constants */
#define		EPS		 0.0000001
#define		TWO_TO_THIRD	 1.2599210498948731648
#define		F_DENOMINATOR    0.5198420997897463296
#define		DF_DENOMINATOR   0.3898815748423097472


#include "md.h"
#include <float.h>
#include <math.h>
#include <stdio.h>


/* Computes exchange correlation energy density */
/* rho_up & rho_dn are the normalized realspace */
/* spin electron charge densities (in Hartrees) */
REAL e_lsd_pz(REAL rho_up, REAL rho_dn)
{

    REAL rho, zeta, rs;
    REAL sqrt_rs, log_rs;
    REAL Ex_U, Ec_U, Ex_P, Ec_P, Exc_U, Exc_P, Exc;
    REAL Bpr_U, Bpr_P, f;
    REAL tmp, tmp1, tmp2;

    /* get total charge */
    rho = rho_up + rho_dn;

    /* make sure rho_up >= rho_dn */
    if (rho_dn > rho_up)
    {
        tmp = rho_dn;
        rho_dn = rho_up;
        rho_up = tmp;
    }

    /* formulate zeta */
    if (rho > ZERO)
    {
        zeta = (rho_up - rho_dn) / rho;
    }
    else
    {
        zeta = ZERO;
    }

    /* check to see if density is valid */
    if (rho < -0.1)
    {

        error_handler("negative total density");

    }
    else if (rho_up < ZERO || rho_dn < ZERO)
    {

        if (rho_up < ZERO)
            rho_up = ZERO;
        if (rho_dn < ZERO)
            rho_dn = ZERO;
        rho = rho_up + rho_dn;

        if (rho > ZERO)
        {
            zeta = (rho_up - rho_dn) / rho;
        }
        else
        {
            zeta = ZERO;
        }                       /*end of if */

    }
    else if (zeta - EPS > ONE || zeta < -EPS)
    {

        tmp = rho_up - rho_dn;
        tmp = sqrt(tmp * tmp);

        if (tmp < EPS)
        {
            zeta = ZERO;
        }
        else
        {
            printf("\n Warning in e_lsd_xc: bad zeta: %f", zeta);
            if (zeta - EPS > ONE)
            {
                zeta = ONE;
            }
            else
            {
                zeta = ZERO;
            }
        }                       /*end of if */

    }                           /* end of if */


    if (rho < 1.0e-10)
    {

        return ZERO;

    }                           /* end of if */


    /* get rs and friends */
    rs = (3.0 / (rho * 4.0 * PI));
    rs = pow(rs, 0.3333333333333333333333333e0);
    sqrt_rs = sqrt(rs);


    /* get the unpolarized exchange energy */
    tmp = NINE * PI / FOUR;
    tmp = -(THREE / (FOUR * PI)) * pow(tmp, 0.3333333333333333333333333e0);
    Ex_U = tmp / rs;


    /* get the unpolarized correlation energy */
    if (rs >= ONE)
    {

        Bpr_U = ONE + beta1_U * sqrt_rs + beta2_U * rs;
        Ec_U = gamma_U / Bpr_U;

    }
    else
    {                           /* if rs < 1 */

        log_rs = log(rs);
        Ec_U = A_U * log_rs + B_U + C_U * rs * log_rs + D_U * rs;

    }                           /* end of if */


    /* now get the polarized exchange energy */
    Ex_P = Ex_U * TWO_TO_THIRD;


    /* get the polarized correlation energy */
    if (rs >= ONE)
    {

        Bpr_P = ONE + beta1_P * sqrt_rs + beta2_P * rs;
        Ec_P = gamma_P / Bpr_P;

    }
    else
    {                           /* if rs < 1 */

        Ec_P = A_P * log_rs + B_P + C_P * rs * log_rs + D_P * rs;

    }                           /* end of if */


    tmp = zeta - ONE;
    tmp = sqrt(tmp * tmp);
    if (tmp < EPS)
        zeta = ONE;

    tmp = sqrt(zeta * zeta);
    if (tmp < EPS)
        zeta = ZERO;

    /* calculate interpolation function f */
    if (zeta == ZERO)
    {
        f = ZERO;
    }
    else if (zeta == ONE)
    {
        f = ONE;
    }
    else
    {
        tmp1 = ONE + zeta;
        tmp1 = pow(tmp1, 1.33333333333333333333e0);
        tmp2 = ONE - zeta;
        tmp2 = pow(tmp2, 1.33333333333333333333e0);
        f = (tmp1 + tmp2 - TWO) / F_DENOMINATOR;
    }

    Exc_U = Ex_U + Ec_U;
    Exc_P = Ex_P + Ec_P;
    Exc = Exc_U + f * (Exc_P - Exc_U);


    return Exc;

}                               /* end of e_lsd_pz */


/* Computes spin up & dn exchange correlation potential */
/* rho_up & rho_dn are the normalized realspace         */
/* spin electron charge densities (in Hartrees)         */
void mu_lsd_pz(REAL rho_up, REAL rho_dn, REAL * Vxc_up, REAL * Vxc_dn)
{

    REAL rho, zeta, rs;
    REAL sqrt_rs, log_rs;
    REAL Ex_U, Ec_U, Ex_P, Ec_P, Exc_U, Exc_P;
    REAL Vx_U, Vc_U, Vx_P, Vc_P, Vxc_U, Vxc_P, VVxc;
    REAL Bpr_U, Bpr_P, f, df, E;
    REAL tmp, tmp1, tmp2;

    /* get total charge */
    rho = rho_up + rho_dn;

    /* make sure rho_up >= rho_dn */
    if (rho_dn > rho_up)
    {
        tmp = rho_dn;
        rho_dn = rho_up;
        rho_up = tmp;
    }

    /* formulate zeta */
    if (rho > ZERO)
    {
        zeta = (rho_up - rho_dn) / rho;
    }
    else
    {
        zeta = ZERO;
    }

    /* check to see if density is valid */
    if (rho < -0.1)
    {

        error_handler("negative total density");

    }
    else if (rho_up < ZERO || rho_dn < ZERO)
    {

        if (rho_up < ZERO)
            rho_up = ZERO;
        if (rho_dn < ZERO)
            rho_dn = ZERO;

        rho = rho_up + rho_dn;

        if (rho > ZERO)
        {
            zeta = (rho_up - rho_dn) / rho;
        }
        else
        {
            zeta = ZERO;
        }                       /*end of if */

    }
    else if (zeta - EPS > ONE || zeta < -EPS)
    {

        tmp = rho_up - rho_dn;
        tmp = sqrt(tmp * tmp);

        if (tmp < EPS)
        {
            zeta = ZERO;
        }
        else
        {
            printf("\n Warning in e_lsd_xc: bad zeta: %f", zeta);
            if (zeta - EPS > ONE)
            {
                zeta = ONE;
            }
            else
            {
                zeta = ZERO;
            }
        }                       /*end of if */

    }                           /* end of if */


    if (rho < 1.0e-10)
    {

        *Vxc_up = ZERO;
        *Vxc_dn = ZERO;
        return;

    }                           /* end of if */


    /* get rs and friends */
    rs = (3.0 / (rho * 4.0 * PI));
    rs = pow(rs, 0.3333333333333333333333333e0);
    sqrt_rs = sqrt(rs);


    /* get the unpolarized exchange energy */
    tmp = NINE * PI / FOUR;
    tmp = -(THREE / (FOUR * PI)) * pow(tmp, 0.3333333333333333333333333e0);
    Ex_U = tmp / rs;
    Vx_U = (FOUR / THREE) * Ex_U;


    /* get the unpolarized correlation energy */
    if (rs >= ONE)
    {

        Bpr_U = ONE + beta1_U * sqrt_rs + beta2_U * rs;
        Ec_U = gamma_U / Bpr_U;
        tmp = ONE + (SEVEN / SIX) * beta1_U * sqrt_rs + (FOUR / THREE) * beta2_U * rs;
        Vc_U = Ec_U * tmp / Bpr_U;

    }
    else
    {                           /* if rs < 1 */

        log_rs = log(rs);
        Ec_U = A_U * log_rs + B_U + C_U * rs * log_rs + D_U * rs;
        Vc_U = Ec_U - (A_U + C_U * rs * (log_rs + ONE) + D_U * rs) / THREE;

    }                           /* end of if */


    /* now get the polarized exchange energy */
    Ex_P = Ex_U * TWO_TO_THIRD;
    Vx_P = (FOUR / THREE) * Ex_P;


    /* get the polarized correlation energy */
    if (rs >= ONE)
    {

        Bpr_P = ONE + beta1_P * sqrt_rs + beta2_P * rs;
        Ec_P = gamma_P / Bpr_P;
        tmp = ONE + (SEVEN / SIX) * beta1_P * sqrt_rs + (FOUR / THREE) * beta2_P * rs;
        Vc_P = Ec_P * tmp / Bpr_P;

    }
    else
    {                           /* if rs < 1 */

        Ec_P = A_P * log_rs + B_P + C_P * rs * log_rs + D_P * rs;
        Vc_P = Ec_P - (A_P + C_P * rs * (log_rs + ONE) + D_P * rs) / THREE;

    }                           /* end of if */


    tmp = zeta - ONE;
    tmp = sqrt(tmp * tmp);
    if (tmp < EPS)
        zeta = ONE;

    tmp = sqrt(zeta * zeta);
    if (tmp < EPS)
        zeta = ZERO;

    if (zeta == ZERO)
    {
        f = ZERO;
    }
    else if (zeta == ONE)
    {
        f = ONE;
    }
    else
    {
        tmp1 = ONE + zeta;
        tmp1 = pow(tmp1, 1.33333333333333333333e0);
        tmp2 = ONE - zeta;
        tmp2 = pow(tmp2, 1.33333333333333333333e0);
        f = (tmp1 + tmp2 - TWO) / F_DENOMINATOR;
    }

    Exc_U = Ex_U + Ec_U;
    Exc_P = Ex_P + Ec_P;

    Vxc_U = Vx_U + Vc_U;
    Vxc_P = Vx_P + Vc_P;

    VVxc = Vxc_U + f * (Vxc_P - Vxc_U);

    if (zeta == ZERO)
    {
        df = ZERO;
    }
    else
    {
        tmp1 = ONE + zeta;
        tmp1 = pow(tmp1, 0.3333333333333333333e0);
        tmp2 = ONE - zeta;
        tmp2 = pow(tmp2, 0.3333333333333333333e0);
        df = (tmp1 - tmp2) / DF_DENOMINATOR;
    }

    E = (Exc_P - Exc_U) * df;
    *Vxc_up = VVxc + (ONE - zeta) * E;
    *Vxc_dn = VVxc - (ONE + zeta) * E;

}                               /* end of mu_lsd_pz */



void xclsd_pz81(REAL * rho, REAL * vxc)
{

    int i;
    REAL *rho_up, *rho_dn;
    REAL *vxc_up, *vxc_dn;

    /* set up the pointers to the up and down rho and vxc */
    rho_up = &rho[0];
    rho_dn = &rho[P0_BASIS];

    vxc_up = &vxc[0];
    vxc_dn = &vxc[P0_BASIS];

    for (i = 0; i < P0_BASIS; i++)
    {

        mu_lsd_pz(rho_up[i], rho_dn[i], &vxc_up[i], &vxc_dn[i]);

    }                           /* end for */

    my_barrier();


}                               /* end xclsda_pz81 */




void exclsd_pz81(REAL * rho, REAL * exc)
{

    int i;
    REAL *rho_up, *rho_dn;

    /* set up the pointers to the up and down rho and vxc */
    rho_up = &rho[0];
    rho_dn = &rho[P0_BASIS];

    for (i = 0; i < P0_BASIS; i++)
    {

        exc[i] = e_lsd_pz(rho_up[i], rho_dn[i]);

    }                           /* end for */

    my_barrier();


}                               /* end exclsda_pz81 */
