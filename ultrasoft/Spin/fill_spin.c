/************************** SVN Revision Information **************************
 **    $Id: fill.c 1066 2009-08-31 18:41:09Z froze $    **
******************************************************************************/

/****f* QMD-MGDFT/fill.c *****
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
 *   REAL fill_spin (STATE *states, REAL width, REAL nel, REAL mix,
 *              int num_st, int occ_flag)
 *   Sets the occupations of the electronic orbitals and stored in
 *   ct.kp[kpt].kstate[st].occupation 
 * INPUTS
 *   states: points to orbital structure
 *   width:  Width of distribution for variable occupations
 *   nel:    Number of electrons
 *   mix:    Linear mixing parameter for the occupations
 *   num_st: Number of states, it is useless now
 *   occ_flag: Type of occupation scheme
 * OUTPUT
 *   Fermi energy is returned
 * PARENTS
 *   scf.c 
 * CHILDREN
 *   functions fd, gs, ef are for Fermi-dirac, Gaussian, 
 *    and error function distribution, respectively.
 * SOURCE
 */

#include <math.h>
#include <stdio.h>
#include "main.h"


static REAL fd_spin (REAL mu, REAL * occ, REAL * occ_oppo, STATE * states, REAL width, REAL nel, int num_st);
static REAL gs_spin (REAL mu, REAL * occ, REAL * occ_oppo, STATE * states, REAL width, REAL nel, int num_st);
static REAL ef_spin (REAL mu, REAL * occ, REAL * occ_oppo, STATE * states, REAL width, REAL nel, int num_st);


/*Spin up and down are just name of processor's own spin and opposite spin for easy implementation,
   not the real spin up and down*/

REAL fill_spin (STATE * states, REAL width, REAL nel, REAL mix, int num_st, int occ_flag)
{
    const int maxit = 50;
    const REAL charge_tol = 1.0e-10;

    int iter, st, st1, kpt;
    STATE *sp;
    REAL mu, dmu, mu1, mu2, f, fmid, eigtemp;
    REAL (*func) () = NULL;
    REAL *occ, *occ_oppo;  /* Occupations calculated from occupation distribution function */


    /* bracket the roots and bisect to find the mu */

    switch (occ_flag % 10)
    {

    case OCC_NONE:
        return 0.0;             /* early return */

    case OCC_FD:
        /* fermi-dirac occupations: f(x) = 2 / (1 + Exp[x/T]) */
        func = fd_spin;
        break;

    case OCC_GS:
        /* Gaussian occupations:
           f(x) = 1 - sign(x) (1 - Exp[-|x|/(8T)(4 +|x|/T)^2]) */
        func = gs_spin;
        break;

    case OCC_EF:
        /* error-function occupations: f(x) = erfc(x/(aT)),
           where a = 4/Sqrt[Pi] */
        func = ef_spin;
        break;

    default:
        error_handler ("unknown filling procedure");

    }                           /* end switch */

    my_malloc (occ, ct.num_kpts * ct.num_states, REAL);
    my_malloc (occ_oppo, ct.num_kpts * ct.num_states_oppo, REAL);
   

    /* find the root by bisection: this algorithm was adapted
       from numerical recipes, 2nd edition */

    /* find max and min eigenvalues */

    /* debug: doesn't handle the case of all energies being equal */

    mu1 = 1.0e30;
    mu2 = -mu1;

    /* own spin */
    for (kpt = 0; kpt < ct.num_kpts; kpt++)
    {
        for (st = 0; st < ct.num_states; st++)
        {
            /* sp = &ct.kp[kpt].kstate[st];  */
	    eigtemp = ct.kp[kpt].kstate[st].eig; 

            mu1 = (eigtemp < mu1) ? eigtemp : mu1;
            mu2 = (eigtemp > mu2) ? eigtemp : mu2;

        }
    }                           /* st and kpt */
    //printf("\n");
    
    /* opposite spin */
    for (kpt = 0; kpt < ct.num_kpts; kpt++)
    {
        for (st = 0; st < ct.num_states_oppo; st++)
    {
            /* sp = &ct.kp[kpt].kstate[st]; */
	    eigtemp = ct.kp[kpt].kstate[st].eig_oppo; 

            mu1 = (eigtemp < mu1) ? eigtemp : mu1;
            mu2 = (eigtemp > mu2) ? eigtemp : mu2;
	    //printf("\n%.4e\t", eigtemp);
    //printf("\n");

        }
    }                           /* st and kpt */
    
    fmid = func (mu2, occ, occ_oppo, states, width, nel, num_st);
    f = func (mu1, occ, occ_oppo, states, width, nel, num_st);

    if (f * fmid >= 0.0)
        error_handler ("root must be bracketed");

    if (f < 0.0)
    {
        mu = mu1;
        dmu = mu2 - mu1;
    }
    else
    {
        mu = mu2;
        dmu = mu1 - mu2;
    }                           /* end if */

    iter = 0;
    do
    {
        iter++;

        dmu *= 0.5;
        mu1 = mu + dmu;
        fmid = func (mu1, occ, occ_oppo, states, width, nel, num_st);

        if (fmid <= 0.0)
        {
            mu = mu1;
            fmid = -fmid;
        }                       /* end if */

    }
    while ((iter < maxit) && (fmid > charge_tol));

    if (iter == maxit)
        error_handler ("too many bisections");

    if (fabs (fmid) > charge_tol)
    {
        if (pct.imgpe == 0)
            printf ("\nfill_spin: \\sum f - n_el= %e", fmid);
        error_handler ("did not converge");
    }                           /* end if */ 


    /* mix occupations */

    fmid = 0.0;
    st = -1;

    /* own spin */
    for (kpt = 0; kpt < ct.num_kpts; kpt++)
    {
        for (st1 = 0; st1 < ct.num_states; st1++)
        {
            st += 1;
            sp = &ct.kp[kpt].kstate[st1];
	    
            sp->occupation = mix * occ[st] + (1.0 - mix) * sp->occupation;
            fmid += sp->occupation * ct.kp[kpt].kweight;
        }
    }                           /* st and kpt */

    /* opposite spin */
    st = -1;
    for (kpt = 0; kpt < ct.num_kpts; kpt++)
    {
        for (st1 = 0; st1 < ct.num_states_oppo; st1++)
        {
            st += 1;
            sp = &ct.kp[kpt].kstate[st1];
	    
            sp->occupation_oppo = mix * occ_oppo[st] + (1.0 - mix) * sp->occupation_oppo;
            fmid += sp->occupation_oppo * ct.kp[kpt].kweight;
        }
    }                           /* st and kpt */

    fmid -= nel; 


    if (fabs (fmid) > charge_tol)
    {
        if (pct.imgpe == 0)
            printf ("\nfill_spin: \\sum f - n_el= %e", fmid);
       // error_handler ("error in mixing occupations");
    }                           /* end if */

    my_free (occ);
    my_free (occ_oppo);

    return (mu);
}
                               /* end fill */

static REAL fd_spin (REAL mu, REAL * occ, REAL * occ_oppo, STATE * states, REAL width, REAL nel, int num_st)
{
    int st, kpt, st1;
    STATE *sp;
    REAL t1, t2, sumf, eigtemp;

    /* fermi-dirac occupations:
       f(x) = 2 / (1 + Exp[x/T]) */

    sumf = 0.0;
    st = -1;

   /* own spin */
    for (kpt = 0; kpt < ct.num_kpts; kpt++)
    {
        for (st1 = 0; st1 < ct.num_states; st1++)
        {
            st += 1;
            eigtemp = ct.kp[kpt].kstate[st1].eig;

            t1 = (eigtemp - mu) / width;

            if (t1 > 0.0)
            {
                t2 = exp (-t1);
                occ[st] = t2 / (1.0 + t2);
            }
            else
            {
                t2 = exp (t1);
                occ[st] = 1.0 / (1.0 + t2);

            }                   /* end if */

            sumf += occ[st] * ct.kp[kpt].kweight;

        }
    }                           /* st1 and kpt */

   /* opposite spin */
    st = -1;
    for (kpt = 0; kpt < ct.num_kpts ; kpt++)
    {
	for (st1 = 0; st1 < ct.num_states_oppo; st1++)
    	{
		st += 1;
	    	eigtemp = ct.kp[kpt].kstate[st1].eig_oppo; 
        	t1 = (eigtemp - mu) / width;
        	if (t1 > 0.0)
        	{

        		t2 = exp (-t1);
                	occ_oppo[st] = t2 / (1.0 + t2);
            
       		 }
        	else
        	{
                	t2 = exp (t1);
               		 occ_oppo[st] = 1.0 / (1.0 + t2);
        	}                   /* end if */
       	        sumf += occ_oppo[st] * ct.kp[kpt].kweight;
        }
    }                           /* st1 and kpt */

    return (sumf - nel);
}                               /* fd_spin */


static REAL gs_spin (REAL mu, REAL * occ, REAL * occ_oppo, STATE * states, REAL width, REAL nel, int num_st)
{
    int st, kpt, st1;
    STATE *sp;
    REAL t1, sumf, eigtemp;

    /* Gaussian occupations:
       f(x) = 1 - sign(x) (1 - Exp[-|x|*(1 + 0.5 * |x|)]) */

    sumf = 0.0;
    st = -1;

    /* own spin */
    for (kpt = 0; kpt < ct.num_kpts; kpt++)
    {
        for (st1 = 0; st1 < ct.num_states; st1++)
        {
            st += 1;
            eigtemp = ct.kp[kpt].kstate[st1].eig;

            t1 = (eigtemp - mu) / (2.0 * width);
            if (t1 > 0.0)
            {
                occ[st] = 0.5 * exp (-t1 * (1.0 + 0.5 * t1));
            }
            else
            {
                t1 = -t1;
                occ[st] = 1.0 - 0.5 * exp (-t1 * (1.0 + 0.5 * t1));

            }                   /* end if */

            sumf += occ[st] * ct.kp[kpt].kweight;

        }
    }                           /* st1 and kpt */


    /* opposite spin */
    st = -1; 
    for (kpt = 0; kpt < ct.num_kpts; kpt++)
    {
        for (st1 = 0; st1 < ct.num_states_oppo; st1++)
        {
            st += 1;
            eigtemp = ct.kp[kpt].kstate[st1].eig_oppo;

            t1 = (eigtemp - mu) / (2.0 * width);
            if (t1 > 0.0)
            {
                occ_oppo[st] = 0.5 * exp (-t1 * (1.0 + 0.5 * t1));
            }
            else
            {
                t1 = -t1;
                occ_oppo[st] = 1.0 - 0.5 * exp (-t1 * (1.0 + 0.5 * t1));

            }                   /* end if */
            sumf += occ_oppo[st] * ct.kp[kpt].kweight;
        }
    }                           /* st1 and kpt */


    return (sumf - nel);
}                               /* gs_spin */



static REAL ef_spin (REAL mu, REAL * occ, REAL * occ_oppo, STATE * states, REAL width, REAL nel, int num_st)
{
    int st, kpt, st1;
    STATE *sp;
    REAL t1, t2, sumf, eigtemp;

    t2 = 4.0 / sqrt (PI);

    /* error-function occupations: f(x) = erfc(x/(aT))
       where a = 4/Sqrt[Pi] */

    sumf = 0.0;

    /* own spin */
    st = -1;
    for (kpt = 0; kpt < ct.num_kpts; kpt++)
    {
        for (st1 = 0; st1 < ct.num_states; st1++)
        {
            st += 1;
            eigtemp = ct.kp[kpt].kstate[st1].eig;

            t1 = (eigtemp - mu) / (t2 * width);

            /* debug: this conditional mayn't be necessary */

            if (t1 > 0.0)
            {
                occ[st] = 0.5 * erfc (t1);
            }
            else
            {
                t1 = -t1;       /* |t1| */
                occ[st] = 0.5 + 0.5 * erf (t1);
            }                   /* end if */

            sumf += occ[st] * ct.kp[kpt].kweight;

        }
    }                           /* st1 and kpt */


    /* opposite spin */
    st = -1;
    for (kpt = 0; kpt < ct.num_kpts; kpt++)
    {
        for (st1 = 0; st1 < ct.num_states_oppo; st1++)
        {
            st += 1;
            eigtemp = ct.kp[kpt].kstate[st1].eig_oppo;

            t1 = (eigtemp - mu) / (t2 * width);

            /* debug: this conditional mayn't be necessary */

            if (t1 > 0.0)
            {
                occ_oppo[st] = 0.5 * erfc (t1);
            }
            else
            {
                t1 = -t1;       /* |t1| */
                occ_oppo[st] = 0.5 + 0.5 * erf (t1);
            }                   /* end if */

            sumf += occ_oppo[st] * ct.kp[kpt].kweight;

        }
    }                           /* st1 and kpt */


    return (sumf - nel);

}                               /* ef_spin */

/******/
