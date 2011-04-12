/************************** SVN Revision Information **************************
 **    $Id$    **
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
 *   REAL fill (STATE *states, REAL width, REAL nel, REAL mix,
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


static REAL fd (REAL mu, REAL * occ, STATE * states, REAL width, REAL nel, int num_st);
static REAL gs (REAL mu, REAL * occ, STATE * states, REAL width, REAL nel, int num_st);
static REAL ef (REAL mu, REAL * occ, STATE * states, REAL width, REAL nel, int num_st);



REAL fill (STATE * states, REAL width, REAL nel, REAL mix, int num_st, int occ_flag)
{

    const int maxit = 50;
    const REAL charge_tol = 1.0e-10;

    int iter, st, st1, kpt, idx, nks, nspin = (ct.spin_flag + 1);
    STATE *sp;
    REAL mu, dmu, mu1, mu2, f, fmid, eig;
    REAL (*func) () = NULL;
    REAL *occ;


    /* bracket the roots and bisect to find the mu */

    switch (occ_flag % 10)
    {

    case OCC_NONE:
        return 0.0;             /* early return */

    case OCC_FD:
        /* fermi-dirac occupations: f(x) = 2 / (1 + Exp[x/T]) */
        func = fd;
        break;

    case OCC_GS:
        /* Gaussian occupations:
           f(x) = 1 - sign(x) (1 - Exp[-|x|/(8T)(4 +|x|/T)^2]) */
        func = gs;
        break;

    case OCC_EF:
        /* error-function occupations: f(x) = erfc(x/(aT)),
           where a = 4/Sqrt[Pi] */
        func = ef;
        break;

    default:
        error_handler ("unknown filling procedure");

    }                           /* end switch */
    
    nks = ct.num_kpts * ct.num_states;
    my_malloc (occ, nspin * nks, REAL);


    /* find the root by bisection: this algorithm was adapted
       from numerical recipes, 2nd edition */

    /* find max and min eigenvalues */

    /* debug: doesn't handle the case of all energies being equal */

    mu1 = 1.0e30;
    mu2 = -mu1; 

    for(idx = 0; idx < nspin; idx++)
    {
    	for (kpt = 0; kpt < ct.num_kpts; kpt++)
    	{
        	for (st = 0; st < ct.num_states; st++)
    		{
            		/* sp = &ct.kp[kpt].kstate[st]; */
	    		eig = ct.kp[kpt].kstate[st].eig[idx]; 
            		mu1 = min (eig, mu1);
            		mu2 = max (eig, mu2); 
		}
    	}                           
    } 


    fmid = func (mu2, occ, states, width, nel, num_st);
    f = func (mu1, occ, states, width, nel, num_st); 

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
        fmid = func (mu1, occ, states, width, nel, num_st);

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
        printf ("\nfill: \\sum f - n_el= %e", fmid);
        error_handler ("did not converge");
    }                           /* end if */

    /* mix occupations */

    fmid = 0.0; 
	
    for (idx = 0; idx < nspin; idx++)
    {
    	st = -1;
    	for (kpt = 0; kpt < ct.num_kpts; kpt++)
   	{
        	for (st1 = 0; st1 < ct.num_states; st1++)
        	{
            		st += 1;
            		sp = &ct.kp[kpt].kstate[st1];
	    
           	 	sp->occupation[idx] = mix * occ[st + idx * nks] + (1.0 - mix) * sp->occupation[idx];
            		fmid += sp->occupation[idx] * ct.kp[kpt].kweight;
        	}
    	}                           /* st and kpt */
    }


    fmid -= nel;

    if (fabs (fmid) > charge_tol)
    {
        printf ("\nfill: \\sum f - n_el= %e", fmid);
        error_handler ("error in mixing occupations");
    }                           /* end if */

    my_free (occ);

    return (mu);

}                               /* end fill */ 



static REAL fd (REAL mu, REAL * occ, STATE * states, REAL width, REAL nel, int num_st)
{
    int st, kpt, st1, idx, nks, nspin = (ct.spin_flag + 1);
    STATE *sp;
    REAL t1, t2, sumf, eig, fac = (2.0 - ct.spin_flag);

    /* fermi-dirac occupations:
       f(x) = 2 / (1 + Exp[x/T]) */

    nks = ct.num_kpts * ct.num_states;
    
    sumf = 0.0; 
    
    for (idx = 0; idx < nspin; idx++)
    {
    	st = -1;
    	for (kpt = 0; kpt < ct.num_kpts ; kpt++)
    	{
		for (st1 = 0; st1 < ct.num_states; st1++)
    		{
			st += 1;
	    		eig = ct.kp[kpt].kstate[st1].eig[idx]; 
        		t1 = (eig - mu) / width;
        		if (t1 > 0.0)
        		{

        			t2 = exp (-t1);
                		occ[st + idx * nks] = fac * t2 / (1.0 + t2);
            
       		 	}
        		else
        		{
                		t2 = exp (t1);
               		 	occ[st + idx * nks] = fac / (1.0 + t2);
        		}                 
       	        	sumf += occ[st + idx * nks] * ct.kp[kpt].kweight;
        	}
    	}                           /* st1 and kpt */

    }    

    return (sumf - nel);

}                               /* fd */


static REAL gs (REAL mu, REAL * occ, STATE * states, REAL width, REAL nel, int num_st)
{
    int st, kpt, st1, nks, idx, nspin = (ct.spin_flag + 1);
    STATE *sp;
    REAL t1, sumf, eig, fac;
    
    fac = (2.0 - ct.spin_flag) * 0.5;

    /* Gaussian occupations:
       f(x) = 1 - sign(x) (1 - Exp[-|x|/(8T)(4 +|x|/T)^2]) */
    
    nks = ct.num_kpts * ct.num_states;

    sumf = 0.0; 
 
    for (idx = 0; idx < nspin; idx++)
    {
    	st = -1; 
    	for (kpt = 0; kpt < ct.num_kpts; kpt++)
    	{
        	for (st1 = 0; st1 < ct.num_states; st1++)
        	{
            	st += 1;
            	eig = ct.kp[kpt].kstate[st1].eig[idx];

            	t1 = (eig - mu) / (2.0 * width);
            	if (t1 > 0.0)
            	{
                	occ[st + idx * nks] = fac * exp (-t1 * (1.0 + 0.5 * t1));
            	}
            	else
            	{
                	t1 = -t1;
                	occ[st + idx * nks] = fac * ( 2.0 - exp (-t1 * (1.0 + 0.5 * t1)) );

            	}                   /* end if */
            	sumf += occ[st + idx * nks] * ct.kp[kpt].kweight;
		}
    	}
    }	

    return (sumf - nel);

}                               /* gs */


static REAL ef (REAL mu, REAL * occ, STATE * states, REAL width, REAL nel, int num_st)
{
    int st, kpt, st1, nks, idx, nspin = (ct.spin_flag + 1);
    STATE *sp;
    REAL t1, t2, sumf, eig, fac;

    fac = (2.0 - ct.spin_flag) * 0.5;

    t2 = 4.0 / sqrt (PI);
    nks = ct.num_kpts * ct.num_states; 

    /* error-function occupations: f(x) = erfc(x/(aT))
       where a = 4/Sqrt[Pi] */

    sumf = 0.0; 

    for (idx = 0; idx < nspin; idx++)
    {
    	st = -1;
    	for (kpt = 0; kpt < ct.num_kpts; kpt++)
    	{
        	for (st1 = 0; st1 < ct.num_states; st1++)
        	{
            		st += 1;
            		eig = ct.kp[kpt].kstate[st1].eig[idx];

            		t1 = (eig - mu) / (t2 * width);

            		/* debug: this conditional mayn't be necessary */

            		if (t1 > 0.0)
            		{
                		occ[st + idx * nks] = fac * erfc (t1);
            		}
            		else
            		{
                		t1 = -t1;       
                		occ[st + idx * nks] = fac * ( 1.0 + erf (t1) );
            		}                   /* end if */

           	 	sumf += occ[st + idx * nks] * ct.kp[kpt].kweight;
		}
    	}
    }	    

    return (sumf - nel);

}                               /* ef */

/******/
