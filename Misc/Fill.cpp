/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

/****f* QMD-MGDFT/fill.c *****
 *
 * FUNCTION
 *   double fill (STATE *states, double width, double nel, double mix,
 *              int num_st, int occ_flag)
 *   Sets the occupations of the electronic orbitals and stored in
 *   ct.kp[kpt].Kstate[st].occupation 
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
 *    and error function distribution, respeget_vel()y.
 * SOURCE
 */

#include <complex>
#include <math.h>
#include "const.h"
#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "State.h"
#include "transition.h"
#include "RmgSumAll.h"
#include "GlobalSums.h"


static double fd (double mu, double * occ, double *eigs, double width, double nel, int num_st, double *weight);
static double gs (double mu, double * occ, double *eigs, double width, double nel, int num_st, double *weight);
static double ef (double mu, double * occ, double *eigs, double width, double nel, int num_st, double *weight);


template double Fill (Kpoint<double> **, double, double, double, int, int);
template double Fill (Kpoint<std::complex<double> > **, double, double, double, int, int);


    template <typename KpointType>
double Fill (Kpoint<KpointType> **Kptr, double width, double nel, double mix, int num_st, int occ_flag)
{

    const int maxit = 50;
    const double charge_tol = 1.0e-10;

    int iter, st, st1, kpt, idx, nks, nspin = (ct.spin_flag + 1);
    State<KpointType> *sp;
    double mu, dmu, mu1, mu2, f, fmid, eig;
    double (*func) (double, double *, double *, double, double, int, double*) = NULL;

    double *occ;
    double eigs[ct.num_states*ct.num_kpts_pe*nspin];
    double weight[ct.num_kpts_pe];
    for(int kpt = 0;kpt < ct.num_kpts_pe;kpt++) 
        weight[kpt] = Kptr[kpt]->kweight;



    nks = ct.num_kpts_pe * ct.num_states;
    int ntot_states = nspin * nks;

    for(int idx = 0; idx < ntot_states; idx++) eigs[idx] = 0.0;

    // Fill eigs
    for(int ispin = 0; ispin<nspin; ispin++)
        for(int kpt = 0;kpt < ct.num_kpts_pe;kpt++) {
            for(int st = 0;st < ct.num_states;st++) {
                eigs[ispin * nks + kpt*ct.num_states + st] = Kptr[kpt]->Kstates[st].eig[ispin];
            }
        }




    if(nel == 1 && ct.num_kpts == 1 && ct.spin_flag == 0)
    {
        sp = &Kptr[0]->Kstates[0];
        sp->occupation[0] = 1.0;
        mu = sp->eig[0];
        for (st1 = 1; st1 < ct.num_states; st1++)
        {
            sp = &Kptr[0]->Kstates[st1];
            sp->occupation[0] = 0.0;
        }

        return(mu);
    }

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
            rmg_error_handler (__FILE__,__LINE__,"unknown filling procedure");

    }                           /* end switch */

    occ = new double[nspin * nks];


    /* find the root by bisection: this algorithm was adapted
       from numerical recipes, 2nd edition */

    /* find max and min eigenvalues */

    /* debug: doesn't handle the case of all energies being equal */

    mu1 = 1.0e30;
    mu2 = -mu1; 

    for(idx = 0; idx < ntot_states; idx++)
    {
        mu1 = std::min (eigs[idx], mu1);
        mu2 = std::max (eigs[idx], mu2); 
    } 

    double mu1_tem = mu1;    
    double mu2_tem = mu2;    
    MPI_Allreduce (&mu1_tem, &mu1, 1, MPI_DOUBLE, MPI_MIN, pct.kpsub_comm);
    MPI_Allreduce (&mu2_tem, &mu2, 1, MPI_DOUBLE, MPI_MAX, pct.kpsub_comm);


    fmid = func (mu2, occ, eigs, width, nel, num_st, weight);
    f = func (mu1, occ, eigs, width, nel, num_st, weight); 

    if (f * fmid >= 0.0)
        rmg_error_handler (__FILE__, __LINE__, "root must be bracketed");

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
        fmid = func (mu1, occ, eigs, width, nel, num_st, weight);

        if (fmid <= 0.0)
        {
            mu = mu1;
            fmid = -fmid;
        }                       /* end if */

    }
    while ((iter < maxit) && (fmid > charge_tol));

    if (iter == maxit)
        rmg_error_handler (__FILE__,__LINE__,"too many bisections");

    if (fabs (fmid) > charge_tol)
    {
        printf ("\nfill: \\sum f - n_el= %e", fmid);
        rmg_error_handler (__FILE__,__LINE__,"did not converge");
    }                           /* end if */

    /* mix occupations */

    fmid = 0.0; 

    for (idx = 0; idx < nspin; idx++)
    {
        st = -1;
        for(int kpt = 0;kpt < ct.num_kpts_pe;kpt++) 
        {
            for (st1 = 0; st1 < ct.num_states; st1++)
            {
                st = kpt * ct.num_states + st1;
                sp = &Kptr[kpt]->Kstates[st1];

                sp->occupation[idx] = mix * occ[st + idx * nks] + (1.0 - mix) * sp->occupation[idx];
                fmid += sp->occupation[idx] * Kptr[kpt]->kweight;
            }
        }                           /* st and kpt */
    }

    fmid = RmgSumAll(fmid, pct.kpsub_comm);

    fmid -= nel;

    if (fabs (fmid) > charge_tol * 10)
    {
        rmg_printf ("\nfill: \\sum f - n_el= %e", fmid);
        rmg_printf ("error in mixing occupations fmid = %e", fmid);
        rmg_error_handler(__FILE__, __LINE__, "Terminating.\n");
    }                           /* end if */

    delete [] occ;

    return (mu1);

}                               /* end fill */ 



static double fd (double mu, double * occ, double *eigs, double width, double nel, int num_st, double *weight)
{
    int st, kpt, st1, idx, nks, nspin = (ct.spin_flag + 1);
    double t1, t2, sumf, eig, fac = (2.0 - ct.spin_flag);

    /* fermi-dirac occupations:
       f(x) = 2 / (1 + Exp[x/T]) */

    nks = ct.num_kpts_pe * ct.num_states;

    sumf = 0.0; 

    for (idx = 0; idx < nspin; idx++)
    {
        st = -1;
        for (kpt = 0; kpt < ct.num_kpts_pe ; kpt++)
        {
            for (st1 = 0; st1 < ct.num_states; st1++)
            {
                st = kpt * ct.num_states + st1;
                eig = eigs[st + idx * nks];
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
                sumf += occ[st + idx * nks] * weight[kpt];
            }
        }                           /* st1 and kpt */

    }    

    sumf = RmgSumAll(sumf, pct.kpsub_comm);

    return (sumf - nel);

}                               /* fd */


static double gs (double mu, double * occ, double *eigs, double width, double nel, int num_st, double *weight)
{
    int st, kpt, st1, nks, idx, nspin = (ct.spin_flag + 1);
    double t1, sumf, eig, fac;

    fac = (2.0 - ct.spin_flag) * 0.5;

    /* Gaussian occupations:
       f(x) = 1 - sign(x) (1 - Exp[-|x|/(8T)(4 +|x|/T)^2]) */

    nks = ct.num_kpts_pe * ct.num_states;

    sumf = 0.0; 

    for (idx = 0; idx < nspin; idx++)
    {
        st = -1; 
        for (kpt = 0; kpt < ct.num_kpts_pe; kpt++)
        {
            for (st1 = 0; st1 < ct.num_states; st1++)
            {
                st = kpt * ct.num_states + st1;
                eig = eigs[idx * nks + st];
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
                sumf += occ[st + idx * nks] * weight[kpt];
            }
        }
    }	

    sumf = RmgSumAll(sumf, pct.kpsub_comm);
    return (sumf - nel);

}                               /* gs */


static double ef (double mu, double * occ, double *eigs, double width, double nel, int num_st, double *weight)
{
    int st, kpt, st1, nks, idx, nspin = (ct.spin_flag + 1);
    double t1, t2, sumf, eig, fac;

    fac = (2.0 - ct.spin_flag) * 0.5;

    t2 = 4.0 / sqrt (PI);
    nks = ct.num_kpts_pe * ct.num_states; 

    /* error-function occupations: f(x) = erfc(x/(aT))
       where a = 4/Sqrt[Pi] */

    sumf = 0.0; 

    for (idx = 0; idx < nspin; idx++)
    {
        for (kpt = 0; kpt < ct.num_kpts_pe; kpt++)
        {
            for (st1 = 0; st1 < ct.num_states; st1++)
            {
                st = kpt * ct.num_states + st1;
                eig = eigs[st + idx * nks];

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

                sumf += occ[st + idx * nks] * weight[kpt];
            }
        }
    }	    

    sumf = RmgSumAll(sumf, pct.kpsub_comm);
    return (sumf - nel);

}                               /* ef */


