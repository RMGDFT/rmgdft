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


static double occ_allstates (double mu, double * occ, double *eigs, double width, double nel, 
    int num_st, double *weight, int occ_flag, int mp_order);
static inline double dist_func(double t1, int occ_flag, int mp_order);


template double Fill (Kpoint<double> **, double, double, double, int, int, int);
template double Fill (Kpoint<std::complex<double> > **, double, double, double, int, int, int);


    template <typename KpointType>
double Fill (Kpoint<KpointType> **Kptr, double width, double nel, double mix, int num_st, int occ_flag, int mp_order)
{
    if(ct.adaptive_convergence)
    {
    //    width = 10.0*width*exp(-0.05*ct.scf_steps) + width;
        double steps = (double)ct.scf_steps/40.0;
        width = 0.5*eV_Ha*exp(-steps*steps) + width;
    }

    if(ct.occ_flag == OCC_TETRA) return Tetra->FillTetra(Kptr);

    const int maxit = 100;
    const double charge_tol = 1.0e-10;

    int iter, st, st1;
    State<KpointType> *sp;
    double mu, dmu, mu1, mu2, f, fmid;

    int nks = ct.num_kpts_pe * ct.num_states;
    int ntot_states = nks;
    if(ct.nspin == 2) ntot_states *=2;

    double *eigs = new double[ntot_states]();
    double *weight = new double[ct.num_kpts_pe];
    double *occ = new double[ntot_states]();

    for(int kpt = 0;kpt < ct.num_kpts_pe;kpt++) weight[kpt] = Kptr[kpt]->kp.kweight;

    // Methfessel Paxton can be unstable on the first step when eigenvalues may be crazy
    // so use FD in that case.
    if((ct.scf_steps < 2) && (ct.exx_steps == 0) && (occ_flag == OCC_MP))
    {
        if((ct.runflag != RESTART) && (ct.runflag != Restart_TDDFT)) occ_flag = OCC_FD;
    }

    // Fill eigs
    double Evbm = -1000.0;
    for(int kpt = 0;kpt < ct.num_kpts_pe;kpt++)
    {
        for(int st = 0;st < ct.num_states;st++)
        {
            eigs[kpt*ct.num_states + st] = Kptr[kpt]->Kstates[st].eig[0];
            if(ct.nspin == 2)
                eigs[nks + kpt*ct.num_states + st] = Kptr[kpt]->Kstates[st].eig[1];
            if(Kptr[kpt]->Kstates[st].occupation[0] > 1.0e-5)
                Evbm = std::max(Evbm, Kptr[kpt]->Kstates[st].eig[0]);
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

    if(occ_flag == OCC_NONE ) 
    {
        MPI_Allreduce (MPI_IN_PLACE, &Evbm, 1, MPI_DOUBLE, MPI_MAX, pct.img_comm);
        return Evbm;
    }



    /* find the root by bisection: this algorithm was adapted
       from numerical recipes, 2nd edition */

    /* find max and min eigenvalues */

    /* debug: doesn't handle the case of all energies being equal */

    mu1 = 1.0e30;
    mu2 = -mu1; 

    for(int idx = 0; idx < ntot_states; idx++)
    {
        mu1 = std::min (eigs[idx], mu1);
        mu2 = std::max (eigs[idx], mu2); 
    } 

    double mu1_tem = mu1;    
    double mu2_tem = mu2;    
    MPI_Allreduce (&mu1_tem, &mu1, 1, MPI_DOUBLE, MPI_MIN, pct.kpsub_comm);
    MPI_Allreduce (&mu2_tem, &mu2, 1, MPI_DOUBLE, MPI_MAX, pct.kpsub_comm);


    fmid = occ_allstates (mu2, occ, eigs, width, nel, num_st, weight, occ_flag, mp_order);
    f = occ_allstates (mu1, occ, eigs, width, nel, num_st, weight, occ_flag, mp_order); 

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
        fmid = occ_allstates (mu1, occ, eigs, width, nel, num_st, weight, occ_flag, mp_order);

        if (fmid <= 0.0)
        {
            mu = mu1;
            fmid = -fmid;
        }                       /* end if */

    }
    while ((iter < maxit) && (std::abs(fmid) > charge_tol));

    if (iter == maxit)
        rmg_error_handler (__FILE__,__LINE__,"too many bisections");

    if (fabs (fmid) > charge_tol)
    {
        printf ("\nfill: \\sum f - n_el= %e", fmid);
        rmg_error_handler (__FILE__,__LINE__,"did not converge");
    }                           /* end if */

    /* mix occupations */

    fmid = 0.0; 

    for(int kpt = 0;kpt < ct.num_kpts_pe;kpt++) 
    {
        for (st1 = 0; st1 < ct.num_states; st1++)
        {
            st = kpt * ct.num_states + st1;
            sp = &Kptr[kpt]->Kstates[st1];

            sp->occupation[0] = mix * occ[st] + (1.0 - mix) * sp->occupation[0];
            fmid += sp->occupation[0] * Kptr[kpt]->kp.kweight;
            if(ct.nspin == 2)
            {
                sp->occupation[1] = mix * occ[st + nks] + (1.0 - mix) * sp->occupation[1];
                fmid += sp->occupation[1] * Kptr[kpt]->kp.kweight;
            }
        }
    }                           /* st and kpt */

    fmid = RmgSumAll(fmid, pct.kpsub_comm);

    fmid -= nel;

    if (fabs (fmid) > charge_tol * 10)
    {
        rmg_printf ("\nfill: \\sum f - n_el= %e", fmid);
        rmg_printf ("error in mixing occupations fmid = %e", fmid);
        rmg_error_handler(__FILE__, __LINE__, "Terminating.\n");
    }                           /* end if */

    delete [] occ;
    delete [] weight;
    delete [] eigs;

    return (mu1);

}                               /* end fill */ 



static double occ_allstates (double mu, double * occ, double *eigs, double width, double nel, 
        int num_st, double *weight, int occ_flag, int mp_order)
{
    int st, kpt, st1, nks;
    double t1, sumf, eig, fac = 1.0;

    if(ct.nspin == 1) fac = 2.0;
    /* fermi-dirac occupations:
       f(x) = 2 / (1 + Exp[x/T]) */

    nks = ct.num_kpts_pe * ct.num_states;

    sumf = 0.0; 

    for (kpt = 0; kpt < ct.num_kpts_pe ; kpt++)
    {
        for (st1 = 0; st1 < ct.num_states; st1++)
        {
            st = kpt * ct.num_states + st1;

            eig = eigs[st];
            t1 = (eig - mu) / width;
            occ[st] = fac * dist_func(t1, occ_flag, mp_order);
            sumf += occ[st] * weight[kpt];

            if(ct.nspin == 2)
            {
                eig = eigs[st + nks];
                t1 = (eig - mu) / width;
                occ[st + nks ] = fac * dist_func(t1, occ_flag, mp_order);
                sumf += occ[st + nks] * weight[kpt];
            }
        }
    }                           /* st1 and kpt */


    sumf = RmgSumAll(sumf, pct.kpsub_comm);
    return (sumf - nel);

}                               /* fd */




static inline double dist_func(double t1, int occ_flag, int mp_order)
{
    double t2, oc, An, h0, h1, hexp;
    switch (occ_flag % 10)
    {

        case OCC_FD:
            if (t1 > 0.0)
            {
                t2 = exp (-t1);
                return t2 / (1.0 + t2);
            }
            else
            {
                t2 = exp (t1);
                return 1.0 / (1.0 + t2);
            }                 
            /* fermi-dirac occupations: f(x) = 2 / (1 + Exp[x/T]) */
            break;

        case OCC_GS:
            /* Gaussian occupations:
               f(x) = 1 - sign(x) (1 - Exp[-|x|/(8T)(4 +|x|/T)^2]) */

            t1 = 0.5 * t1;
            if (t1 > 0.0)
            {
                return exp (-t1 * (1.0 + 0.5 * t1));
            }
            else
            {
                t1 = -t1;
                return ( 2.0 - exp (-t1 * (1.0 + 0.5 * t1)) );

            }                   /* end if */
            break;

        case OCC_EF:
            /* error-function occupations: f(x) = erfc(x/(aT)),
               where a = 4/Sqrt[Pi] */

            t2 = 4.0 / sqrt (PI);
            t1 = t1/ t2;

            return erfc (t1);

            break;
        case OCC_MV:
            t2 = -t1 -1.0/sqrt(2.0);
            return 0.5 * erfc(-t2) + 1.0/sqrt(2.0*PI) * exp(-t2*t2);
            break;

        case OCC_MP:
            oc = 0.5 * erfc(t1);
            hexp = exp(-t1*t1);
            An = 1.0/sqrt(PI);
            h0 = 1.0;
            h1 = 2.0 * t1;
            for(int n = 1; n <= ct.mp_order; n++)
            {
                An = -An/(4.0 * n);
                oc += An * h1 * hexp;
                h0 = 2.0 * t1 * h1 - 2.0 * (2*n-1) * h0;    // h0 = H_2n(x)
                h1 = 2.0* t1 * h0 - 2.0 * (2.0*n) * h1; // h1 = H_2n+1(x)
            }

            return oc; 
            break;

        default:
            rmg_error_handler (__FILE__,__LINE__,"unknown filling procedure");

    }                           /* end switch */

    return 0.0;
}
