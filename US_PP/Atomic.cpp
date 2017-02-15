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


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include "transition.h"
#include "GlobalSums.h"
#include "RmgException.h"
#include "Atomic.h"
#include "BesselRoots.h"
#include <boost/math/special_functions/math_fwd.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include "blas.h"




int Atomic::Log_grid_initialized = false;
double Atomic::r_filtered[MAX_LOGGRID];
double Atomic::log_r_filtered[MAX_LOGGRID];
double Atomic::gvec[RADIAL_GVECS];

// Interpolation parameters
double Atomic::a;
double Atomic_inv_a;
double Atomic::b;
double Atomic_inv_b;
double GlogGrid_inv_a;
double GlogGrid_inv_b;

using boost::math::policies::policy;
using boost::math::policies::promote_double;
typedef policy<promote_double<false> > bessel_policy;

// constructor just sets up logarithmic grid for first instance
// which is used by all later instances
Atomic::Atomic(void)
{

    if(!Log_grid_initialized) {

        r_filtered[0] = LOGGRID_START;
        log_r_filtered[0] = log(r_filtered[0]);
        for(int idx = 1;idx < MAX_LOGGRID;idx++) {
            r_filtered[idx] = LOGGRID_MESHPARM * r_filtered[idx-1];
            log_r_filtered[idx] = log(r_filtered[idx]);
        }


        b = log((r_filtered[2] - r_filtered[1]) / (r_filtered[1] - r_filtered[0]));
        Atomic_inv_b = 1.0 / b;
        a = r_filtered[0];
        Atomic_inv_a = 1.0 / a;
        Log_grid_initialized = true;

    }

}

Atomic::~Atomic(void)
{
}

//  This function is used to filter the high frequencies from a radial function 
//  defined on a logarithmic grid. It computes a DFT, applies a cutoff function
//  defined in Gcutoff to the transform and then reconstructs the filtered function
//  on a standard logarithmic grid which is set up in this class.
//
//  The radial function must be short ranged.
//
void Atomic::RftToLogGrid (
                   double cparm,        // IN:  filtering parameter
                   double * f,          // IN:  function to be filtered defined on pseudopotential grid
                   double * r,          // IN:  pseudopotential grid dimensioned r[rg_points], logarithmic
                   double * ffil,       // OUT: filtered function defined on RMG standard log grid
                   double *rab,         // IN:  radial volume element from the pseudopotential file
                   int rg_points,       // IN:  number of points in pseudopotential radial grid
                   int lval,            // IN:  momentum l (0=s, 1=p, 2=d)
                   double width)        // IN:  width of gaussian cutoff in g-space
{


    int istep;
    double t1;
    int npes = Rmg_G->get_PE_X() * Rmg_G->get_PE_Y() * Rmg_G->get_PE_Z();

    /* Get some temporary memory */
    int alloc = rg_points;
    if(alloc < RADIAL_GVECS) alloc = RADIAL_GVECS;
    if (alloc < MAX_LOGGRID)
        alloc = MAX_LOGGRID;
    if (alloc < MAX_RGRID)
        alloc = MAX_RGRID;

    double *work1 = new double[alloc]();
    double *work2 = new double[alloc]();
    double *gcof = new double[alloc]();
    double *gvec = new double[alloc]();

    int gnum = RADIAL_GVECS;

    /* G-vectors are defined on a log grid with the smallest value set */
    /* by the largest real-space value of r.                           */
    gvec[0] = PI / (r[rg_points - 1]);

    /* The smallest g-vector that we generate is defined by the global */
    /* grid spacing.                                                  */
    double gcut = PI / (cparm * ct.hmaxgrid);

    // The largest g-vector we use corresponds to an energy cutoff of 5483 Rydbergs.
    double gmax = PI / 0.03;

    double gmesh = (log (gmax) - log (gvec[0])) / gnum;
    t1 = exp (gmesh);


    /* Generate g-vectors */
    for (int idx = 1; idx < gnum; idx++)
    {

        gvec[idx] = gvec[0] * pow (t1, (double) idx);

    }                           /* end for */


    istep = gnum / npes;
    double alpha = (double)lval + 0.5;

    for (int ift = istep * pct.gridpe; ift < istep * pct.gridpe + istep; ift++)
    {
        for (int idx = 0; idx < rg_points; idx++)
        {

            double jarg = r[idx] * gvec[ift];
            work1[idx] = f[idx] * sqrt(PI/(2.0*jarg)) * boost::math::cyl_bessel_j(alpha, jarg, bessel_policy());

        }                   /* end for */
        gcof[ift] = radint1 (work1, r, rab, rg_points);
    }
    istep = npes * istep;
    GlobalSums (gcof, istep, pct.grid_comm);

    for (int ift = istep; ift < gnum; ift++)
    {
        for (int idx = 0; idx < rg_points; idx++)
        {

            double jarg = r[idx] * gvec[ift];
            work1[idx] = f[idx] * sqrt(PI/(2.0*jarg)) * boost::math::cyl_bessel_j(alpha, jarg, bessel_policy());

        }                   /* end for */
        gcof[ift] = radint1 (work1, r, rab, rg_points);
    }


    /* Fourier Filter the transform and store in work2 */
    for (int idx = 0; idx < gnum; idx++)
    {

        work2[idx] = gcof[idx] * Gcutoff (gvec[idx], gcut, width);

    }                           /* end for */


    /* Zero out the array in which the filtered function is to be returned */
    for (int idx = 0; idx < MAX_LOGGRID; idx++) ffil[idx] = ZERO;

    /* Reconstruct the filtered function onto our standard log grid */
    istep = MAX_LOGGRID / npes;

    for (int idx = istep * pct.gridpe; idx < istep * pct.gridpe + istep; idx++)
    {

        for (int ift = 0; ift < gnum; ift++)
        {

            double jarg = r_filtered[idx] * gvec[ift];
            work1[ift] = work2[ift] * sqrt(PI/(2.0*jarg)) * boost::math::cyl_bessel_j(alpha, jarg, bessel_policy());

        }                   /* end for */


        /* Integrate it */
        t1 = radint (work1, gvec, gnum, gmesh);
        ffil[idx] = t1 * TWO / PI;

    }                       /* end for */


    istep = npes * istep;
    GlobalSums (ffil, istep, pct.grid_comm); 


    for (int idx = istep; idx < MAX_LOGGRID; idx++)
    {

        for (int ift = 0; ift < gnum; ift++)
        {

            double jarg = r_filtered[idx] * gvec[ift];
            work1[ift] = work2[ift] * boost::math::cyl_bessel_j(alpha, jarg, bessel_policy());

        }                   /* end for */

        /* Integrate it */
        t1 = radint (work1, gvec, gnum, gmesh);
        ffil[idx] = t1 * TWO / PI;

    }                       /* end for */


    /* Release memory */
    delete [] gvec;
    delete [] gcof;
    delete [] work2;
    delete [] work1;

} // end RftToLogGrid


// G-vector cutoff function. Applies a gaussian cutoff function above
// the specified cutoff.
double Atomic::Gcutoff (double g1, double gcut, double width)
{

    double t1;

    if (g1 < gcut)
        return ONE;

    t1 = (g1 - gcut) / gcut;
    return exp (-width * (t1 * t1));

}


// This routine computes a finite spherical Bessel function expansion of f, defined on
// the pseudopotential log grid and then uses the coefficients of the expansion to 
// regenerate a filtered version of f on our internal log grid. The filtering is based
// on the number of zeros in each spherical bessel function.
double Atomic::BesselToLogGrid (
                   double cparm,        // IN:  filtering parameter
                   double * f,          // IN:  function to be filtered defined on pseudopotential grid
                   double * r,          // IN:  pseudopotential grid dimensioned r[rg_points], logarithmic
                   double * ffil,       // OUT: filtered function defined on RMG standard log grid
                   double *rab,         // IN:  radial volume element from the pseudopotential file
                   int rg_points,       // IN:  number of points in pseudopotential radial grid
                   int lval,            // IN:  momentum l (0=s, 1=p, 2=d)
                   double rcut,         // IN:  f is defined on [0,rcut] with f(r>=rcut) = 0.0
                   double hmin)         // IN:  target 3D grid spacing

{


    // get minimum grid spacing then find the largest root of the Bessel function such that
    // the normalized distance between any two roots is larger than the minimum grid spacing.
    int N = NUM_BESSEL_ROOTS;
    for(int i = 0;i < NUM_BESSEL_ROOTS;i++) {
        for(int j = 0;j < i;j++) {
            double h = cparm*rcut * (J_roots[lval][i] - J_roots[lval][j]) / J_roots[lval][i];
            if(h <= hmin) {
                N = j;
                goto Bessel1;
                break;
            }
        }
    }
Bessel1:
    //rcut = hmin * J_roots[lval][N] / (cparm * (J_roots[lval][N] - J_roots[lval][N-1])) + hmin/ratio;
    if(ct.localize_projectors)rcut = hmin * J_roots[lval][N] / (cparm * (J_roots[lval][N] - J_roots[lval][N-1]));
    //if(pct.gridpe == 0) printf("Using %d Bessel roots in radial expansion with rcut = %12.6f  hmin = %12.6f\n",N, rcut, hmin);

    // Normalization coefficient
    double JNorm = 2.0 / (rcut*rcut*rcut);

    /* Get some temporary memory */
    int alloc = rg_points;
    if(alloc < RADIAL_GVECS) alloc = RADIAL_GVECS;
    if (alloc < MAX_LOGGRID)
        alloc = MAX_LOGGRID;
    if (alloc < MAX_RGRID)
        alloc = MAX_RGRID;

    double *work1 = new double[alloc]();
    double *bcof = new double[alloc]();

    double alpha = (double)lval;
    for (int i = 0;i < N;i++)
    {
        double JN_i = sqrt(PI/(2.0*J_roots[lval][i])) * boost::math::cyl_bessel_j(alpha + 1.5, J_roots[lval][i]);
        for (int idx = 0; idx < rg_points; idx++)
        {
            double jarg = r[idx] * J_roots[lval][i] / rcut;
            work1[idx] = f[idx] * sqrt(PI/(2.0*jarg)) * boost::math::cyl_bessel_j(alpha + 0.5, jarg);
        }
        bcof[i] = JNorm * radint1 (work1, r, rab, rg_points) / (JN_i * JN_i);
    }


    /* Now we reconstruct the filtered function */
    /* Zero out the array in which the filtered function is to be returned */
    for (int idx = 0; idx < MAX_LOGGRID; idx++)
    {
        ffil[idx] = ZERO;
    }                           /* end for */


    for (int idx = 0;idx < MAX_LOGGRID;idx++)
    {
        if(r_filtered[idx] < rcut) {
            for (int i = 0;i < N;i++)
            {
                double jarg = r_filtered[idx] * J_roots[lval][i] / rcut;
                ffil[idx] += bcof[i] * sqrt(PI/(2.0*jarg)) * boost::math::cyl_bessel_j(alpha + 0.5, jarg);
            }
        }
    }

    /* Release memory */
    delete [] bcof;
    delete [] work1;

    return rcut;

} // end BesselToLogGrid


void Atomic::FilterPotential (
    double *potential,       // IN:  potential to be filtered
    double *r,               // IN:  radial grid for the potential
    int rg_points,           // IN:  number of points in the radial grid
    double rmax,             // IN:  mask function parameter
    double parm,             // IN:  filtering parameter
    double* potential_lgrid, // OUT: potential on standard log grid
    double *rab,             // IN:  radial volume element from the pseudopotential file
    int l_value,             // IN:  momentum l (0=s, 1=p, 2=d)
    double gwidth,           // IN:  gaussian parameter controlling speed of g-space damping
    double rcut,             // IN:  filtered potential is damped in real space starting here
    double rwidth,           // IN:  gaussian parameter controlling speed of real space damping
    double hmin)             // IN:  target 3D grid spacing
{

    // rmax is updated in BesselToLogGrid
    double nrmax = BesselToLogGrid (parm, potential, r, potential_lgrid, rab, rg_points, l_value, rmax, hmin);

    /* Transform to g-space and filter it */
    //RftToLogGrid (parm, potential, r, potential_lgrid, rab, rg_points, l_value, gwidth);

    /*Fix up first point in filtered potential*/
    //potential_lgrid[0] = 2.0 * potential_lgrid[1] - potential_lgrid[2];

    // Damp oscillatory tails in real space
    for (int idx = 0; idx < MAX_LOGGRID; idx++)
    {
        double rdist = r_filtered[idx];

#if 1
        if (rdist > rcut)
        {
//            double t1 = (rdist - rcut) / rcut;
//            double exp_fac = exp (-rwidth * t1 * t1);
//            potential_lgrid[idx] *= exp_fac;
              double t1 = (rdist - rcut) / (nrmax - rcut);
              if(t1 > 1.0) t1 = 1.0;
              potential_lgrid[idx] = (1.0 - t1)*potential_lgrid[idx];
        }               /* end if */

#endif
    }
} 


// Returns a pointer to the shared logarithmic grid
double *Atomic::GetRgrid(void)
{
    return Atomic::r_filtered;
}


// Interpolates function f that is defined on the shared logarithmic grid
// if you modify this remember to modify the inline version in the header file
// AtomicInterpolate.h
double Atomic::Interpolate(double *f, double r)
{
    double d0, d1, dm;
    double f0, g0, g1, g2, h1, h2, i2;
    int ic;

    // truncate to nearest integer index into r_filtered array
    if((r < LOGGRID_START)) {
        r = LOGGRID_START;
        if(fabs(f[0]) < 1.0e-5) return 0.0;
    }

    d0 = (log (r*Atomic_inv_a) * Atomic_inv_b);
    ic = (int)d0;
    ic = (ic > 0) ? ic : 1;

    /* cubic interpolation using forward differences */
    d0 -= (double) (ic);
    d1 = (d0 - 1.0) * 0.5;
    dm = (d0 - 2.0) / 3.0;

    f0 = f[ic];
    g0 = f[ic] - f[ic - 1];
    g1 = f[ic + 1] - f[ic];
    g2 = f[ic + 2] - f[ic + 1];
    h1 = g1 - g0;
    h2 = g2 - g1;
    i2 = h2 - h1;

    return f0 + d0 * (g1 + d1 * (h2 + dm * i2));
    
}


// This function is used to determine the range of the Beta and Q functions in
// real space. It first computes the integral of the square of the function
// and then determines what r includes most of the function.
double Atomic::GetRange(double *f, double *r, double *rab, int rg_points)
{

    int i;

    /* Simpson's rule weights */
    //double w0 = 1.0 / 3.0;
    double w1 = 4.0 / 3.0;
    double w2 = 2.0 / 3.0;
    double *f2 = new double[rg_points];

    for(i = 0;i < rg_points;i++) f2[i] = sqrt(f[i] * f[i]);

    double fsum = radint1 (f2, r, rab, rg_points);

    double tsum = fsum;

    for(i = rg_points - 1;i > 1;i--) {

        tsum -= f2[i] * r[i] * r[i] * rab[i] * (( (i) % 2 == 1 ) ? w1 : w2);
        double ratio = tsum/fsum;
        if(ratio <= 0.999999999) break;
    }

    //if(pct.gridpe==0)printf("BREAK = %d  %20.12f  %20.12f  %20.12f  %20.12f\n",i, fsum, f[rg_points-1], r[i], f[i]);
    return r[i];
    delete [] f2;
}


void Atomic::PackFine2Rhogrid (std::complex<double> *gwptr, int ngrid_fine, std::complex<double> *gbptr, int ngrid) 
{

    //  pack G-spaced atomic local potential, projectors, rhoc etc from fine grid to rho_grid 
    //  gwptr as input, forward function in fine grid
   //   gbptr as output, forward function in rho_grid.

    int idx1, idx2, i, j, k, size;
    int i1, i2, j1, j2, k1, k2;

    size = ngrid_fine * ngrid_fine * ngrid_fine;
    int icut = (ngrid / 2) * (ngrid/ 2)*1000;

    for (i = -ngrid / 2; i <= ngrid / 2; i++)
    {
        int isq = i * i;
        for (j = -ngrid / 2; j <= ngrid / 2; j++)
        {
            int jsq = j * j;
            for (k = -ngrid / 2; k <= ngrid / 2; k++)
            {
                int ksq = k * k;
                if (i < 0)
                {
                    i1 = i + ngrid;
                    i2 = i + ngrid_fine;
                }
                else
                {
                    i1 = i;
                    i2 = i;
                }
                if (j < 0)
                {
                    j1 = j + ngrid;
                    j2 = j + ngrid_fine;
                }
                else
                {
                    j1 = j;
                    j2 = j;
                }
                if (k < 0)
                {
                    k1 = k + ngrid;
                    k2 = k + ngrid_fine;
                }
                else
                {
                    k1 = k;
                    k2 = k;
                }
                idx1 = i1 * ngrid * ngrid + j1 * ngrid + k1;
                idx2 = i2 * ngrid_fine * ngrid_fine + j2 * ngrid_fine + k2;
                
                if(icut > (isq + jsq + ksq)) {
                    gbptr[idx1] =  gwptr[idx2] / (double) size;
                }
                else {
                    gbptr[idx1] = 0.0 ;
                }
            }
        }
    }

}



void Atomic::RLogGridToGLogGrid (
                   double * f,          // IN:  function to be filtered defined on pseudopotential grid
                   double * r,          // IN:  pseudopotential grid dimensioned r[rg_points], logarithmic
                   double *rab,         // IN:  radial volume element from the pseudopotential file
                   double * f_g,        // OUT: on radial g grid
                   int rg_points,       // IN:  number of points in pseudopotential radial grid
                   int lval,
                   double *bessel_rg)       
{



    double *work1 = new double[rg_points]();

    int gnum = RADIAL_GVECS;


    for(int ift = 0; ift< gnum; ift++) f_g[ift] = 0.0;

    double *bessel_rg_l = &bessel_rg[lval * RADIAL_GVECS * rg_points ] ;
    for (int ift =  0; ift < gnum; ift++)
    {
        for(int idx = 0; idx < rg_points; idx++)
            work1[idx] = f[idx] * bessel_rg_l[ift * rg_points + idx];
        f_g[ift] = 4.0 * PI * radint1 (work1, r, rab, rg_points);
    
    }


    /* Fourier Filter the transform and store in work2 */
    for (int idx = 0; idx < gnum; idx++)
    {

    //    printf("\n %e %e adad", gvec[idx], f_g[idx]);
        // f_g[idx] = gcof[idx] * Gcutoff (gvec[idx], gcut, width);

    }                           /* end for */



    /* Release memory */
    delete [] work1;

} // end RftToLogGrid

void Atomic::InitBessel(
                   double * r,       // IN:  pseudopotential grid dimensioned r[rg_points], logarithmic
                   int rg_points,    // IN:  number of points in pseudopotential radial grid
                   int lmax,
                   double *bessel_rg)   // Bessel_rg store lmax * RADIAL_GVECS * rg_points
{


    double t1;

    /* Get some temporary memory */
    int alloc = rg_points;
    if(alloc < RADIAL_GVECS) alloc = RADIAL_GVECS;
    if (alloc < MAX_LOGGRID)
        alloc = MAX_LOGGRID;
    if (alloc < MAX_RGRID)
        alloc = MAX_RGRID;

    int gnum = RADIAL_GVECS;

    /* G-vectors are defined on a log grid with the smallest value set */
    /* by the largest real-space value of r.                           */
    gvec[0] = LOGGRID_START;
    // The largest g-vector we use corresponds to an energy cutoff of 5483 Rydbergs.
    double gmax = PI / 0.03;

    double gmesh = (log (gmax) - log (gvec[0])) / gnum;
    t1 = exp (gmesh);


    /* Generate g-vectors */
    for (int idx = 1; idx < gnum; idx++)
    {
        gvec[idx] = gvec[idx-1] * t1;
    }                           /* end for */

    GlogGrid_inv_b = log((gvec[2] - gvec[1]) / (gvec[1] - gvec[0]));
    GlogGrid_inv_b = 1.0 / GlogGrid_inv_b;
    GlogGrid_inv_a = 1.0 / gvec[0];

    double alpha ;

    for(int idx = 0; idx < (lmax+1) * gnum * rg_points; idx++) bessel_rg[idx] = 0.0;
    for(int lval = 0; lval <= lmax; lval++)
    {
        alpha = (double)lval + 0.5;
        for (int ift = pct.gridpe; ift < gnum; ift+=pct.grid_npes)
        {
            int idx1 = lval * gnum * rg_points + ift * rg_points;
            for (int idx = 0; idx < rg_points; idx++)
            {

                double jarg = r[idx] * gvec[ift];
                bessel_rg[idx1 + idx] = sqrt(PI/(2.0*jarg)) * boost::math::cyl_bessel_j(alpha, jarg, bessel_policy());

            }                   /* end for */
        }
    }

    int size = (lmax+1) * gnum * rg_points; 
    GlobalSums (bessel_rg, size, pct.grid_comm); 


} 
