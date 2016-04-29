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





int Atomic::Log_grid_initialized = false;
double Atomic::r_filtered[MAX_LOGGRID];
double Atomic::log_r_filtered[MAX_LOGGRID];

// Interpolation parameters
double Atomic::a;
double Atomic_inv_a;
double Atomic::b;
double Atomic_inv_b;

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
    double t1, t2;
    int npes = Rmg_G->get_PE_X() * Rmg_G->get_PE_Y() * Rmg_G->get_PE_Z();
    const double small = 1.e-35;

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

    /* The largest g-vector that we generate is defined by the global */
    /* grid spacing.                                                  */
    double gcut = PI / (cparm * ct.hmaxgrid);
    double gmax = PI / r[0];

    double gmesh = (log (gmax) - log (gvec[0])) / gnum;
    t1 = exp (gmesh);


    /* Generate g-vectors */
    for (int idx = 1; idx < gnum; idx++)
    {

        gvec[idx] = gvec[0] * pow (t1, (double) idx);

    }                           /* end for */


    istep = gnum / npes;
    switch (lval)
    {

        case S_STATE:

            for (int ift = istep * pct.gridpe; ift < istep * pct.gridpe + istep; ift++)
            {
                for (int idx = 0; idx < rg_points; idx++)
                {

                    t1 = r[idx] * gvec[ift];

                    if (t1 > small)
                        work1[idx] = f[idx] * sin (t1) / t1;
                    else
                        work1[idx] = f[idx];
                    /*t2 = sin(t1) / t1;
                       if(t1 < 1.0e-8) t2 = 1.0;
                       work1[idx] = f[idx] * t2; */

                }                   /* end for */
                gcof[ift] = radint1 (work1, r, rab, rg_points);
            }
            istep = npes * istep;
            GlobalSums (gcof, istep, pct.grid_comm);

            for (int ift = istep; ift < gnum; ift++)
            {
                for (int idx = 0; idx < rg_points; idx++)
                {

                    t1 = r[idx] * gvec[ift];

                    if (t1 > small)
                        work1[idx] = f[idx] * sin (t1) / t1;
                    else
                        work1[idx] = f[idx];
                    /*t2 = sin(t1) / t1;
                       if(t1 < 1.0e-8) t2 = 1.0;
                       work1[idx] = f[idx] * t2; */

                }                   /* end for */
                gcof[ift] = radint1 (work1, r, rab, rg_points);
            }
            break;

        case P_STATE:

            for (int ift = istep * pct.gridpe; ift < istep * pct.gridpe + istep; ift++)
            {
                for (int idx = 0; idx < rg_points; idx++)
                {

                    t1 = r[idx] * gvec[ift];
                    if (t1 > small)
                        t2 = cos (t1) / t1 - sin (t1) / (t1 * t1);
                    else
                        t2 = 0.0;
                    /*if(t1 < 1.0e-8) t2 = 0.0; */
                    work1[idx] = f[idx] * t2;

                }                   /* end for */
                gcof[ift] = radint1 (work1, r, rab, rg_points);
            }
            istep = npes * istep;
            GlobalSums (gcof, istep, pct.grid_comm);

            for (int ift = istep; ift < gnum; ift++)
            {
               for (int idx = 0; idx < rg_points; idx++)
                {

                    t1 = r[idx] * gvec[ift];
                    if (t1 > small)
                        t2 = cos (t1) / t1 - sin (t1) / (t1 * t1);
                    else
                        t2 = 0.0;
                    /*if(t1 < 1.0e-8) t2 = 0.0; */
                    work1[idx] = f[idx] * t2;

                }                   /* end for */
                gcof[ift] = radint1 (work1, r, rab, rg_points);
            }

            break;

        case D_STATE:

            for (int ift = istep * pct.gridpe; ift < istep * pct.gridpe + istep; ift++)
            {
                for (int idx = 0; idx < rg_points; idx++)
                {

                    t1 = r[idx] * gvec[ift];
                    if (t1 > small)
                    {
                        t2 = (THREE / (t1 * t1) - ONE) * sin (t1) / t1;
                        t2 -= THREE * cos (t1) / (t1 * t1);
                    }
                    else
                        t2 = 0.0;
                    /*if(t1 < 1.0e-8) t2 = 0.0; */
                    work1[idx] = f[idx] * t2;

                }                   /* end for */
                gcof[ift] = radint1 (work1, r, rab, rg_points);
            }
            istep = npes * istep;
            GlobalSums (gcof, istep, pct.grid_comm);

            for (int ift = istep; ift < gnum; ift++)
            {
                for (int idx = 0; idx < rg_points; idx++)
                {

                    t1 = r[idx] * gvec[ift];
                    if (t1 > small)
                    {
                        t2 = (THREE / (t1 * t1) - ONE) * sin (t1) / t1;
                        t2 -= THREE * cos (t1) / (t1 * t1);
                    }
                    else
                        t2 = 0.0;
                    /*if(t1 < 1.0e-8) t2 = 0.0; */
                    work1[idx] = f[idx] * t2;

                }                   /* end for */
                gcof[ift] = radint1 (work1, r, rab, rg_points);
            }

            break;

        case F_STATE:

            for (int ift = 0; ift < gnum; ift++)
            {
                for (int idx = 0; idx < rg_points; idx++)
                {
                    t1 = r[idx] * gvec[ift];
                    if (t1 > small)
                    {
                        t2 = (15.0 / (t1 * t1) - 6.0) * sin (t1) / (t1 * t1);
                        t2 += (1.0 - 15.0 / (t1 * t1)) * cos (t1) / t1;
                    }
                    else
                        t2 = 0.0;
                    /*if(t1 < 1.0e-8) t2 = 0.0; */
                    work1[idx] = f[idx] * t2;
                }                   /* end for */
                gcof[ift] = radint1 (work1, r, rab, rg_points);
            }
            break;

        case G_STATE:

            for (int ift = 0; ift < gnum; ift++)
            {
                for (int idx = 0; idx < rg_points; idx++)
                {
                    t1 = r[idx] * gvec[ift];
                    if (t1 > small)
                    {
                        t2 = (105.0 / (t1 * t1 * t1 * t1) - 45.0 / (t1 * t1) + 1.0) * sin (t1) / t1;
                        t2 += (10.0 - 105.0 / (t1 * t1)) * cos (t1) / (t1 * t1);
                    }
                    else
                        t2 = 0.0;
                    /*if(t1 < 1.0e-8) t2 = 0.0; */
                    work1[idx] = f[idx] * t2;
                }                   /*end for */
                gcof[ift] = radint1 (work1, r, rab, rg_points);
            }
            break;

        default:

            throw RmgFatalException() << "angular momentum l=" << lval << "is not programmed in " << __FILE__ << " at line " << __LINE__ << "\n";

    }                       /* end switch */


    /* Fourier Filter the transform and store in work2 */
    for (int idx = 0; idx < gnum; idx++)
    {

        work2[idx] = gcof[idx] * gcutoff (gvec[idx], gcut, width);

    }                           /* end for */


    /* Zero out the array in which the filtered function is to be returned */
    for (int idx = 0; idx < MAX_LOGGRID; idx++)
    {

        ffil[idx] = ZERO;

    }                           /* end for */



    /* Now we reconstruct the filtered function */
    istep = MAX_LOGGRID / npes;

    switch (lval)
    {

    case S_STATE:

        for (int idx = istep * pct.gridpe; idx < istep * pct.gridpe + istep; idx++)
        {


            for (int ift = 0; ift < gnum; ift++)
            {

                t1 = r_filtered[idx] * gvec[ift];
                if (t1 > small)
                    t2 = sin (t1) / t1;
                else
                    t2 = 1.0;
                /*if(t1 < 1.0e-8) t2 = 1.0; */

                /* G-space cutoff */
                work1[ift] = work2[ift] * t2;

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

                t1 = r_filtered[idx] * gvec[ift];
                if (t1 > small)
                    t2 = sin (t1) / t1;
                else
                    t2 = 1.0;
                /*if(t1 < 1.0e-8) t2 = 1.0; */

                /* G-space cutoff */
                work1[ift] = work2[ift] * t2;

            }                   /* end for */


            /* Integrate it */
            t1 = radint (work1, gvec, gnum, gmesh);
            ffil[idx] = t1 * TWO / PI;

        }                       /* end for */




        break;

    case P_STATE:

        for (int idx = istep * pct.gridpe; idx < istep * pct.gridpe + istep; idx++)
        {

            for (int ift = 0; ift < gnum; ift++)
            {

                t1 = r_filtered[idx] * gvec[ift];
                if (t1 > small)
                    t2 = cos (t1) / t1 - sin (t1) / (t1 * t1);
                else
                    t2 = 0.0;
                /*if(t1 < 1.0e-8) t2 = 0.0; */

                /* G-space cutoff */
                work1[ift] = work2[ift] * t2;

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

                t1 = r_filtered[idx] * gvec[ift];
                if (t1 > small)
                    t2 = cos (t1) / t1 - sin (t1) / (t1 * t1);
                else
                    t2 = 0.0;
                /*if(t1 < 1.0e-8) t2 = 0.0; */

                /* G-space cutoff */
                work1[ift] = work2[ift] * t2;

            }                   /* end for */


            /* Integrate it */
            t1 = radint (work1, gvec, gnum, gmesh);
            ffil[idx] = t1 * TWO / PI;

        }                       /* end for */

        break;

    case D_STATE:

        for (int idx = istep * pct.gridpe; idx < istep * pct.gridpe + istep; idx++)
        {

            for (int ift = 0; ift < gnum; ift++)
            {

                t1 = r_filtered[idx] * gvec[ift];
                if (t1 > small)
                {
                    t2 = (THREE / (t1 * t1) - ONE) * sin (t1) / t1;
                    t2 -= THREE * cos (t1) / (t1 * t1);
                }
                else
                    t2 = 0.0;
                /*if(t1 < 1.0e-8) t2 = 0.0; */

                /* G-space cutoff */
                work1[ift] = work2[ift] * t2;

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

                t1 = r_filtered[idx] * gvec[ift];
                if (t1 > small)
                {
                    t2 = (THREE / (t1 * t1) - ONE) * sin (t1) / t1;
                    t2 -= THREE * cos (t1) / (t1 * t1);
                }
                else
                    t2 = 0.0;
                /*if(t1 < 1.0e-8) t2 = 0.0; */

                /* G-space cutoff */
                work1[ift] = work2[ift] * t2;

            }                   /* end for */

            /* Integrate it */
            t1 = radint (work1, gvec, gnum, gmesh);
            ffil[idx] = t1 * TWO / PI;

        }                       /* end for */

        break;

    case F_STATE:

        for (int idx = istep * pct.gridpe; idx < istep * pct.gridpe + istep; idx++)
        {

            for (int ift = 0; ift < gnum; ift++)
            {

                t1 = r_filtered[idx] * gvec[ift];
                if (t1 > small)
                {
                    t2 = (15.0 / (t1 * t1) - 6.0) * sin (t1) / (t1 * t1);
                    t2 += (1.0 - 15.0 / (t1 * t1)) * cos (t1) / t1;
                }
                else
                    t2 = 0.0;
                /*if(t1 < 1.0e-8) t2 = 0.0; */

                /* G-space cutoff */
                work1[ift] = work2[ift] * t2;

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

                t1 = r_filtered[idx] * gvec[ift];
                if (t1 > small)
                {
                    t2 = (15.0 / (t1 * t1) - 6.0) * sin (t1) / (t1 * t1);
                    t2 += (1.0 - 15.0 / (t1 * t1)) * cos (t1) / t1;
                }
                else
                    t2 = 0.0;
                /*if(t1 < 1.0e-8) t2 = 0.0; */

                /* G-space cutoff */
                work1[ift] = work2[ift] * t2;

            }                   /* end for */

            /* Integrate it */
            t1 = radint (work1, gvec, gnum, gmesh);
            ffil[idx] = t1 * TWO / PI;

        }                       /* end for */

        break;

    case G_STATE:

        for (int idx = istep * pct.gridpe; idx < istep * pct.gridpe + istep; idx++)
        {

            for (int ift = 0; ift < gnum; ift++)
            {

                t1 = r_filtered[idx] * gvec[ift];
                if (t1 > small)
                {
                    t2 = (105.0 / (t1 * t1 * t1 * t1) - 45.0 / (t1 * t1) + 1.0) * sin (t1) / t1;
                    t2 += (10.0 - 105.0 / (t1 * t1)) * cos (t1) / (t1 * t1);
                }
                else
                    t2 = 0.0;
                /*if(t1 < 1.0e-8) t2 = 0.0; */

                /* G-space cutoff */
                work1[ift] = work2[ift] * t2;

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

                t1 = r_filtered[idx] * gvec[ift];
                if (t1 > small)
                {
                    t2 = (105.0 / (t1 * t1 * t1 * t1) - 45.0 / (t1 * t1) + 1.0) * sin (t1) / t1;
                    t2 += (10.0 - 105.0 / (t1 * t1)) * cos (t1) / (t1 * t1);
                }
                else
                    t2 = 0.0;
                /*if(t1 < 1.0e-8) t2 = 0.0; */

                /* G-space cutoff */
                work1[ift] = work2[ift] * t2;

            }                   /* end for */

            /* Integrate it */
            t1 = radint (work1, gvec, gnum, gmesh);
            ffil[idx] = t1 * TWO / PI;

        }                       /* end for */

        break;

    default:

        throw RmgFatalException() << "angular momentum l=" << lval << "is not programmed in " << __FILE__ << " at line " << __LINE__ << "\n";

    }                           /* end switch */


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




void Atomic::FilterPotential (
    double *potential,       // IN:  potential to be filtered
    double *r,               // IN:  radial grid for the potential
    int rg_points,           // IN:  number of points in the radial grid
    double rmax,             // IN:  mask function parameter
    double offset,           // IN:  mask function parameter
    double parm,             // IN:  filtering parameter
    double* potential_lgrid, // OUT: potential on standard log grid
    double *rab,             // IN:  radial volume element from the pseudopotential file
    int l_value,             // IN:  momentum l (0=s, 1=p, 2=d)
    double gwidth,           // IN:  gaussian parameter controlling speed of g-space damping
    double rcut,             // IN:  filtered potential is damped in real space starting here
    double rwidth,           // IN:  gaussian parameter controlling speed of real space damping
    double *drpotential_lgrid) // OUT: derivative of potential on standard log grid
{

    bool der_flag = false;
    const double small = 1.e-8;


    if(drpotential_lgrid)
	der_flag = true;

//    if (ct.mask_function)
//	apply_mask_function(potential, r, rg_points, rmax, offset);


    /* Transform to g-space and filter it */
    RftToLogGrid (parm, potential, r, potential_lgrid, rab, rg_points, l_value, gwidth);

//    if (ct.mask_function)
//	backout_mask_function (potential_lgrid, dr, MAX_LOGGRID,  rmax);

    /*Evaluate radial derivative, if requested*/
    if (der_flag)
    {

	radiff (potential_lgrid, drpotential_lgrid, r_filtered, MAX_LOGGRID, log(LOGGRID_MESHPARM));

	/* Fix up the first point */
	drpotential_lgrid[1] = 2.0 * drpotential_lgrid[2] - drpotential_lgrid[3];
	drpotential_lgrid[0] = 2.0 * drpotential_lgrid[1] - drpotential_lgrid[2];

    }

    /*Fix up first point in filtered potential*/
    potential_lgrid[0] = 2.0 * potential_lgrid[1] - potential_lgrid[2];


    /* Without mask function, result needs to be dampened by gaussian starting at rcut*/
    if (!ct.mask_function)
    {
	for (int idx = 0; idx < MAX_LOGGRID; idx++)
	{
            double rdist = r_filtered[idx];

	    if (rdist > rcut)
	    {
		double t1 = (rdist - rcut) / rcut;
		double exp_fac = exp (-rwidth * t1 * t1);

		potential_lgrid[idx] *= exp_fac;

		if (fabs (potential_lgrid[idx]) < small)
		    potential_lgrid[idx] = 0.0;

		if (der_flag)
		{
		    drpotential_lgrid[idx] *= exp_fac;

		    if (fabs (drpotential_lgrid[idx]) < small)
			drpotential_lgrid[idx] = 0.0;
		}

	    }               /* end if */

	}
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

