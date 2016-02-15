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
        Log_grid_initialized = true;

    }

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


    /* Get some temporary memory */
    int alloc = rg_points;
    if(alloc < RADIAL_GVECS) alloc = RADIAL_GVECS;
    if (alloc < MAX_LOGGRID)
        alloc = MAX_LOGGRID;

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



    /* Loop over frequency components */
    for (int ift = 0; ift < gnum; ift++)
    {

        switch (lval)
        {

        case S_STATE:

            for (int idx = 0; idx < rg_points; idx++)
            {

                t1 = r[idx] * gvec[ift];

                if (t1 > SMALL)
                    work1[idx] = f[idx] * sin (t1) / t1;
                else
                    work1[idx] = f[idx];
                /*t2 = sin(t1) / t1;
                   if(t1 < 1.0e-8) t2 = 1.0;
                   work1[idx] = f[idx] * t2; */

            }                   /* end for */

            break;

        case P_STATE:

            for (int idx = 0; idx < rg_points; idx++)
            {

                t1 = r[idx] * gvec[ift];
                if (t1 > SMALL)
                    t2 = cos (t1) / t1 - sin (t1) / (t1 * t1);
                else
                    t2 = 0.0;
                /*if(t1 < 1.0e-8) t2 = 0.0; */
                work1[idx] = f[idx] * t2;

            }                   /* end for */

            break;

        case D_STATE:

            for (int idx = 0; idx < rg_points; idx++)
            {

                t1 = r[idx] * gvec[ift];
                if (t1 > SMALL)
                {
                    t2 = (THREE / (t1 * t1) - ONE) * sin (t1) / t1;
                    t2 -= THREE * cos (t1) / (t1 * t1);
                }
                else
                    t2 = 0.0;
                /*if(t1 < 1.0e-8) t2 = 0.0; */
                work1[idx] = f[idx] * t2;

            }                   /* end for */

            break;

        case F_STATE:

            for (int idx = 0; idx < rg_points; idx++)
            {
                t1 = r[idx] * gvec[ift];
                if (t1 > SMALL)
                {
                    t2 = (15.0 / (t1 * t1) - 6.0) * sin (t1) / (t1 * t1);
                    t2 += (1.0 - 15.0 / (t1 * t1)) * cos (t1) / t1;
                }
                else
                    t2 = 0.0;
                /*if(t1 < 1.0e-8) t2 = 0.0; */
                work1[idx] = f[idx] * t2;
            }                   /* end for */

            break;

        case G_STATE:

            for (int idx = 0; idx < rg_points; idx++)
            {
                t1 = r[idx] * gvec[ift];
                if (t1 > SMALL)
                {
                    t2 = (105.0 / (t1 * t1 * t1 * t1) - 45.0 / (t1 * t1) + 1.0) * sin (t1) / t1;
                    t2 += (10.0 - 105.0 / (t1 * t1)) * cos (t1) / (t1 * t1);
                }
                else
                    t2 = 0.0;
                /*if(t1 < 1.0e-8) t2 = 0.0; */
                work1[idx] = f[idx] * t2;
            }                   /*end for */

            break;

        default:

            throw RmgFatalException() << "angular momentum l=" << lval << "is not programmed in " << __FILE__ << " at line " << __LINE__ << "\n";

        }                       /* end switch */


        /* Get coefficients */
        gcof[ift] = radint1 (work1, r, rab, rg_points);

    }                           /* end for */

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
                if (t1 > SMALL)
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
                if (t1 > SMALL)
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
                if (t1 > SMALL)
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
                if (t1 > SMALL)
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
                if (t1 > SMALL)
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
                if (t1 > SMALL)
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
                if (t1 > SMALL)
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
                if (t1 > SMALL)
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
                if (t1 > SMALL)
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
                if (t1 > SMALL)
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

	radiff (potential_lgrid, drpotential_lgrid, r_filtered, MAX_LOGGRID, LOGGRID_MESHPARM);

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

		if (fabs (potential_lgrid[idx]) < SMALL)
		    potential_lgrid[idx] = 0.0;

		if (der_flag)
		{
		    drpotential_lgrid[idx] *= exp_fac;

		    if (fabs (drpotential_lgrid[idx]) < SMALL)
			drpotential_lgrid[idx] = 0.0;
		}

	    }               /* end if */

	}
    }



} 

