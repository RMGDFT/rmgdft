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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "transition.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "RmgException.h"
#include "Atomic.h"


/*For a quantity localized around ionic positions, this function finds radius in number of grid points
 * given a radius in a.u.*/
int Radius2grid (double radius, double mingrid_spacing)
{

    double scale, t1, t2;
    int it1, dim, ibrav;

    ibrav = Rmg_L.get_ibrav_type();
    
    /* Set the scaling factor for determining the radius of the local grids */
    scale = 1.0;
    if (ibrav == CUBIC_BC)
	scale = 1.1;
    if (ibrav == CUBIC_FC)
	scale = 1.3;
        
	t1 = 2.0 * scale * radius / mingrid_spacing;
        t1 = modf (t1, &t2);
        it1 = (int) t2;
        if (t1 > 0.5)
            it1++;
        if (!(it1 % 2))
            it1++;
	dim = it1;

	return dim;
}


void InitPseudo (std::unordered_map<std::string, InputKey *>& ControlMap)
{

    int isp, ip;
    double Zv, rc, t1;
    char newname[MAX_PATH];
    FILE *psp = NULL;
    FILE *psp2 = NULL;
    Atomic *A = new Atomic();
    double *rgrid = A->GetRgrid();
    int ibrav = Rmg_L.get_ibrav_type();

    double *work = new double[MAX_LOGGRID];

    bool write_flag = Verify ("write_pseudopotential_plots", true, ControlMap);

    /*Initialize max_nlpoints and max_nlfpoints */
    ct.max_nlpoints = 0;
    ct.max_nlfpoints = 0;


    /* Loop over species */
    for (isp = 0; isp < ct.num_species; isp++)
    {
        SPECIES *sp = &ct.sp[isp];
        if (pct.gridpe == 0 && write_flag)
        {
            snprintf (newname, MAX_PATH, "local_%s.xmgr",  sp->atomic_symbol);
            if(NULL == (psp = fopen (newname, "w+")))
               throw RmgFatalException() << "Unable to open pseudopotential graph file " << " in " << __FILE__ << " at line " << __LINE__ << "\n";
 
        }


        /*Get nldim */
        int done = false;
        while(!done) {

            sp->nldim = Radius2grid (sp->nlradius, ct.hmingrid);
            sp->nldim = sp->nldim/2*2 +1;
            sp->nlfdim = ct.nxfgrid * sp->nldim;

            if ((sp->nldim >= get_NX_GRID()) || (sp->nldim >= get_NY_GRID()) || (sp->nldim >= get_NZ_GRID())) {
                printf("Warning: diameter of non-local projectors %8.4f exceeds cell size. Reducing.\n", sp->nlradius);
                sp->nlradius *= 0.95;
            }
            else {
                done = true;
            }

        }


        /*ct.max_nlpoints is max of nldim*nldim*nldim for all species */
        if (ct.max_nlpoints < (sp->nldim * sp->nldim * sp->nldim))
            ct.max_nlpoints = sp->nldim * sp->nldim * sp->nldim;

        if (ct.max_nlfpoints < (sp->nlfdim * sp->nlfdim * sp->nlfdim))
            ct.max_nlfpoints = sp->nlfdim * sp->nlfdim * sp->nlfdim;


        /*Filter and interpolate local potential into the internal log grid*/
        Zv = sp->zvalence;
        rc = sp->rc;

        /* Generate the difference potential */
        for (int idx = 0; idx < sp->rg_points; idx++)
            work[idx] = sp->vloc0[idx] + Zv * erf (sp->r[idx] / rc) / sp->r[idx];


        if (pct.gridpe == 0 && write_flag)
        {
            for (int idx = 0; idx < sp->rg_points; idx++)
                fprintf (psp, "%e  %e\n", sp->r[idx], work[idx]);

            fprintf (psp, "\n&&\n");
        }



        // Transform to g-space and filter it with filtered function returned on standard log grid
        A->FilterPotential(work, sp->r, sp->rg_points, sp->lradius, 0.25, ct.cparm, sp->localig,
                           sp->rab, 0, sp->gwidth, sp->lrcut, sp->rwidth, sp->drlocalig);

        /*Write local projector into a file if requested*/
        if ((pct.gridpe == 0) && write_flag)
        {
            for (int idx = 0; idx < MAX_LOGGRID; idx++)
                fprintf (psp, "%e  %e \n", rgrid[idx], sp->localig[idx]);

            /* output xmgr data separator */
            fprintf (psp, "\n&&\n");

            for (int idx = 0; idx < MAX_LOGGRID; idx++)
                fprintf (psp, "%e  %e \n", rgrid[idx], sp->drlocalig[idx]);

            fclose (psp);
        }


        /*Open file for writing beta function*/
        if (pct.gridpe == 0 && write_flag)
        {
            snprintf (newname, MAX_PATH, "beta_%s.xmgr", sp->atomic_symbol);
            if(NULL == (psp = fopen (newname, "w+")))
               throw RmgFatalException() << "Unable to open beta function graph file " << " in " << __FILE__ << " at line " << __LINE__ << "\n";

            snprintf (newname, MAX_PATH, "drbeta_%s.xmgr", sp->atomic_symbol);
            if(NULL == (psp2 = fopen (newname, "w+")))
               throw RmgFatalException() << "Unable to open drbeta function graph file " << " in " << __FILE__ << " at line " << __LINE__ << "\n";
        }

        /* Write raw beta function into file if requested*/
        for (ip = 0; ip < sp->nbeta; ip++)
        {

            if (pct.gridpe == 0 && write_flag)
            {
                for (int idx = 0; idx < sp->kkbeta; idx++)
                    fprintf (psp, "%e  %e\n", sp->r[idx], sp->beta[ip][idx]);
                fprintf (psp, "\n&&\n");
            }

            A->FilterPotential(&sp->beta[ip][0], sp->r, sp->rg_points, sp->nlradius, 0.25, ct.betacparm, &sp->betalig[ip][0],
            sp->rab, sp->llbeta[ip], sp->gwidth, sp->nlrcut[sp->llbeta[ip]], sp->rwidth, &sp->drbetalig[ip][0]);


            /* Is this necessary ??? */
            if (sp->llbeta[ip])
                sp->betalig[ip][0] = 0.0;


            /* output filtered non-local projector to a file  if requested */
            if (pct.gridpe == 0 && write_flag)
            {
                for (int idx = 0; idx < MAX_LOGGRID; idx++)
                {
                    {
                        fprintf (psp, "%e  %e\n", rgrid[idx], sp->betalig[ip][idx]);
                        fprintf (psp2, "%e  %e\n", rgrid[idx], sp->drbetalig[ip][idx]);
                    }
                }                   /* end for */
            }

            /* output xmgr data separator */
            if (pct.gridpe == 0 && write_flag)
            {
                fprintf (psp, "\n&&\n");
                fprintf (psp2, "\n&&\n");
            }

        }                       /* end for ip */

        /* Now take care of the core charge if nlcc is being used */
        if (sp->nlccflag)
        {

            for (int idx = 0; idx < sp->rg_points; idx++)
                work[idx] = sp->rspsco[idx] / (4.0 * PI);


            /*Write raw (pre-filtered) data to a file if requested */
            if (pct.gridpe == 0 && write_flag)
            {
                for (int idx = 0; idx < sp->rg_points; idx++)
                    fprintf (psp, "%e  %e\n", sp->r[idx], work[idx]);
                fprintf (psp, "\n&&\n");
            }

            A->FilterPotential(work, sp->r, sp->rg_points, sp->lradius, 0.25, ct.cparm, &sp->rhocorelig[0],
                           sp->rab, 0, sp->gwidth, sp->lrcut, sp->rwidth, NULL);

            /*Oscilations at the tail end of filtered function may cause rhocore to be negative
             * but I am not sure if this is the right solution, it may be better to fix charge density
             * this rhocore less smooth*/
            for (int idx = 0; idx < MAX_LOGGRID; idx++)
                if (sp->rhocorelig[idx] < 0.0)
                    sp->rhocorelig[idx] = 0.0;



            /*Write filtered data to a file if requested */
            if (pct.gridpe == 0 && write_flag)
            {
                for (int idx = 0; idx < MAX_LOGGRID; idx++)
                    fprintf (psp, "%e  %e\n", rgrid[idx], sp->rhocorelig[idx]);
                fprintf (psp, "\n&&\n");
            }

        }                       /* end if */


        if (pct.gridpe == 0 && write_flag)
        {
            fclose (psp);
            fclose (psp2);
        }

    }                           /* end for */


    delete A;
    delete [] work;

} // end InitPseudo
