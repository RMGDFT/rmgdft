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


void InitPseudo (std::unordered_map<std::string, InputKey *>& ControlMap)
{

    int isp, ip;
    double Zv, rc;
    char newname[MAX_PATH];
    FILE *psp = NULL;
    Atomic *A = new Atomic();
    double *rgrid = A->GetRgrid();

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


        // Might need to adjust this depending on filtering changes. Also assumes that all beta have roughly the same range
        sp->nlradius = 4.5 * A->GetRange(&sp->beta[0][0], sp->r, sp->rab, sp->rg_points);

        /*Get nldim */
        bool done = false;
        bool reduced = false;
        
        while(!done) {

            sp->nldim = Radius2grid (sp->nlradius, ct.hmingrid);
            sp->nldim = sp->nldim/4*4;
            sp->nlfdim = ct.nxfgrid * sp->nldim;

            if ((sp->nldim >= get_NX_GRID()) || (sp->nldim >= get_NY_GRID()) || (sp->nldim >= get_NZ_GRID())) {
                sp->nlradius *= 0.99;
                reduced = true;
            }
            else {
                done = true;
            }

        }
        sp->nlradius = 0.5 * ct.hmingrid * (double)(sp->nldim);
        sp->nlradius -= 0.5 * ct.hmingrid / (double)ct.nxfgrid;
        if(reduced && ct.localize_projectors) rmg_printf("Warning: diameter of non-local projectors exceeds cell size. Reducing. New radius = %12.6f\n", sp->nlradius);

        // If projectors will span the full wavefunction grid then use a larger value for the nlradius for all remaining operations
        if(!ct.localize_projectors) {
            sp->nlradius = 10.0;
            sp->nldim = Radius2grid (sp->nlradius, ct.hmingrid);
            sp->nldim = sp->nldim/4*4;
            sp->nlfdim = ct.nxfgrid * sp->nldim;
        }

        /*ct.max_nlpoints is max of nldim*nldim*nldim for all species */
        if (ct.max_nlpoints < (sp->nldim * sp->nldim * sp->nldim))
            ct.max_nlpoints = sp->nldim * sp->nldim * sp->nldim;

        if (ct.max_nlfpoints < (sp->nlfdim * sp->nlfdim * sp->nlfdim))
            ct.max_nlfpoints = sp->nlfdim * sp->nlfdim * sp->nlfdim;

        if(!ct.localize_projectors) {
            ct.max_nlpoints = get_NX_GRID() * get_NY_GRID() * get_NZ_GRID();
            ct.max_nlfpoints = ct.max_nlpoints * ct.nxfgrid * ct.nxfgrid * ct.nxfgrid;
        }

        sp->ldim = Radius2grid (sp->lradius, ct.hmingrid / (double)Rmg_G->default_FG_RATIO);
        sp->ldim = sp->ldim/2 * 2 +1;
        sp->lradius = 0.5 * ct.hmingrid * (double)(sp->ldim-1) / (double)Rmg_G->default_FG_RATIO;
        sp->lradius -= 0.5 * ct.hmingrid / (double)Rmg_G->default_FG_RATIO;
        sp->ldim_fine = sp->ldim *ct.nxfgrid/Rmg_G->default_FG_RATIO;
        if(ct.nxfgrid %Rmg_G->default_FG_RATIO != 0)
        {
            sp->ldim_fine = sp->ldim *2;
            printf("\n WARNING::   FG_RATIO %d is not a factor of ct.nxfgrid %d\n",
                    Rmg_G->default_FG_RATIO , ct.nxfgrid);
            fflush(NULL);
        }



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
        double parm = ct.rhocparm / ct.FG_RATIO;
        A->FilterPotential(work, sp->r, sp->rg_points, sp->lradius, parm, sp->localig,
                sp->rab, 0, sp->gwidth, 0.66*sp->lradius, sp->rwidth, ct.hmingrid / (double)Rmg_G->default_FG_RATIO);

        /*Write local projector into a file if requested*/
        if ((pct.gridpe == 0) && write_flag)
        {
            for (int idx = 0; idx < MAX_LOGGRID; idx++)
                fprintf (psp, "%e  %e \n", rgrid[idx], sp->localig[idx]);

            /* output xmgr data separator */
            fprintf (psp, "\n&&\n");
            fclose (psp);
        }

        // Next we want to fourier filter the input atomic charge density and transfer
        // it to the interpolation grid for use by LCAO starts and force correction routines
        A->FilterPotential(sp->atomic_rho, sp->r, sp->rg_points, sp->lradius, ct.rhocparm, sp->arho_lig,
                sp->rab, 0, sp->gwidth, 0.66*sp->lradius, sp->rwidth, ct.hmingrid);


        /*Open file for writing beta function*/
        if (pct.gridpe == 0 && write_flag)
        {
            snprintf (newname, MAX_PATH, "beta_%s.xmgr", sp->atomic_symbol);
            if(NULL == (psp = fopen (newname, "w+")))
                throw RmgFatalException() << "Unable to open beta function graph file " << " in " << __FILE__ << " at line " << __LINE__ << "\n";
        }

        /* Write raw beta function into file if requested*/
        for (ip = 0; ip < sp->nbeta; ip++)
        {

            if (pct.gridpe == 0 && write_flag)
            {
                for (int idx = 0; idx < sp->kkbeta; idx++) fprintf (psp, "%e  %e\n", sp->r[idx], sp->beta[ip][idx]);
                fprintf (psp, "\n&&\n");
            }

            sp->beta_g[ip] = new double[RADIAL_GVECS];
            // Filtering for wavefunction grid
            A->FilterPotential(&sp->beta[ip][0], sp->r, sp->rg_points, sp->nlradius, ct.betacparm, &sp->betalig[ip][0],
                    sp->rab, sp->llbeta[ip], sp->gwidth, 0.66*sp->nlradius, sp->rwidth, ct.hmingrid);
            A->RLogGridToGLogGrid(ct.betacparm, &sp->beta[ip][0], sp->r, sp->rab, sp->beta_g[ip],
                    sp->rg_points, sp->llbeta[ip], sp->gwidth);

            /* output filtered non-local projector to a file  if requested */
            if (pct.gridpe == 0 && write_flag)
            {
                for (int idx = 0; idx < MAX_LOGGRID; idx++)
                {
                    {
                        fprintf (psp, "%e  %e\n", rgrid[idx], sp->betalig[ip][idx]);
                    }
                }                   /* end for */
            }

            /* output xmgr data separator */
            if (pct.gridpe == 0 && write_flag)
            {
                fprintf (psp, "\n&&\n");
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

            double nlccradius = 3.5 * A->GetRange(work, sp->r, sp->rab, sp->rg_points);

            // Make adjustments so radii terminates on a grid point
            int nlccdim = Radius2grid (nlccradius, ct.hmingrid/(double)Rmg_G->default_FG_RATIO);
            nlccdim = nlccdim/2*2 + 1;
            nlccradius = 0.5 * ct.hmingrid * (double)(nlccdim-1) / (double)Rmg_G->default_FG_RATIO;
            nlccradius -= 0.5 * ct.hmingrid / (double)Rmg_G->default_FG_RATIO;
            double nlcccut = 0.66 * nlccradius;

            A->FilterPotential(work, sp->r, sp->rg_points, nlccradius, ct.rhocparm, &sp->rhocorelig[0],
                           sp->rab, 0, sp->gwidth, nlcccut, sp->rwidth, ct.hmingrid/(double)Rmg_G->default_FG_RATIO);

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
        }

    }                           /* end for */


    delete A;
    delete [] work;

} // end InitPseudo
