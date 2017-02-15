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



#include "portability.h"
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
#include <boost/math/special_functions/erf.hpp>
#include "Atomic.h"

using boost::math::policies::policy;
using boost::math::policies::promote_double;
typedef policy<promote_double<false> > bessel_policy;

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
        sp->nlradius = ct.projector_expansion_factor * 4.5 * A->GetRange(&sp->beta[0][0], sp->r, sp->rab, sp->rg_points);

        /*Get nldim */
        bool done = false;
        bool reduced = false;
        
        while(!done) {

            sp->nldim = Radius2grid (sp->nlradius, ct.hmingrid);
            sp->nldim = sp->nldim/2*2 + 1;
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
            sp->nlradius = 7.0;
            sp->nldim = Radius2grid (sp->nlradius, ct.hmingrid);
            sp->nldim = sp->nldim/2*2 +1;
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
        sp->ldim_fine = sp->ldim ;



        /*Filter and interpolate local potential into the internal log grid*/
        Zv = sp->zvalence;
        rc = sp->rc;

        /* Generate the difference potential */
        for (int idx = 0; idx < sp->rg_points; idx++)
            work[idx] = sp->vloc0[idx] + Zv * boost::math::erf (sp->r[idx] / rc) / sp->r[idx];

        /* Write raw beta function into file if requested*/
        double *bessel_rg = new double[(ct.max_l+1) * RADIAL_GVECS * sp->rg_points];
            RmgTimer *RT1 = new RmgTimer("radial beta");
        A->InitBessel(sp->r, sp->rg_points, ct.max_l, bessel_rg);
            delete RT1;
        // get the local pseudopotential in G space
        sp->localpp_g = new double[RADIAL_GVECS];
        sp->arho_g = new double[RADIAL_GVECS];
        sp->rhocore_g = new double[RADIAL_GVECS];
        A->RLogGridToGLogGrid(work, sp->r, sp->rab, sp->localpp_g,
                sp->rg_points, 0, bessel_rg);
        A->RLogGridToGLogGrid(sp->atomic_rho, sp->r, sp->rab, sp->arho_g,
                sp->rg_points, 0, bessel_rg);

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

#if 1
        // Get the G=0 component
        double fac = Zv * 4.0 * PI;
        for (int idx = 0; idx < sp->rg_points; idx++)
            work[idx] = sp->vloc0[idx] + Zv/sp->r[idx];
        double G0 = 4.0*PI*radint1 (work, sp->r, sp->rab, sp->rg_points);
        sp->localpp_g[0] = G0;

        // Subtract off analytic transform of erf that was added above
        for(int idx=1;idx < RADIAL_GVECS;idx++)
        {
            double gval = A->gvec[idx]*A->gvec[idx];
            sp->localpp_g[idx] -= fac * exp ( - 0.25*gval*rc*rc) / gval;
        }
#endif
        // Next we want to fourier filter the input atomic charge density and transfer
        // it to the interpolation grid for use by LCAO starts and force correction routines
        
        for (int idx = 0; idx < sp->rg_points; idx++)
        if(sp->atomic_rho[idx] < 0.0) work[idx] = 0.0;
        else work[idx] = sqrt(sp->atomic_rho[idx]);
        A->FilterPotential(work, sp->r, sp->rg_points, sp->lradius, ct.rhocparm, sp->arho_lig,
                sp->rab, 0, sp->gwidth, 0.66*sp->lradius, sp->rwidth, ct.hmingrid);
       for (int idx = 0; idx < MAX_LOGGRID; idx++)
           sp->arho_lig[idx] *= sp->arho_lig[idx];


        /*Open file for writing beta function*/
        if (pct.gridpe == 0 && write_flag)
        {
            snprintf (newname, MAX_PATH, "beta_%s.xmgr", sp->atomic_symbol);
            if(NULL == (psp = fopen (newname, "w+")))
                throw RmgFatalException() << "Unable to open beta function graph file " << " in " << __FILE__ << " at line " << __LINE__ << "\n";
        }

        for (ip = 0; ip < sp->nbeta; ip++)
        {

            if (pct.gridpe == 0 && write_flag)
            {
                for (int idx = 0; idx < sp->kkbeta; idx++) fprintf (psp, "%e  %e\n", sp->r[idx], sp->beta[ip][idx]);
                fprintf (psp, "\n&&\n");
            }

            sp->beta_g[ip] = new double[RADIAL_GVECS];
            A->RLogGridToGLogGrid(&sp->beta[ip][0], sp->r, sp->rab, sp->beta_g[ip],
                    sp->rg_points, sp->llbeta[ip], bessel_rg);


        }                       /* end for ip */



        /* Now take care of the core charge if nlcc is being used */
        if (sp->nlccflag)
        {

            for (int idx = 0; idx < sp->rg_points; idx++)
                work[idx] = sp->rspsco[idx] / (4.0 * PI);

            A->RLogGridToGLogGrid(work, sp->r, sp->rab, sp->rhocore_g,
                sp->rg_points, 0, bessel_rg);

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

        delete [] bessel_rg;

    }                           /* end for */


    delete A;
    delete [] work;

} // end InitPseudo
