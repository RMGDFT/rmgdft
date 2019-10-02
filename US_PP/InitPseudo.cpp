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
#include <boost/math/special_functions/erf.hpp>
#include "Atomic.h"

using boost::math::policies::policy;
using boost::math::policies::promote_double;
typedef policy<promote_double<false> > bessel_policy;

void InitPseudo (std::unordered_map<std::string, InputKey *>& ControlMap)
{

    int isp;
    double Zv, rc;
    char newname[MAX_PATH];
    FILE *psp = NULL;
    Atomic *A = new Atomic();
    double *rgrid = A->GetRgrid();


    bool write_flag = Verify ("write_pseudopotential_plots", true, ControlMap);

    /*Initialize max_nlpoints */
    ct.max_nlpoints = 0;


    /* Loop over species */
    for (isp = 0; isp < ct.num_species; isp++)
    {
        SPECIES *sp = &Species[isp];
        double *work = new double[std::max(MAX_LOGGRID, sp->rg_points)];
        if(!std::strcmp(sp->atomic_symbol, "DLO")) continue;
        if (pct.gridpe == 0 && write_flag)
        {
            snprintf (newname, MAX_PATH, "local_%s.xmgr",  sp->atomic_symbol);
            if(NULL == (psp = fopen (newname, "w+")))
               throw RmgFatalException() << "Unable to open pseudopotential graph file " << " in " << __FILE__ << " at line " << __LINE__ << "\n";
 
        }


        // Might need to adjust this depending on filtering changes. Also assumes that all beta have roughly the same range
        sp->nlradius = ct.projector_expansion_factor * 4.5 * A->GetRange(&sp->beta[0][0], sp->r, sp->rab, sp->rg_points, 0.999999999);
        sp->nlradius = std::max(sp->nlradius, ct.min_nlradius);
        sp->nlradius = std::min(sp->nlradius, ct.max_nlradius);


        /*Get nldim */
        sp->nldim = Radius2grid (sp->nlradius, ct.hmingrid, Rmg_L.get_ibrav_type(), ct.localize_projectors);
        sp->nldim = sp->nldim/2*2 + 1;

        int printed = 0;
        while ((sp->nldim >= get_NX_GRID()) || (sp->nldim >= get_NY_GRID()) || (sp->nldim >= get_NZ_GRID()))
        {
            sp->nlradius -= ct.hmingrid;
            sp->nldim = Radius2grid (sp->nlradius, ct.hmingrid, Rmg_L.get_ibrav_type(), ct.localize_projectors);
            sp->nldim = sp->nldim/2*2 + 1;
            if(ct.localize_projectors)
            {
                if(ct.rmg_branch == RMG_BASE)
                {
                    if((pct.gridpe == 0) && (printed == 0))
                    {
                        rmg_printf("Warning: localized projectors selected but their diameter exceeds cell size. Switching to delocalized.\n");
                        printf("Warning: localized projectors selected but their diameter exceeds cell size. Switching to delocalized.\n");
                        printed++;
                    }
                    ct.localize_projectors = false;
                }
                else if(printed == 0)
                {
                    rmg_printf("Warning: localized projectors selected but their diameter exceeds cell size.\n");
                    printf("Warning: localized projectors selected but their diameter exceeds cell size.\n");
                    printed++;
                }
            }
        }

        sp->nlradius = 0.5 * ct.hmingrid * (double)(sp->nldim);

        // If projectors will span the full wavefunction grid then use a larger value for the nlradius for all remaining operations
        if(!ct.localize_projectors) {
            sp->nlradius = std::max(sp->nlradius, 7.0);
            sp->nldim = Radius2grid (sp->nlradius, ct.hmingrid, Rmg_L.get_ibrav_type(), ct.localize_projectors);
            sp->nldim = sp->nldim/2*2 +1;
        }

        /*ct.max_nlpoints is max of nldim*nldim*nldim for all species */
        if (ct.max_nlpoints < (sp->nldim * sp->nldim * sp->nldim))
            ct.max_nlpoints = sp->nldim * sp->nldim * sp->nldim;

        if(!ct.localize_projectors) {
            ct.max_nlpoints = get_NX_GRID() * get_NY_GRID() * get_NZ_GRID();
        }

        if(ct.localize_projectors)
        {
            // Grid object local to this MPI process
            sp->OG = (void *)new BaseGrid(sp->nldim, sp->nldim, sp->nldim, 1, 1, 1, 0, 1);
            BaseGrid *OG = (BaseGrid *)sp->OG;
            OG->set_rank(0, pct.my_comm);

            // Lattice object for the localized projector region. Global vectors need to be
            // scaled so that they correspond to the local region.
            Lattice *L = new Lattice();
            double a0[3], a1[3], a2[3], s1, celldm[6], omega;
            s1 = (double)sp->nldim / (double)Rmg_G->get_NX_GRID(1);
            a0[0] = s1*Rmg_L.get_a0(0);a0[1] = s1*Rmg_L.get_a0(1);a0[2] = s1*Rmg_L.get_a0(2);
            s1 = (double)sp->nldim / (double)Rmg_G->get_NY_GRID(1);
            a1[0] = s1*Rmg_L.get_a1(0);a1[1] = s1*Rmg_L.get_a1(1);a1[2] = s1*Rmg_L.get_a1(2);
            s1 = (double)sp->nldim / (double)Rmg_G->get_NZ_GRID(1);
            a2[0] = s1*Rmg_L.get_a2(0);a2[1] = s1*Rmg_L.get_a2(1);a2[2] = s1*Rmg_L.get_a2(2);

            L->set_ibrav_type(None);
            L->latgen(celldm, &omega, a0, a1, a2);

            sp->prj_pwave = new Pw(*OG, *L, 1, false);
        }
        else
        {
            sp->prj_pwave = new Pw(*Rmg_G, Rmg_L, 1, false);
        }


        sp->ldim = Radius2grid (sp->lradius, ct.hmingrid / (double)Rmg_G->default_FG_RATIO, Rmg_L.get_ibrav_type(), ct.localize_localpp);
        sp->ldim = sp->ldim/2 * 2 +1;
        if(!ct.localize_localpp)
        {
            sp->lradius = 10.0;
            sp->ldim = Radius2grid (sp->lradius, ct.hmingrid / (double)Rmg_G->default_FG_RATIO, Rmg_L.get_ibrav_type(), ct.localize_localpp);
            sp->ldim = sp->nldim/2*2 + 1;
        }
        else
        {
            sp->lradius = 0.5 * ct.hmingrid * (double)(sp->ldim-1) / (double)Rmg_G->default_FG_RATIO;
            sp->lradius -= 0.5 * ct.hmingrid / (double)Rmg_G->default_FG_RATIO;
        }
        sp->ldim_fine = sp->ldim ;



        /*Filter and interpolate local potential into the internal log grid*/
        Zv = sp->zvalence;
        rc = sp->rc;

        /* Generate the difference potential */
        for (int idx = 0; idx < sp->rg_points; idx++)
            work[idx] = sp->vloc0[idx] + Zv * boost::math::erf (sp->r[idx] / rc) / sp->r[idx];

        /* Write raw beta function into file if requested*/
        double *bessel_rg = new double[(ct.max_l+2) * RADIAL_GVECS * sp->rg_points];

        // Reset sp->rg_points to correspond to the closest value to 10.0
        int ii;
        double max_r = std::max(sp->nlradius, sp->lradius);
        max_r = std::max(max_r, 10.0);
        for(ii = 0;ii < sp->rg_points;ii++) if(sp->r[ii] >= max_r) break;
        sp->rg_points = ii;


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
        A->FilterPotential(work, sp->r, sp->rg_points, sp->lradius, ct.rhocparm, sp->localig,
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

        // Next we want to fourier filter the input atomic charge density and transfer
        // it to the interpolation grid for use by LCAO starts and force correction routines
        for (int idx = 0; idx < sp->rg_points; idx++)
        {
            if(sp->atomic_rho[idx] < 0.0)
                work[idx] = 0.0;
            else 
                work[idx] = sqrt(sp->atomic_rho[idx]);
        }
        A->FilterPotential(work, sp->r, sp->rg_points, sp->lradius, ct.rhocparm, sp->arho_lig,
                sp->rab, 0, sp->gwidth, 0.66*sp->lradius, sp->rwidth, ct.hmingrid/(double)ct.FG_RATIO);
        for (int idx = 0; idx < MAX_LOGGRID; idx++) sp->arho_lig[idx] *= sp->arho_lig[idx];


        /*Open file for writing beta function*/
        if (pct.gridpe == 0 && write_flag)
        {
            snprintf (newname, MAX_PATH, "beta_%s.xmgr", sp->atomic_symbol);
            if(NULL == (psp = fopen (newname, "w+")))
                throw RmgFatalException() << "Unable to open beta function graph file " << " in " << __FILE__ << " at line " << __LINE__ << "\n";
        }

        for (int ip = 0; ip < sp->nbeta; ip++)
        {

            if (pct.gridpe == 0 && write_flag)
            {
                for (int idx = 0; idx < sp->kkbeta; idx++) fprintf (psp, "%e  %e\n", sp->r[idx], sp->beta[ip][idx]);
                fprintf (psp, "\n&&\n");
            }

            sp->beta_g[ip] = new double[RADIAL_GVECS];
            A->RLogGridToGLogGrid(&sp->beta[ip][0], sp->r, sp->rab, sp->beta_g[ip],
                    sp->rg_points, sp->llbeta[ip], bessel_rg);

            // Raw beta function from pp is no longer used so free it's memory
            delete [] sp->beta[ip];

        }                       /* end for ip */

        if (pct.gridpe == 0 && write_flag) fclose (psp);


        /* Now take care of the core charge if nlcc is being used */
        if (sp->nlccflag)
        {
            /*Open file for writing core charge */
            if (pct.gridpe == 0 && write_flag)
            {
                snprintf (newname, MAX_PATH, "nlcc_%s.xmgr", sp->atomic_symbol);
                if(NULL == (psp = fopen (newname, "w+")))
                    throw RmgFatalException() << "Unable to open core charge graph file " << " in " << __FILE__ << " at line " << __LINE__ << "\n";
            }

            for (int idx = 0; idx < sp->rg_points; idx++)
                work[idx] = sp->rspsco[idx] / (4.0 * PI);

            A->RLogGridToGLogGrid(work, sp->r, sp->rab, sp->rhocore_g,
                    sp->rg_points, 0, bessel_rg);

            /*Write raw (pre-filtered) data to a file if requested */
            if (pct.gridpe == 0 && write_flag)
            {
                snprintf (newname, MAX_PATH, "nlcc_%s.xmgr", sp->atomic_symbol);
                if(NULL == (psp = fopen (newname, "w+")))
                    throw RmgFatalException() << "Unable to open beta function graph file " << " in " << __FILE__ << " at line " << __LINE__ << "\n";
                for (int idx = 0; idx < sp->rg_points; idx++)
                    fprintf (psp, "%e  %e\n", sp->r[idx], work[idx]);
                fprintf (psp, "\n&&\n");
            }

            double nlccradius = 3.5 * A->GetRange(work, sp->r, sp->rab, sp->rg_points, 0.999999999);

            // Make adjustments so radii terminates on a grid point
            int nlccdim = Radius2grid (nlccradius, ct.hmingrid/(double)Rmg_G->default_FG_RATIO, Rmg_L.get_ibrav_type(), ct.localize_localpp);
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
                fclose (psp);
            }

        }                       /* end if */


        /*Open file for writing atomic orbitals */
        if (pct.gridpe == 0 && write_flag)
        {
            snprintf (newname, MAX_PATH, "atomicwf_%s.xmgr", sp->atomic_symbol);
            if(NULL == (psp = fopen (newname, "w+")))
                throw RmgFatalException() << "Unable to open atomic orbital graph file " << " in " << __FILE__ << " at line " << __LINE__ << "\n";

        }


        int ldaU_orbitals = 0;
        int ldaU_l = 0;
        sp->awave_is_ldaU = new bool[sp->num_atomic_waves_m]();

        int lm_index = 0;
        for (int ip = 0; ip < sp->num_atomic_waves; ip++)
        {

            // First get the g-space representation of the radial part of the atomic wave
            sp->atomic_wave_g[ip] = new double[RADIAL_GVECS];
            A->RLogGridToGLogGrid(sp->atomic_wave[ip], sp->r, sp->rab, sp->atomic_wave_g[ip],
                    sp->rg_points, sp->atomic_wave_l[ip], bessel_rg);


            // The range over which an atomic orbital is non-zero can vary widely so we
            // compute a range here for when we are using them directly in real space
            // (e.g. initialization and possibley LDA+U stuff)
            sp->aradius[ip]  = A->GetRange(sp->atomic_wave[ip], sp->r, sp->rab, sp->rg_points, 0.9999);
            if(sp->aradius[ip] < 6.0) sp->aradius[ip] = 6.0;

            /* Get adim_wave */
            sp->adim_wave = Radius2grid (sp->aradius[ip], ct.hmingrid, Rmg_L.get_ibrav_type(), (ct.atomic_orbital_type == LOCALIZED));
            sp->adim_wave = sp->adim_wave/2*2 - 1;

            if ((sp->adim_wave >= get_NX_GRID()) || (sp->adim_wave >= get_NY_GRID()) || (sp->adim_wave >= get_NZ_GRID()))
            {
                if(ct.atomic_orbital_type == LOCALIZED)
                {
                    if(pct.gridpe == 0) 
                    {
                        rmg_printf("Warning: localized atomic orbitals selected but their diameter exceeds cell size. Switching to delocalized.\n");
                        printf("Warning: localized atomic orbitals selected but their diameter exceeds cell size. Switching to delocalized.\n");
                    }
                    ct.atomic_orbital_type = DELOCALIZED;
                }
            }

            sp->aradius[ip] = 0.5 * ct.hmingrid * (double)(sp->adim_wave);

            double rcut = 0.8 * sp->aradius[ip];
            for(int idx = 0; idx < sp->rg_points; idx++)
            {

                if(sp->r[idx] > rcut) 
                {
                    double t1 = (sp->r[idx] - rcut) / (sp->aradius[ip] - rcut);
                    if(t1 > 1.0) t1 = 1.0; 
                    sp->atomic_wave[ip][idx] *= (1.0-t1);
                }
            }
            if (pct.gridpe == 0 && write_flag)
            {
                for (int idx = 0; idx < sp->rg_points; idx++) fprintf (psp, "%e  %e\n", sp->r[idx], sp->atomic_wave[ip][idx]);
                fprintf (psp, "\n&&\n");
            }


            // Now lda+u stuff
            int l = sp->atomic_wave_l[ip];
            int m = 2*l + 1;
            if(ct.ldaU_mode != LDA_PLUS_U_NONE)
            {
                if(l > 1)  // d or f state
                {
                    // Ignore filled shell
                    double full_occ = 2.0 * (double)m;
                    double tol = fabs(full_occ - sp->atomic_wave_oc[ip]);

                    // Only use one set of localized orbitals per species
                    if((tol > 1.0e-5) && (ldaU_l == 0))
                    {
                        ldaU_orbitals += m;                        
                        for(int im=lm_index;im < lm_index+m;im++) sp->awave_is_ldaU[im] = true;
                        ldaU_l = l;
                    }
                } 
            }
            lm_index += 2*l + 1;
        }

        ct.max_ldaU_orbitals = std::max(ct.max_ldaU_orbitals, ldaU_orbitals);
        ct.max_ldaU_l = std::max(ct.max_ldaU_l, ldaU_l);
        sp->num_ldaU_orbitals = ldaU_orbitals;
        sp->ldaU_l = ldaU_l;
        if (pct.gridpe == 0 && write_flag) fclose (psp);

        delete [] bessel_rg;
        delete [] work;

    }                           /* end for */
    delete A;

} // end InitPseudo
