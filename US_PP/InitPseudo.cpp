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
#include "RmgException.h"
#include <boost/math/special_functions/erf.hpp>

using boost::math::policies::policy;
using boost::math::policies::promote_double;
typedef policy<promote_double<false> > bessel_policy;

void SPECIES::InitPseudo (Lattice &L, BaseGrid *G, bool write_flag)
{

    if(!std::strcmp(this->atomic_symbol, "DLO")) return;

    RmgTimer *RT0 = new RmgTimer("2-Init: radial potentials");

    char newname[MAX_PATH];
    FILE *psp = NULL;
    double *rgrid = GetRgrid();

    double *work = new double[std::max(MAX_LOGGRID, this->rg_points)];
    if (pct.gridpe == 0 && write_flag)
    {
        snprintf (newname, MAX_PATH, "local_%s.xmgr",  this->atomic_symbol);
        if(NULL == (psp = fopen (newname, "w+")))
           throw RmgFatalException() << "Unable to open pseudopotential graph file " << " in " << __FILE__ << " at line " << __LINE__ << "\n";

    }


    // Might need to adjust this depending on filtering changes. Also assumes that all beta have roughly the same range
    this->nlradius = ct.projector_expansion_factor * 4.5 * GetRange(&this->beta[0][0], this->r, this->rab, this->rg_points, 0.999999999);
    this->nlradius = std::max(this->nlradius, ct.min_nlradius);
    this->nlradius = std::min(this->nlradius, ct.max_nlradius);


    /*Get nldim */
    this->nldim = Radius2grid (this->nlradius, ct.hmingrid, L.get_ibrav_type(), ct.localize_projectors);
    this->nldim = this->nldim/2*2 + 1;

    int printed = 0;
    while ((this->nldim >= get_NX_GRID()) || (this->nldim >= get_NY_GRID()) || (this->nldim >= get_NZ_GRID()))
    {
        this->nlradius -= ct.hmingrid;
        this->nldim = Radius2grid (this->nlradius, ct.hmingrid, L.get_ibrav_type(), ct.localize_projectors);
        this->nldim = this->nldim/2*2 + 1;
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

    this->nlradius = 0.5 * ct.hmingrid * (double)(this->nldim);

    // If projectors will span the full wavefunction grid then use a larger value for the nlradius for all remaining operations
    if(!ct.localize_projectors) {
        this->nlradius = std::max(this->nlradius, 7.0);
        this->nldim = Radius2grid (this->nlradius, ct.hmingrid, L.get_ibrav_type(), ct.localize_projectors);
        this->nldim = this->nldim/2*2 +1;
    }

    /*ct.max_nlpoints is max of nldim*nldim*nldim for all species */
    if (ct.max_nlpoints < (this->nldim * this->nldim * this->nldim))
        ct.max_nlpoints = this->nldim * this->nldim * this->nldim;

    if(!ct.localize_projectors) {
        ct.max_nlpoints = get_NX_GRID() * get_NY_GRID() * get_NZ_GRID();
    }

    if(ct.localize_projectors)
    {
        // Grid object local to this MPI process
        this->OG = new BaseGrid(this->nldim, this->nldim, this->nldim, 1, 1, 1, 0, 1);
        BaseGrid *OG = (BaseGrid *)this->OG;
        OG->set_rank(0, pct.my_comm);

        // Lattice object for the localized projector region. Global vectors need to be
        // scaled so that they correspond to the local region.
        Lattice *LL = new Lattice();
        double a0[3], a1[3], a2[3], s1, celldm[6], omega;
        s1 = (double)this->nldim / (double)G->get_NX_GRID(1);
        a0[0] = s1*L.get_a0(0);a0[1] = s1*L.get_a0(1);a0[2] = s1*L.get_a0(2);
        s1 = (double)this->nldim / (double)G->get_NY_GRID(1);
        a1[0] = s1*L.get_a1(0);a1[1] = s1*L.get_a1(1);a1[2] = s1*L.get_a1(2);
        s1 = (double)this->nldim / (double)G->get_NZ_GRID(1);
        a2[0] = s1*L.get_a2(0);a2[1] = s1*L.get_a2(1);a2[2] = s1*L.get_a2(2);

        LL->set_ibrav_type(None);
        LL->latgen(celldm, &omega, a0, a1, a2, true);

        this->prj_pwave = new Pw(*OG, *LL, 1, false);
    }
    else
    {
        this->prj_pwave = new Pw(*G, L, 1, false);
    }


    this->ldim = Radius2grid (this->lradius, ct.hmingrid / (double)G->default_FG_RATIO, L.get_ibrav_type(), ct.localize_localpp);
    this->ldim = this->ldim/2 * 2 +1;
    if(!ct.localize_localpp)
    {
        this->lradius = 10.0;
        this->ldim = Radius2grid (this->lradius, ct.hmingrid / (double)G->default_FG_RATIO, L.get_ibrav_type(), ct.localize_localpp);
        this->ldim = this->nldim/2*2 + 1;
    }
    else
    {
        this->lradius = 0.5 * ct.hmingrid * (double)(this->ldim-1) / (double)G->default_FG_RATIO;
        this->lradius -= 0.5 * ct.hmingrid / (double)G->default_FG_RATIO;
    }
    this->ldim_fine = this->ldim ;


    /*Filter and interpolate local potential into the internal log grid*/
    double Zv = this->zvalence;
    double rc = this->rc;

    /* Generate the difference potential */
    for (int idx = 0; idx < this->rg_points; idx++)
        work[idx] = this->vloc0[idx] + Zv * boost::math::erf (this->r[idx] / rc) / this->r[idx];

    /* Write raw beta function into file if requested*/
    double *bessel_rg = new double[(ct.max_l+2) * RADIAL_GVECS * this->rg_points];

    // Reset this->rg_points to correspond to the closest value to 10.0
    int ii;
    double max_r = std::max(this->nlradius, this->lradius);
    max_r = std::max(max_r, 10.0);
    for(ii = 0;ii < this->rg_points;ii++) if(this->r[ii] >= max_r) break;
    this->rg_points = ii;


    RmgTimer *RT1 = new RmgTimer("radial beta");
    InitBessel(this->r, this->rg_points, ct.max_l + 1, bessel_rg);
    delete RT1;

    // get the local pseudopotential in G space
    this->localpp_g = new double[RADIAL_GVECS];
    this->der_localpp_g = new double[RADIAL_GVECS];
    this->arho_g = new double[RADIAL_GVECS];
    this->rhocore_g = new double[RADIAL_GVECS];
    RLogGridToGLogGrid(work, this->r, this->rab, this->localpp_g,
            this->rg_points, 0, bessel_rg);
    RLogGridToGLogGrid(this->atomic_rho, this->r, this->rab, this->arho_g,
            this->rg_points, 0, bessel_rg);
    Der_Localpp_g(work, this->r, this->rab, this->der_localpp_g, this->rg_points);

    if (pct.gridpe == 0 && write_flag)
    {
        for (int idx = 0; idx < this->rg_points; idx++)
            fprintf (psp, "%e  %e\n", this->r[idx], work[idx]);

        fprintf (psp, "\n&&\n");
    }


    // Transform to g-space and filter it with filtered function returned on standard log grid
    FilterPotential(work, this->r, this->rg_points, this->lradius, ct.rhocparm, this->localig,
            this->rab, 0, this->gwidth, 0.66*this->lradius, this->rwidth, ct.hmingrid / (double)G->default_FG_RATIO);

    /*Write local projector into a file if requested*/
    if ((pct.gridpe == 0) && write_flag)
    {
        for (int idx = 0; idx < MAX_LOGGRID; idx++)
            fprintf (psp, "%e  %e \n", rgrid[idx], this->localig[idx]);

        /* output xmgr data separator */
        fprintf (psp, "\n&&\n");
        fclose (psp);
    }

    // Get the G=0 component
    double fac = Zv * 4.0 * PI;
    for (int idx = 0; idx < this->rg_points; idx++)
        work[idx] = this->vloc0[idx] + Zv/this->r[idx];
    double G0 = 4.0*PI*radint1 (work, this->r, this->rab, this->rg_points);
    this->localpp_g[0] = G0;

    // Subtract off analytic transform of erf that was added above
    for(int idx=1;idx < RADIAL_GVECS;idx++)
    {
        double gval = gvec[idx]*gvec[idx];
        this->localpp_g[idx] -= fac * exp ( - 0.25*gval*rc*rc) / gval;

        this->der_localpp_g[idx] += fac * exp(-0.25*gval*rc*rc)/(gval * gval) * (0.25 * gval * rc * rc + 1.0);
    }

    // Next we want to fourier filter the input atomic charge density and transfer
    // it to the interpolation grid for use by LCAO starts and force correction routines
    for (int idx = 0; idx < this->rg_points; idx++)
    {
        if(this->atomic_rho[idx] < 0.0)
            work[idx] = 0.0;
        else 
            work[idx] = sqrt(this->atomic_rho[idx]);
    }
    FilterPotential(work, this->r, this->rg_points, this->lradius, ct.rhocparm, this->arho_lig,
            this->rab, 0, this->gwidth, 0.66*this->lradius, this->rwidth, ct.hmingrid/(double)ct.FG_RATIO);
    for (int idx = 0; idx < MAX_LOGGRID; idx++) this->arho_lig[idx] *= this->arho_lig[idx];


    /*Open file for writing beta function*/
    if (pct.gridpe == 0 && write_flag)
    {
        snprintf (newname, MAX_PATH, "beta_%s.xmgr", this->atomic_symbol);
        if(NULL == (psp = fopen (newname, "w+")))
            throw RmgFatalException() << "Unable to open beta function graph file " << " in " << __FILE__ << " at line " << __LINE__ << "\n";
    }

    this->rbeta_g.resize(boost::extents[this->nbeta][MAX_L+1]);
    for (int ip = 0; ip < this->nbeta; ip++)
    {
        if (pct.gridpe == 0 && write_flag)
        {
            for (int idx = 0; idx < this->kkbeta; idx++) fprintf (psp, "%e  %e\n", this->r[idx], this->beta[ip][idx]);
            fprintf (psp, "\n&&\n");
        }

        this->beta_g.emplace_back(new double[RADIAL_GVECS]);

        RLogGridToGLogGrid(&this->beta[ip][0], this->r, this->rab, this->beta_g[ip].get(),
                this->rg_points, this->llbeta[ip], bessel_rg);

        for(int idx = 0; idx < this->rg_points; idx++)
            work[idx] = this->beta[ip][idx] * this->r[idx];

        // generate beta * x, beta * y, beta * z for stress calculation 
        // equally, the beta * r will have the angular momentum +1
        // l1 = llbeta[ip], l2 = 1 (for x, y, z)
        // L = l1 + l2 ... |l1-l2|
        for(int L = std::abs(this->llbeta[ip] - 1); L <= this->llbeta[ip] +1; L++)
        {
            
            this->rbeta_g[ip][L] = new double[RADIAL_GVECS];
            RLogGridToGLogGrid(work, this->r, this->rab, this->rbeta_g[ip][L],
                    this->rg_points, L, bessel_rg);
        }

    }                       /* end for ip */

    if (pct.gridpe == 0 && write_flag) fclose (psp);

    // Raw beta functions from pp are no longer used so free  memory
    this->beta.clear();


    /* Now take care of the core charge if nlcc is being used */
    if (this->nlccflag)
    {
        /*Open file for writing core charge */
        if (pct.gridpe == 0 && write_flag)
        {
            snprintf (newname, MAX_PATH, "nlcc_%s.xmgr", this->atomic_symbol);
            if(NULL == (psp = fopen (newname, "w+")))
                throw RmgFatalException() << "Unable to open core charge graph file " << " in " << __FILE__ << " at line " << __LINE__ << "\n";
        }

        for (int idx = 0; idx < this->rg_points; idx++)
            work[idx] = this->rspsco[idx] / (4.0 * PI);

        RLogGridToGLogGrid(work, this->r, this->rab, this->rhocore_g,
                this->rg_points, 0, bessel_rg);

        /*Write raw (pre-filtered) data to a file if requested */
        if (pct.gridpe == 0 && write_flag)
        {
            snprintf (newname, MAX_PATH, "nlcc_%s.xmgr", this->atomic_symbol);
            if(NULL == (psp = fopen (newname, "w+")))
                throw RmgFatalException() << "Unable to open beta function graph file " << " in " << __FILE__ << " at line " << __LINE__ << "\n";
            for (int idx = 0; idx < this->rg_points; idx++)
                fprintf (psp, "%e  %e\n", this->r[idx], work[idx]);
            fprintf (psp, "\n&&\n");
        }

        double nlccradius = 3.5 * GetRange(work, this->r, this->rab, this->rg_points, 0.999999999);

        // Make adjustments so radii terminates on a grid point
        int nlccdim = Radius2grid (nlccradius, ct.hmingrid/(double)G->default_FG_RATIO, L.get_ibrav_type(), ct.localize_localpp);
        nlccdim = nlccdim/2*2 + 1;
        nlccradius = 0.5 * ct.hmingrid * (double)(nlccdim-1) / (double)G->default_FG_RATIO;
        nlccradius -= 0.5 * ct.hmingrid / (double)G->default_FG_RATIO;
        double nlcccut = 0.66 * nlccradius;

        FilterPotential(work, this->r, this->rg_points, nlccradius, ct.rhocparm, &this->rhocorelig[0],
                this->rab, 0, this->gwidth, nlcccut, this->rwidth, ct.hmingrid/(double)G->default_FG_RATIO);

        /*Oscilations at the tail end of filtered function may cause rhocore to be negative
         * but I am not sure if this is the right solution, it may be better to fix charge density
         * this rhocore less smooth*/
        for (int idx = 0; idx < MAX_LOGGRID; idx++)
            if (this->rhocorelig[idx] < 0.0)
                this->rhocorelig[idx] = 0.0;

        /*Write filtered data to a file if requested */
        if (pct.gridpe == 0 && write_flag)
        {
            for (int idx = 0; idx < MAX_LOGGRID; idx++)
                fprintf (psp, "%e  %e\n", rgrid[idx], this->rhocorelig[idx]);
            fprintf (psp, "\n&&\n");
            fclose (psp);
        }

    }                       /* end if */


    /*Open file for writing atomic orbitals */
    if (pct.gridpe == 0 && write_flag)
    {
        snprintf (newname, MAX_PATH, "atomicwf_%s.xmgr", this->atomic_symbol);
        if(NULL == (psp = fopen (newname, "w+")))
            throw RmgFatalException() << "Unable to open atomic orbital graph file " << " in " << __FILE__ << " at line " << __LINE__ << "\n";

    }


    this->awave_is_ldaU = new bool[this->num_atomic_waves_m]();

    int lm_index = 0;
    this->num_ldaU_orbitals = 0;
    this->atomic_wave_g.resize(this->num_atomic_waves);
    for (int ip = 0; ip < this->num_atomic_waves; ip++)
    {

        // First get the g-space representation of the radial part of the atomic wave
        this->atomic_wave_g[ip] = new double[RADIAL_GVECS];
        RLogGridToGLogGrid(this->atomic_wave[ip], this->r, this->rab, this->atomic_wave_g[ip],
                this->rg_points, this->atomic_wave_l[ip], bessel_rg);


        // The range over which an atomic orbital is non-zero can vary widely so we
        // compute a range here for when we are using them directly in real space
        // (e.g. initialization and possibley LDA+U stuff)
        this->aradius[ip]  = GetRange(this->atomic_wave[ip], this->r, this->rab, this->rg_points, 0.9999);
        if(this->aradius[ip] < 6.0) this->aradius[ip] = 6.0;

        /* Get adim_wave */
        this->adim_wave = Radius2grid (this->aradius[ip], ct.hmingrid, L.get_ibrav_type(), (ct.atomic_orbital_type == LOCALIZED));
        this->adim_wave = this->adim_wave/2*2 - 1;

        if ((this->adim_wave >= get_NX_GRID()) || (this->adim_wave >= get_NY_GRID()) || (this->adim_wave >= get_NZ_GRID()))
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

        this->aradius[ip] = 0.5 * ct.hmingrid * (double)(this->adim_wave);

        double rcut = 0.8 * this->aradius[ip];
        for(int idx = 0; idx < this->rg_points; idx++)
        {

            if(this->r[idx] > rcut) 
            {
                double t1 = (this->r[idx] - rcut) / (this->aradius[ip] - rcut);
                if(t1 > 1.0) t1 = 1.0; 
                this->atomic_wave[ip][idx] *= (1.0-t1);
            }
        }
        if (pct.gridpe == 0 && write_flag)
        {
            for (int idx = 0; idx < this->rg_points; idx++) fprintf (psp, "%e  %e\n", this->r[idx], this->atomic_wave[ip][idx]);
            fprintf (psp, "\n&&\n");
        }


        // Now lda+u stuff
        int l = this->atomic_wave_l[ip];
        int m = 2*l + 1;
        if(ct.ldaU_mode != LDA_PLUS_U_NONE)
        {
            if(!strcasecmp(this->ldaU_label.c_str(), this->atomic_wave_label[ip].c_str() )) 
            {
                this->ldaU_l = l;
                for(int im=lm_index;im < lm_index+m;im++) 
                {
                    this->awave_is_ldaU[im] = true;
                }
                this->num_ldaU_orbitals = m;
            } 
        }
        lm_index += 2*l + 1;
    }

    ct.max_ldaU_orbitals = std::max(ct.max_ldaU_orbitals, this->num_ldaU_orbitals);
    ct.max_ldaU_l = std::max(ct.max_ldaU_l, this->ldaU_l);
    if (pct.gridpe == 0 && write_flag) fclose (psp);

    delete RT0;

    // If this is an USPP we need to setup the radial qfunction stuff
    if(!ct.norm_conserving_pp)
    {
        RmgTimer RT1("2-Init: qfunct");
        FILE *fqq = NULL;
        
        if (write_flag)
        {
            snprintf (newname, MAX_PATH, "q_%s.xmgr", this->atomic_symbol);
            if (pct.gridpe == 0)
            {
                if(NULL == (fqq = fopen (newname, "w+")))
                    throw RmgFatalException() << "Unable to open pseudopotential graph file " << " in " << __FILE__ << " at line " << __LINE__ << "\n";

            }
        }

        //  change from qnm to qnm_l
        if(!this->q_with_l)
        {
            for (int i = 0; i < this->nbeta; i++)
            {
                int il = this->llbeta[i];
                for (int j = i; j < this->nbeta; j++)
                {
                    int jl = this->llbeta[j];
                    int idx = j * (j + 1) / 2 + i;
                    double *qnm_tpr = this->qnm + idx * MAX_RGRID;
                    for (int ll = abs (il - jl); ll <= il + jl; ll = ll + 2)
                    {

                        for (int k = 0; k < MAX_RGRID; k++)
                            this->qnm_l[k + (idx * this->nlc + ll) * MAX_RGRID ] = 0.0;

                        for (int k = 0; k < this->rg_points; k++)
                        {
                            this->qnm_l[k + (idx * this->nlc + ll) * MAX_RGRID ] = qnm_tpr[k]/this->r[k]/this->r[k];
                            if(this->nqf > 0 && this->r[k] < this->rinner[ll])
                                this->qnm_l[k + (idx * this->nlc + ll) * MAX_RGRID ] = get_QnmL (idx, ll, this->r[k], this);
                        }
                    }
                }
            }
        }
        else
        {
            for(int i = 0; i < (this->nlc * this->nbeta *(this->nbeta +1)) /2; i++)
            for (int k = 0; k < this->rg_points; k++)
            {
                this->qnm_l[k + i * MAX_RGRID ] /= this->r[k]*this->r[k];
            }
        }

        double qradius = ct.min_qradius;
        for(int i = 0; i < (this->nlc * this->nbeta *(this->nbeta +1)) /2; i++)
        {
            double *qnm_tpr = &this->qnm_l[i * MAX_RGRID];
            this->qradius = 2.5 * GetRange(qnm_tpr, this->r, this->rab, this->rg_points, 0.999999999);
            this->qradius = std::max(this->qradius, qradius);
            qradius = this->qradius;

        }

        this->qradius = std::min(this->qradius, ct.max_qradius);

        // Make adjustments so radii terminates on a grid point
        this->qdim = Radius2grid (this->qradius, ct.hmingrid/(double)Rmg_G->default_FG_RATIO, Rmg_L.get_ibrav_type(), false);
        this->qdim = this->qdim/2*2 + 1;

        this->qcut = 0.5 * this->qradius;

        if (this->qdim >= get_FNX_GRID()) this->qdim = get_FNX_GRID();
        if (this->qdim >= get_FNY_GRID()) this->qdim = get_FNY_GRID();
        if (this->qdim >= get_FNZ_GRID()) this->qdim = get_FNZ_GRID();

        if (ct.max_Qpoints < (this->qdim * this->qdim * this->qdim))
            ct.max_Qpoints = this->qdim * this->qdim * this->qdim;

        this->qnmlig.clear();
        int idx = 0;
        for (int i = 0; i < this->nbeta; i++)
        {
            int il = this->llbeta[i];
            for (int j = i; j < this->nbeta; j++)
            {
                int jl = this->llbeta[j];
                idx = j * (j + 1) / 2 + i;
                for (int ll = abs (il - jl); ll <= il + jl; ll = ll + 2)
                {
                    double *qnm_tpr = this->qnm_l + (idx * this->nlc + ll) *  MAX_RGRID;
                    this->qnmlig[qnm_key(i, j, ll)] = {};
                    this->qnmlig[qnm_key(i, j, ll)].resize(MAX_LOGGRID);

                    if (pct.gridpe == 0 && write_flag)
                    {
                        for (int k = 0; k < this->kkbeta; k++)
                            fprintf (fqq, "%e  %e\n", this->r[k], qnm_tpr[k]);
                        fprintf (fqq, "&\n");
                    }

                    FilterPotential(qnm_tpr, this->r, this->rg_points, this->qradius, ct.rhocparm, this->qnmlig[qnm_key(i, j, ll)].data(),
                             this->rab, ll, this->gwidth, this->qcut, 1.0, ct.hmingrid/(double)Rmg_G->default_FG_RATIO);

                    /*Write final filtered Q function if requested*/
                    if (pct.gridpe == 0 && write_flag)
                    {
                        for (int k = 0; k < MAX_LOGGRID; k++)
                        {
                            fprintf (fqq, "%e  %e\n", rgrid[k], this->qnmlig[qnm_key(i, j, ll)][k]);
                        }
                        fprintf (fqq, "&\n");
                    }

                }               /*end for ll */
            }                   /*end for j */
        }                       /*end for i */

        // Raw q function from pp is no longer needed so free it's memory
        delete []  this->qnm;

        if (pct.gridpe == 0 && write_flag)
        {
            fclose (fqq);
        }

    } // end if for !ct.norm_conserving_pp

    delete [] bessel_rg;
    delete [] work;

    if(ct.ldaU_mode != LDA_PLUS_U_NONE && ct.max_ldaU_orbitals == 0)
         throw RmgFatalException() << "LDA+U: no U assigned" << " in " << __FILE__ << " at line " << __LINE__ << "\n";

} // end InitPseudo
