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


#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "portability.h"
#include "transition.h"
#include "const.h"
#include "RmgTimer.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "RmgException.h"
#include "Kpoint.h"
#include "transition.h"
#include "Plots.h"
#include "../Headers/prototypes.h"
#include "../Headers/macros.h"

// Instantiate gamma and non-gamma versions
template bool Quench<double> (double *, double *, double *, double *, double *, double *, double *, Kpoint<double> **Kptr);
template bool Quench<std::complex<double> > (double *, double *, double *, double *, double *, double *, double *, Kpoint<std::complex<double>> **Kptr);

template <typename OrbitalType> bool Quench (double * vxc, double * vh, double * vnuc, double * rho,
             double * rho_oppo, double * rhocore, double * rhoc, Kpoint<OrbitalType> **Kptr)
{

    bool CONVERGED;
    std::vector<double> RMSdV;

    int numacc = 1, ic;
    /*int ist, ik;
       double KE; */

    /* ---------- begin scf loop ---------- */
    
    double start_time = my_crtc ();
    double elapsed_time;

    for (ct.scf_steps = 0, CONVERGED = false;
         ct.scf_steps < ct.max_scf_steps && !CONVERGED; ct.scf_steps++, ct.total_scf_steps++)
    {

        elapsed_time = my_crtc() - start_time;
        if (pct.imgpe == 0)
            rmg_printf ("\n\nquench: ----- [md: %3d/%-d  scf: %3d/%-d  scf time: %8.2f secs  RMS[dV]: %8.2e ] -----\n",
                    ct.md_steps, ct.max_md_steps, ct.scf_steps, ct.max_scf_steps, elapsed_time, ct.rms);


        /* perform a single self-consistent step */
        CONVERGED = Scf (vxc, vh, ct.vh_ext, vnuc, rho, rho_oppo, rhocore, rhoc, ct.spin_flag, ct.hartree_min_sweeps, ct.hartree_max_sweeps, ct.boundaryflag, Kptr, RMSdV);


	/* output the eigenvalues with occupations */
	if (ct.write_eigvals_period)
	{
	    if (ct.scf_steps % ct.write_eigvals_period == 0)
	    {
		if (pct.imgpe == 0)
		{
                    OutputEigenvalues (Kptr, 0, ct.scf_steps);
		    rmg_printf ("\nTotal charge in supercell = %16.8f\n", ct.tcharge);
		}
	    }
	}

#if PLPLOT_LIBS
        if(pct.imgpe == 0) {
            std::vector<double> x;
            std::string ConvergencePlot(ct.basename);
            ConvergencePlot = ConvergencePlot + ".rmsdv.png";
            std::string ConvergenceTitle("RMG convergence:");
            if(CONVERGED) {
                ConvergenceTitle = ConvergenceTitle + " quench completed.";
            }
            else {
                ConvergenceTitle = ConvergenceTitle + " quenching electrons.";
            }
            LinePlotLog10y(ConvergencePlot.c_str(), "SCF Steps", "log10(RMS[dV])", ConvergenceTitle.c_str(), x, RMSdV);
        }
#endif

    }

    /* ---------- end scf loop ---------- */

    if (CONVERGED)
    {
	rmg_printf ("\n");
	//progress_tag ();
        rmg_printf ("\n\nquench: ----- [md: %3d/%-d  scf: %3d/%-d  scf time: %8.2f secs  RMS[dV]: %8.2e ] -----\n",
                    ct.md_steps, ct.max_md_steps, ct.scf_steps, ct.max_scf_steps, elapsed_time, ct.rms);
	rmg_printf ("potential convergence has been achieved. stopping ...\n");
	    
	/*Write PDOS if converged*/
//	if (ct.pdos_flag)
//	    get_pdos (Kptr[0]->kstates, ct.Emin, ct.Emax, ct.E_POINTS);
    }

    rmg_printf ("\n");
    progress_tag ();
    rmg_printf ("final total energy = %14.7f Ha\n", ct.TOTAL);



    /* output final eigenvalues with occupations */

    OutputEigenvalues (Kptr, 0, ct.scf_steps);
    rmg_printf ("\nTotal charge in supercell = %16.8f\n", ct.tcharge);

//    wvfn_residual (Kptr[0]->kstates);


    // Output RMSdV for convergence analysis
    if(pct.imgpe == 0) {
        // Write out convergence file
        std::string ConvergenceFile(ct.basename);
        ConvergenceFile = ConvergenceFile + ".rmsdv.xmgr";
        mode_t mode = O_CREAT |O_TRUNC |O_RDWR;
        if(ct.md_steps > 0) mode = O_RDWR | O_APPEND;

        int fhand = open(ConvergenceFile.c_str(), mode, S_IREAD | S_IWRITE);
        if (fhand < 0)
            throw RmgFatalException() <<  "Unable to write file in " << __FILE__ << " at line " << __LINE__ << "\n";
        char tbuf[100];

        int idx = 0;
        snprintf(tbuf, sizeof(tbuf), "@    title \"RMG convergence plot\"\n@    xaxis  label \"SCF steps\"@    yaxis  label \"log10(RMSdV)\"\n");

        for(auto it = RMSdV.begin();it != RMSdV.end();it++) {
            snprintf(tbuf, sizeof(tbuf), "%d %12.6f\n", idx, log10(*it));
            write(fhand, tbuf, strlen(tbuf));
            idx++;
        }

    }
        


    /*When running MD, force pointers need to be rotated before calculating new forces */
    if ((ct.forceflag == MD_CVE) || (ct.forceflag == MD_CVT) || (ct.forceflag == MD_CPT))
    {

	/* rotate the force pointers */
	switch (ct.mdorder)
	{
	    case ORDER_2:
		numacc = 1;
		break;
	    case ORDER_3:
		numacc = 2;
		break;
	    case ORDER_5:
		numacc = 4;
		break;
	}
	for (ic = (numacc - 1); ic > 0; ic--)
	{
	    ct.fpt[ic] = ct.fpt[ic - 1];
	}
	ct.fpt[0] = ct.fpt[numacc - 1] + 1;
	if (ct.fpt[0] > (numacc - 1) || numacc == 1)
	    ct.fpt[0] = 0;

    }


    /* compute the forces */
    /* Do not calculate forces for quenching when we are not converged */
//    if (CONVERGED || (ct.forceflag != MD_QUENCH))
	Force (rho, rho_oppo, rhoc, vh, vxc, vnuc, Kptr);

    /* output the forces */
    if (pct.imgpe == 0)
	write_force ();

    return CONVERGED;


}                               /* end quench */




/******/
