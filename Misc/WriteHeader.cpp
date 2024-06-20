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


#include <stdio.h>
#include <time.h>
#include <math.h>

#include "common_prototypes.h"
#include "main.h"
#include "Functional.h"
#include "RmgParallelFft.h"
#include "transition.h"
#include "Functional.h"

static void init_write_pos (void);


char *lattice_type[] = {
    "",
    "Cubic_primitive",
    "Cubic_FC",
    "Cubic_BC",
    "Hexagonal gamma = 120",
    "Trigonal_primitive",
    "Tetragonal_primitive",
    "Tetragonal_BC",
    "Orthorhombic_primitive",
    "Orthorhombic_base_centred",
    "Orthorhombic_BC",
    "Orthorhombic_FC",
    "Monoclinic_primitive",
    "Monoclinic_base_centred",
    "Triclinic_primitive",
    "Hexagonal gamma = 60"
};


/* Writes out header information */
void WriteHeader (void)
{

    int kpt, i;
    time_t tt;
    double crho_fract;
    int max_funct_length, funct_legend_length, funct_spacing, funct_padding_left, funct_padding_right;

    char *timeptr;
    time (&tt);
    timeptr = ctime (&tt);


    if(pct.imgpe==0) fprintf(ct.logfile, "\n\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "                     * * * * * * * * * *\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "                     *    R   M   G    *\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "                     * * * * * * * * * *\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "     -- A Real Space Multigrid Electronic structure code --\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "     --      More information at www.rmgdft.org          --\n");
        

    if(pct.imgpe==0) fprintf(ct.logfile, "\nCode Revision:     %s", RMG_REVISION);
#ifdef RMG_BRANCH
    if(pct.imgpe==0) fprintf(ct.logfile, "-%s", RMG_BRANCH);
#endif
#ifdef RMG_COMMIT
    if(pct.imgpe==0) fprintf(ct.logfile, "-%s", RMG_COMMIT);
#endif
    if(pct.imgpe==0) fprintf(ct.logfile, "\nBuild Date:        %s  %s", __DATE__, __TIME__);
    if(pct.imgpe==0) fprintf(ct.logfile, "\nStart time:        %s", timeptr);

    if(pct.imgpe==0) fprintf(ct.logfile, "\nNOTICE: RMG internal pseudopotentials have switched to");
    if(pct.imgpe==0) fprintf(ct.logfile, "\nONCVP from Ultrasoft. You can revert to Ultrasoft by");
    if(pct.imgpe==0) fprintf(ct.logfile, "\nadding the input tag internal_pseudo_type=\"ultrasoft\" to");
    if(pct.imgpe==0) fprintf(ct.logfile, "\nyour input files.\n\n");

    if(pct.imgpe==0) fprintf(ct.logfile, "\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "Files\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "   Control input file:        %s\n", ct.cfile);
    if (ct.runflag == 1)
	if(pct.imgpe==0) fprintf(ct.logfile, "   Data input file:           %s\n", ct.infile);
    if(pct.imgpe==0) fprintf(ct.logfile, "   Data output file:          %s\n", ct.outfile);
    
    if(pct.imgpe==0) fprintf(ct.logfile, "\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "Run Setup\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "    Calculation type:         ");
    switch (ct.forceflag)
    {
    case MD_QUENCH:
        if(pct.imgpe==0) fprintf(ct.logfile, "Quench electrons - Fixed ions SCF calculation\n");
        break;

    case MD_FASTRLX:
        if(pct.imgpe==0) fprintf(ct.logfile, "Structure Optimization.\n");
        break;

    case BAND_STRUCTURE:
        if(pct.imgpe==0) fprintf(ct.logfile, "Band structure calculation.\n");
        break;

    case MD_CVE:
        if(pct.imgpe==0) fprintf(ct.logfile, "Molecular dynamics - CVE\n");
        break;

    case MD_CVT:
        if(pct.imgpe==0) fprintf(ct.logfile, "Molecular dynamics - CVT\n");
        break;

    case MD_CPT:
        if(pct.imgpe==0) fprintf(ct.logfile, "Molecular dynamics - CPT\n");
        break;

    case PLOT:
        if(pct.imgpe==0) fprintf(ct.logfile, "Plot density in DX form.\n");
        break;

    case PSIPLOT:
        if(pct.imgpe==0) fprintf(ct.logfile, "Plot Psi^2 in DX form.\n");
        break;

    case NEB_RELAX:
        if(pct.imgpe==0) fprintf(ct.logfile, "Molecular dynamics using Nudged Elastic Band.\n");
        break;

    case TDDFT:
        if(pct.imgpe==0) fprintf(ct.logfile, "Time dependent DFT (TDDFT) calculation \n");
        break;
    case Exx_only:
        if(pct.imgpe==0) fprintf(ct.logfile, "calculate Exx integral's from saveed wave functions \n");
        break;
    case STM:
        if(pct.imgpe==0) fprintf(ct.logfile, "calculate STM charge density \n");
        break;
    case NSCF:
        if(pct.imgpe==0) fprintf(ct.logfile, "NSCF calculate \n");
        break;


    default:
        rmg_error_handler (__FILE__,__LINE__,"Unknown molecular dynamics method.");
    }
    if(pct.imgpe==0) fprintf(ct.logfile, "    Description:              %s\n", ct.description.c_str());
    if(pct.imgpe==0) fprintf(ct.logfile, "    Orbital Initialization:   ");
    switch (ct.runflag)
    {

	case RESTART:
	    if(pct.imgpe==0) fprintf(ct.logfile, "Read from %s \n", ct.infile);
	    break;

	case RANDOM_START:
	    if(pct.imgpe==0) fprintf(ct.logfile, "Random\n");
	    break;
	
	case LCAO_START:
	    if(pct.imgpe==0) fprintf(ct.logfile, "LCAO (%d LCAO and %d random orbitals)\n",  ct.total_atomic_orbitals, ct.extra_random_lcao_states);
	    break;

	case MODIFIED_LCAO_START:
	    if(pct.imgpe==0) fprintf(ct.logfile, "LCAO (%d MODIFIED LCAO and %d random orbitals)\n",  ct.init_states, ct.extra_random_lcao_states);
	    break;
	
	case INIT_FIREBALL:
	    if(pct.imgpe==0) fprintf(ct.logfile, "Fireball\n");
	    break;
	
	case INIT_GAUSSIAN:
	    if(pct.imgpe==0) fprintf(ct.logfile, "Gaussian\n");
	    break;
	
	case Start_TDDFT:
	    if(pct.imgpe==0) fprintf(ct.logfile, "Initial start TD DFT calculation\n");
	    break;
	
	case Restart_TDDFT:
	    if(pct.imgpe==0) fprintf(ct.logfile, "Restart TD DFT calculation\n");
	    break;
        
	default:
            if(pct.imgpe==0) fprintf(ct.logfile, "Unknown start mode\n");
    }
    const std::string tstr = Functional::get_dft_name_rmg();
    if(pct.imgpe==0) fprintf(ct.logfile, "    Exchange Correlation:     %s\n", tstr.c_str());
    if(pct.imgpe==0) fprintf(ct.logfile, "    Spin Polarization:        ");
    if (ct.spin_flag)
    {
	if(pct.imgpe==0) fprintf(ct.logfile, "ON\n");
    }
    else
    {
	if(pct.imgpe==0) fprintf(ct.logfile, "OFF\n");
    }
    if (fabs(ct.background_charge) > 1e-6)
    {
	if(pct.imgpe==0) fprintf(ct.logfile, "    System charge:            %6.2f\n", ct.background_charge);
    }
    else
    {
	if(pct.imgpe==0) fprintf(ct.logfile, "    System charge:            Neutral\n");
    }


    
    if(pct.imgpe==0) fprintf(ct.logfile, "\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "Processor Topology\n");  
    if(pct.imgpe==0) fprintf(ct.logfile, "   Total PEs:                 %d\n", (get_PE_X() * get_PE_Y() * get_PE_Z()));
    if(pct.imgpe==0) fprintf(ct.logfile, "   X direction:               %d\n", get_PE_X());
    if(pct.imgpe==0) fprintf(ct.logfile, "   Y direction:               %d\n", get_PE_Y());
    if(pct.imgpe==0) fprintf(ct.logfile, "   Z direction:               %d\n", get_PE_Z());
    if(pct.imgpe==0) fprintf(ct.logfile, "   MG Threads/PE:             %d\n", ct.MG_THREADS_PER_NODE);
    if(pct.imgpe==0) fprintf(ct.logfile, "   OMP Threads/PE:            %d\n", ct.OMP_THREADS_PER_NODE);
    
    
    
    if(pct.imgpe==0) fprintf(ct.logfile, "\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "Grid Points");
    if(pct.imgpe==0) fprintf(ct.logfile, "  (Linear Anisotropy: %5.3f)", get_anisotropy());
    if(pct.imgpe==0) fprintf(ct.logfile, "\n");
    
    if(pct.imgpe==0) fprintf(ct.logfile, "    X:  Total: %d   Per PE: %d   Spacing:%5.3f a0  \n", get_NX_GRID(), get_PX0_GRID(),   get_hxgrid() * get_xside());
    if(pct.imgpe==0) fprintf(ct.logfile, "    Y:  Total: %d   Per PE: %d   Spacing:%5.3f a0  \n", get_NY_GRID(), get_PY0_GRID(),   get_hygrid() * get_yside());
    if(pct.imgpe==0) fprintf(ct.logfile, "    Z:  Total: %d   Per PE: %d   Spacing:%5.3f a0  \n", get_NZ_GRID(), get_PZ0_GRID(),   get_hzgrid() * get_zside());
    if(pct.imgpe==0) fprintf(ct.logfile, "\n");

    if(ct.coalesce_states)
    {
        if(pct.imgpe==0) fprintf(ct.logfile, "    Coalescing states in X with factor %d\n", pct.coalesce_factor);
    }

    /* We compute the equivalent energy cutoff using the density of grid
     * points in the cell with a correction for the grid anisotropy.
     */
    double tpiba2 = 4.0 * PI * PI / (Rmg_L.celldm[0] * Rmg_L.celldm[0]);
    double t1 = coarse_pwaves->gcut * tpiba2;
    double t2 = fine_pwaves->gcut * tpiba2;
    int ibrav = get_ibrav_type();
    if(ibrav < 0) ibrav *= -1;
    if(pct.imgpe==0) fprintf(ct.logfile, "    Equivalent energy cutoffs (psi,rho):  %8.3f   %8.3f Ry\n", t1, t2);
    if(pct.imgpe==0) fprintf(ct.logfile, "\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "    Charge density grid:         %d times finer\n", get_FG_RATIO());


    double density[3];
    double planar_anisotropy = GetPlanarAnisotropy(density);
    if(pct.imgpe==0) fprintf(ct.logfile, "\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "Coordinate planes\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "  Planar Anisotropy: %5.3f\n", planar_anisotropy);
    if(pct.imgpe==0) fprintf(ct.logfile, "  A0-A1 density: %9.3f\n", density[0]);
    if(pct.imgpe==0) fprintf(ct.logfile, "  A0-A2 density: %9.3f\n", density[1]);
    if(pct.imgpe==0) fprintf(ct.logfile, "  A1-A2 density: %9.3f\n", density[2]);

    if(pct.imgpe==0) fprintf(ct.logfile, "\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "Lattice (Cell) Setup\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "    Type:                       %s\n", lattice_type[ibrav]);
    if(pct.imgpe==0) fprintf(ct.logfile, "    Volume (a0^3):              %8.2f\n", get_vel() * ct.psi_nbasis);
    if(pct.imgpe==0) fprintf(ct.logfile, "    Boundary conditions:        ");
    switch (ct.boundaryflag)
    {

    case PERIODIC:
        if(pct.imgpe==0) fprintf(ct.logfile, "Periodic\n");
        break;

    case CLUSTER:
        if(pct.imgpe==0) fprintf(ct.logfile, "Cluster\n");
        break;

    case SURFACE:
        if(pct.imgpe==0) fprintf(ct.logfile, "Surface\n");
        break;
	
    default:
	if(pct.imgpe==0) fprintf(ct.logfile, "Unknown boundary conditions\n");

    }                           /* end switch */
    if(pct.imgpe==0) fprintf(ct.logfile, "\n");

    // Print out original lattice vectors
    double *a0 = Rmg_L.a0i;
    double *a1 = Rmg_L.a1i;
    double *a2 = Rmg_L.a2i;
    if(pct.imgpe==0) fprintf(ct.logfile, "    X Basis Vector:  %10.3f  %10.3f  %10.3f a0\n", a0[0], a0[1], a0[2]);
    if(pct.imgpe==0) fprintf(ct.logfile, "    Y Basis Vector:  %10.3f  %10.3f  %10.3f a0\n", a1[0], a1[1], a1[2]);
    if(pct.imgpe==0) fprintf(ct.logfile, "    Z Basis Vector:  %10.3f  %10.3f  %10.3f a0\n", a2[0], a2[1], a2[2]);
    
    double *b0 = Rmg_L.b0;
    double *b1 = Rmg_L.b1;
    double *b2 = Rmg_L.b2;
    if(pct.imgpe==0) fprintf(ct.logfile, "    X Reci Vector:  %10.3f  %10.3f  %10.3f a0\n", b0[0], b0[1], b0[2]);
    if(pct.imgpe==0) fprintf(ct.logfile, "    Y Reci Vector:  %10.3f  %10.3f  %10.3f a0\n", b1[0], b1[1], b1[2]);
    if(pct.imgpe==0) fprintf(ct.logfile, "    Z Reci Vector:  %10.3f  %10.3f  %10.3f a0\n", b2[0], b2[1], b2[2]);
    
    if(pct.imgpe==0) fprintf(ct.logfile, "\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "K-points\n");
    if(ct.is_gamma)
    {
	if(pct.imgpe==0) fprintf(ct.logfile, "    Gamma Point Only (real orbitals)\n");
    }
    else
    {
	if(pct.imgpe==0) fprintf(ct.logfile, "    Brillouin Zone sampling with %d K-points (orbitals are complex)\n", ct.num_kpts);
	if(pct.imgpe==0) fprintf(ct.logfile, "\n");
	if(pct.imgpe==0) fprintf(ct.logfile, "         Kx      Ky        Kz     Weight in crystal unit\n");
	for (kpt = 0; kpt < ct.num_kpts; kpt++)
	{
	    if(pct.imgpe==0) fprintf(ct.logfile, "    %8.4f   %8.4f   %8.4f   %5.3f\n",
		    ct.kp[kpt].kpt[0], ct.kp[kpt].kpt[1], ct.kp[kpt].kpt[2], ct.kp[kpt].kweight);

	}

	if(pct.imgpe==0) fprintf(ct.logfile, "\n");
	if(pct.imgpe==0) fprintf(ct.logfile, "         Kx      Ky        Kz     Weight in 2PI/a \n");
        double kunit = twoPI /Rmg_L.celldm[0];
	for (kpt = 0; kpt < ct.num_kpts; kpt++)
	{
	    if(pct.imgpe==0) fprintf(ct.logfile, "    %8.4f   %8.4f   %8.4f   %5.3f\n",
		    ct.kp[kpt].kvec[0]/kunit, ct.kp[kpt].kvec[1]/kunit, ct.kp[kpt].kvec[2]/kunit, ct.kp[kpt].kweight);

	}
    }

    
    if(pct.imgpe==0) fprintf(ct.logfile, "\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "Atoms and States\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "    Number of ions:                          %lu\n", Atoms.size());
    if(pct.imgpe==0) fprintf(ct.logfile, "    Number of species:                       %lu\n", Species.size());
    if (ct.spin_flag)
    {
	if(pct.imgpe==0) fprintf(ct.logfile, "    Number of spin up states:                %d\n", ct.run_states);
	if(pct.imgpe==0) fprintf(ct.logfile, "    Number of spin down states:              %d\n", ct.run_states);
    }
    else
    {
	if(pct.imgpe==0) fprintf(ct.logfile, "    Number of states:                        %d\n", ct.run_states);
    }	
    
    
    



    if(pct.imgpe==0) fprintf(ct.logfile, "\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "Run Parameters\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "    SCF Convergence criterion (potential):   %4.2e\n", ct.thr_rms);
    if(pct.imgpe==0) fprintf(ct.logfile, "    SCF Convergence criterion (energy):      %4.2e\n", ct.thr_energy);
    if(pct.imgpe==0) fprintf(ct.logfile, "    Max SCF steps:                           %d\n", ct.max_scf_steps);
    if (ct.forceflag == MD_FASTRLX)
	if(pct.imgpe==0) fprintf(ct.logfile, "    Structural optimization force criterion: %d\n", ct.max_md_steps);
    if (ct.forceflag != MD_QUENCH)
    {
	if(pct.imgpe==0) fprintf(ct.logfile, "    Max MD steps                             %d\n", ct.max_md_steps);
        if(pct.imgpe==0) fprintf(ct.logfile, "    Timestep for molecular dynamics:         %12.8f\n", ct.iondt);
	if(pct.imgpe==0) fprintf(ct.logfile, "    Restart data write frequency:            %d MD steps\n", ct.checkpoint);
    }
    
    
    
    if(pct.imgpe==0) fprintf(ct.logfile, "\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "SCF Cycle Settings\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "    Charge density mixing:                   ");
    switch(ct.charge_mixing_type) {
        case 0:
            if(pct.imgpe==0) fprintf(ct.logfile, "Linear  (Mixing constant %4.2f)\n", ct.mix);
            break;
	case 1:
	    if(pct.imgpe==0) fprintf(ct.logfile, "Pulay  (Order:%d Scale:%4.2f Refresh:%d)\n", ct.charge_pulay_order, ct.mix, ct.charge_pulay_refresh);
	    break;
	case 2:
	    if(pct.imgpe==0) fprintf(ct.logfile, "Broyden\n");
	    break;
        default:
            if(pct.imgpe==0) fprintf(ct.logfile, "Unknown charge mixing method\n");
    }
    
    if(pct.imgpe==0) fprintf(ct.logfile, "    Hartree Solver:                          ");
    switch(ct.poisson_solver) {
        case POISSON_PFFT_SOLVER:
            if(pct.imgpe==0) fprintf(ct.logfile, "PFFT\n");
            break;
	case MULTIGRID_SOLVER:
	    if(pct.imgpe==0) fprintf(ct.logfile, "Multigrid\n");
	    break;
        default:
            if(pct.imgpe==0) fprintf(ct.logfile, "Unknown Hartree solver\n");
    }
    
    if(pct.imgpe==0) fprintf(ct.logfile, "    Interpolation type:                      ");
    switch(ct.interp_flag) {
        case BSPLINE_INTERPOLATION:
            if(pct.imgpe==0) fprintf(ct.logfile, "B-spline   (Order %d  using trade_image %d)\n",
                ct.interp_order, ct.interp_trade);
            break;
        case PROLONG_INTERPOLATION:
            if(pct.imgpe==0) fprintf(ct.logfile, "Prolong\n");
            break;
        case FFT_INTERPOLATION:
            if(pct.imgpe==0) fprintf(ct.logfile, "FFT\n");
            break;
        default:
            if(pct.imgpe==0) fprintf(ct.logfile, "Cubic\n");
    }

    if(pct.imgpe==0) fprintf(ct.logfile, "\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "Adaptive Interpolation\n");
    if(ct.cmix != 0.0)
    {
        if(pct.imgpe==0) fprintf(ct.logfile, "    On: cmix                                %12.9f\n", ct.cmix);
    }
    else
    {
        if(pct.imgpe==0) fprintf(ct.logfile, "    Off\n");
    }

    FiniteDiff FD(&Rmg_L);

    if(pct.imgpe==0) fprintf(ct.logfile, "\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "Finite Difference parameters\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "    Kohn Sham FD Order:                      %-4d\n", ct.kohn_sham_fd_order);
    if(ct.force_grad_order > 0)
    {
        if(pct.imgpe==0) fprintf(ct.logfile, "    Force Gradient Order:                    %-4d\n", ct.force_grad_order);
    }
    else
    {
        if(pct.imgpe==0) fprintf(ct.logfile, "    Force Gradient Order:                    FFT\n");
    }
    if(ct.alt_laplacian)
    {
        if(pct.imgpe==0) fprintf(ct.logfile, "    Adaptive Finite Differencing:            On\n");
        if(pct.imgpe==0) fprintf(ct.logfile, "    Adaptive CFAC:                          %12.9f\n", FD.cfac[0]);
    }
    else
    {
        if(pct.imgpe==0) fprintf(ct.logfile, "    Adaptive Finite Differencing:            Off\n");
    }
    
    if(pct.imgpe==0) fprintf(ct.logfile, "\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "Subspace Diagonalization Options\n");
    
    if(pct.imgpe==0) fprintf(ct.logfile, "    Frequency:                               every %d SCF step(s)\n", ct.diag);
    
    if(pct.imgpe==0) fprintf(ct.logfile, "    Driver:                                  ");
    switch(ct.subdiag_driver) {
        case SUBDIAG_SCALAPACK:
            if(pct.imgpe==0) fprintf(ct.logfile, "ScaLapack\n");
            break;
        case SUBDIAG_LAPACK:
#if CUDA_ENABLED
            if(pct.imgpe==0) fprintf(ct.logfile, "Lapack changed to MAGMA\n");
#else
            if(pct.imgpe==0) fprintf(ct.logfile, "Lapack\n");
#endif
            break;
        case SUBDIAG_MAGMA:
            if(pct.imgpe==0) fprintf(ct.logfile, "MAGMA\n");
            break;
        case SUBDIAG_CUSOLVER:
            if(pct.imgpe==0) fprintf(ct.logfile, "Cusolver\n");
            break;
        case SUBDIAG_ROCSOLVER:
            if(pct.imgpe==0) fprintf(ct.logfile, "Rocsolver\n");
            break;
        default:
            if(pct.imgpe==0) fprintf(ct.logfile, "Unknown diagonalization method");
    }
    
    if(pct.imgpe==0) fprintf(ct.logfile, "    Initial diagonalization:                 ");
    if (ct.initdiag)
    {
	if(pct.imgpe==0) fprintf(ct.logfile, "ON\n");
    }
    else
    {
	if(pct.imgpe==0) fprintf(ct.logfile, "OFF\n");
    }
    
    if(pct.imgpe==0) fprintf(ct.logfile, "    Folded spectrum:                         ");
    if (ct.use_folded_spectrum)
    {
	if(pct.imgpe==0) fprintf(ct.logfile, "ON\n");
    }
    else
    {
	if(pct.imgpe==0) fprintf(ct.logfile, "OFF\n");
    }
    

    if(pct.imgpe==0) fprintf(ct.logfile, "\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "Filtering Cutoff  Parameters\n");  
    if(pct.imgpe==0) fprintf(ct.logfile, "    Wavefunction grid (cparm):               %5.3f\n", ct.cparm);
    if(pct.imgpe==0) fprintf(ct.logfile, "    Charge density grid (rhocparm):          %5.3f\n", ct.rhocparm);
    if(pct.imgpe==0) fprintf(ct.logfile, "    Filter factor:                           %5.3f\n", ct.filter_factor);


    if(pct.imgpe==0) fprintf(ct.logfile, "\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "Multigrid (MG) Parameters\n");
    if (ct.poisson_solver == MULTIGRID_SOLVER) {
	if(pct.imgpe==0) fprintf(ct.logfile, "\n");
	if(pct.imgpe==0) fprintf(ct.logfile, "    Poisson MG levels:                   %d\n", ct.poi_parm.levels);
	if(pct.imgpe==0) fprintf(ct.logfile, "    Poisson global step:                 %-6.3f\n", ct.poi_parm.gl_step);
	if(pct.imgpe==0) fprintf(ct.logfile, "    Poisson pre:                         %d\n", ct.poi_parm.gl_pre);
	if(pct.imgpe==0) fprintf(ct.logfile, "    Poisson post:                        %d\n", ct.poi_parm.gl_pst);
    }

    if(pct.imgpe==0) fprintf(ct.logfile, "    Psi MG levels:                           %d\n", ct.eig_parm.levels);
    if(pct.imgpe==0) fprintf(ct.logfile, "    Psi global step:                         %-6.3f\n", ct.eig_parm.gl_step);
    if(pct.imgpe==0) fprintf(ct.logfile, "    Psi pre:                                 %d\n", ct.eig_parm.gl_pre);
    if(pct.imgpe==0) fprintf(ct.logfile, "    Psi post:                                %d\n", ct.eig_parm.gl_pst);

    if(pct.imgpe==0) fprintf(ct.logfile, "\n");
    if (ct.kohn_sham_solver == DAVIDSON_SOLVER) {
        if(pct.imgpe==0) fprintf(ct.logfile, "Davidson Parameters\n");
	if(pct.imgpe==0) fprintf(ct.logfile, "    Davidson multiplier:                     %d\n", ct.davidx);
	if(pct.imgpe==0) fprintf(ct.logfile, "    Davidson max step:                       %d\n", ct.david_max_steps);
	if(pct.imgpe==0) fprintf(ct.logfile, "    Davidson unocc tol factor:               %-6.3f\n", ct.unoccupied_tol_factor);
    }

    if(pct.imgpe==0) fprintf(ct.logfile, "\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "Blas Libraries\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "    CPU support with %s\n", RMG_BLAS_MSG);
#if HIP_ENABLED
    if(pct.imgpe==0) fprintf(ct.logfile, "    GPU support with hipblas\n");
#endif
#if CUDA_ENABLED
    if(pct.imgpe==0) fprintf(ct.logfile, "    GPU support with cublas\n");
#endif
    std::string serial("serial"), unknown("unknown"), msg(RMG_BLAS_MSG);
    size_t found = msg.find(serial);
    if (found!=std::string::npos)
    {
        if(pct.imgpe==0) fprintf(ct.logfile, "    WARNING!!!  serial blas library detected. Performance in hybrid mode may be suboptimal\n    with a CPU only build.");
        std::cout << "    WARNING!!!  serial blas library detected. Performance in hybrid mode may be suboptimal\n    with a CPU only build." << std::endl;

    }
    found = msg.find(unknown);
    if (found!=std::string::npos)
    {
        if(pct.imgpe==0) fprintf(ct.logfile, "    WARNING!!!  unknown blas library detected. Performance may be suboptimal\n    with a CPU only build.");
        std::cout << "WARNING!!!  unknown blas library detected. Performance may be suboptimal\n    with a CPU only build." << std::endl;
    }

#if 0
    /* Forces are updated under normalized constraint field */
    if (verify ("atom_constraints", NULL))
    {
        if(pct.imgpe==0) fprintf(ct.logfile, "    Constrained per atom dynamics vector field.\n");
        for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
        {
            if(pct.imgpe==0) fprintf(ct.logfile, "       % 10f % 10f % 10f % 10f\n",
					Atoms[ion].constraint.setA_coord[0],
					Atoms[ion].constraint.setA_coord[1],
					Atoms[ion].constraint.setA_coord[2],
                    Atoms[ion].constraint.setA_weight);
        }
    }
#endif
    
    /**********  Begin Species Table  ************/
    if(pct.imgpe==0) fprintf(ct.logfile, "\n");


    /*Determine spacing for column about functional, string is 24 characters long, but usually is shorter
     * Having excessive amount of space look strange, so we only use as much as is needed*/
    funct_legend_length = strlen("Functional");
    max_funct_length = 0;
    for (auto& sp : Species)
    {
	if ((int)sp.functional.length() > max_funct_length)
	    max_funct_length = sp.functional.length();
    }


    /* The column cannot be shorter than its legend*/
    funct_spacing = max_funct_length;
    if ( funct_spacing < funct_legend_length) funct_spacing = funct_legend_length;

    funct_padding_left = (funct_spacing - funct_legend_length) / 2;
    funct_padding_right = funct_spacing - funct_padding_left - funct_legend_length;

    


    if(pct.imgpe==0) fprintf(ct.logfile, "\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "Atomic Species Information\n(PP = Pseudopotential, US = Ultrasoft, NC = Norm Conserving)\n");
    
    /*Begin table printout, vertical line first*/
    if(pct.imgpe==0) fprintf(ct.logfile, "---------------------------------------------------------------");
    for (i=0; i < funct_spacing; i++)
	if(pct.imgpe==0) fprintf(ct.logfile, "-");
    if(pct.imgpe==0) fprintf(ct.logfile, "\n");
    
    /*Table legend, line 1*/
    if(pct.imgpe==0) fprintf(ct.logfile, "|Index|Symbol| Mass|Valence| PP |  Comp  |Local| Local|Nlocal|");
     
    for (i=0; i < funct_padding_left + 4; i++)
	if(pct.imgpe==0) fprintf(ct.logfile, " ");
    if(pct.imgpe==0) fprintf(ct.logfile, "PP");
    for (i=0; i < funct_padding_right + 4; i++)
	if(pct.imgpe==0) fprintf(ct.logfile, " ");
    if(pct.imgpe==0) fprintf(ct.logfile, "|\n");
    
    
    /*Table legend, line 2*/
    if(pct.imgpe==0) fprintf(ct.logfile, "|     |      |     |       |Type|Gaussian|  l  |Radius|Radius|");
    
    for (i=0; i < funct_padding_left; i++)
	if(pct.imgpe==0) fprintf(ct.logfile, " ");
    if(pct.imgpe==0) fprintf(ct.logfile, "Functional");
    for (i=0; i < funct_padding_right; i++)
	if(pct.imgpe==0) fprintf(ct.logfile, " ");
    if(pct.imgpe==0) fprintf(ct.logfile, "|\n");

    if(pct.imgpe==0) fprintf(ct.logfile, "---------------------------------------------------------------");
    for (i=0; i < funct_spacing; i++)
	if(pct.imgpe==0) fprintf(ct.logfile, "-");
    if(pct.imgpe==0) fprintf(ct.logfile, "\n");


    for (auto& sp : Species)
    {
	if(pct.imgpe==0) fprintf(ct.logfile, "|%5d",sp.index + 1);
	if(pct.imgpe==0) fprintf(ct.logfile, "|%6.6s", sp.atomic_symbol);
	if(pct.imgpe==0) fprintf(ct.logfile, "|%5.1lf", sp.atomic_mass);
	if(pct.imgpe==0) fprintf(ct.logfile, "|%7.2lf", sp.zvalence);
	if (sp.is_norm_conserving)
        {
	    if(pct.imgpe==0) fprintf(ct.logfile, "|  NC");     
        }
	else
        {
	    if(pct.imgpe==0) fprintf(ct.logfile, "|  US");     
        }
	if(pct.imgpe==0) fprintf(ct.logfile, "|%8.2lf", sp.rc);
	if(pct.imgpe==0) fprintf(ct.logfile, "|%5d", sp.local);
	if(pct.imgpe==0) fprintf(ct.logfile, "|%6.2lf", sp.lradius); 
	if(pct.imgpe==0) fprintf(ct.logfile, "|%6.2lf", sp.nlradius); 
	if(pct.imgpe==0) fprintf(ct.logfile, "|%*s", funct_spacing,sp.functional.c_str());
	if(pct.imgpe==0) fprintf(ct.logfile, "|\n");
    }
    
    if(pct.imgpe==0) fprintf(ct.logfile, "---------------------------------------------------------------");
    for (i=0; i < funct_spacing; i++)
	if(pct.imgpe==0) fprintf(ct.logfile, "-");
    if(pct.imgpe==0) fprintf(ct.logfile, "\n");
    /**********  End Species Table  ************/
    
    if(pct.imgpe==0) fprintf(ct.logfile, "\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "Pseudopotential generation information:\n");
    for (auto& sp : Species)
    {
        if(pct.imgpe==0) fprintf(ct.logfile, "  %2s pseudopotential file: %s\n", sp.atomic_symbol, sp.pseudo_filename.c_str());
        if(pct.imgpe==0) fprintf(ct.logfile, "      Generation info     : %s\n", sp.generated.c_str());
        if(pct.imgpe==0) fprintf(ct.logfile, "      Author info         : %s\n", sp.author.c_str());
    }

    if(pct.imgpe==0) fprintf(ct.logfile, "\n\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "Memory usage (Mbytes):     Min        Max       Total\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "    wave functions      %8.2f   %8.2f   %8.2f\n",
                            (double)ct.psi_alloc[1] / 1000.0 / 1000.0,
                            (double)ct.psi_alloc[2] / 1000.0 / 1000.0,
                            (double)ct.psi_alloc[0] / 1000.0 / 1000.0);
    if(pct.imgpe==0) fprintf(ct.logfile, "    beta-functions      %8.2f   %8.2f   %8.2f\n",
                            (double)ct.beta_alloc[1] / 1000.0 / 1000.0,
                            (double)ct.beta_alloc[2] / 1000.0 / 1000.0,
                            (double)ct.beta_alloc[0] / 1000.0 / 1000.0);
    if(!ct.norm_conserving_pp)
    {
        if(pct.imgpe==0) fprintf(ct.logfile, "    q-functions         %8.2f   %8.2f   %8.2f\n",
                                (double)ct.q_alloc[1] / 1000.0 / 1000.0,
                                (double)ct.q_alloc[2] / 1000.0 / 1000.0,
                                (double)ct.q_alloc[0] / 1000.0 / 1000.0);
    }
    if(ct.xc_is_hybrid)
    {
        if(pct.imgpe==0) fprintf(ct.logfile, "    vexx                %8.2f   %8.2f   %8.2f\n",
                                (double)ct.vexx_alloc[1] / 1000.0 / 1000.0,
                                (double)ct.vexx_alloc[2] / 1000.0 / 1000.0,
                                (double)ct.vexx_alloc[0] / 1000.0 / 1000.0);
    }

    if(pct.imgpe==0) fprintf(ct.logfile, "\n");

    /* Write out the ionic positions and displacements */
    init_write_pos ();

#if 0
    if ((pct.imgpe == 0) && (verify ("pdb_atoms", NULL)))
        write_pdb ();
#endif

    crho_fract = ct.crho - ct.ionic_charge;
    if((fabs(crho_fract) > 1.0e-08) && ct.localize_localpp) {
        if(pct.imgpe==0) fprintf(ct.logfile, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
        if(pct.imgpe==0) fprintf(ct.logfile, "    crho %e  %e", crho_fract, ct.crho - ct.nel);
         
        if(pct.imgpe==0) fprintf(ct.logfile, "    WARNING: FRACTIONAL PART OF COMPENSATING CHARGES IS LARGER THAN TOLERANCE!!!\n");
        if(pct.imgpe==0) fprintf(ct.logfile, "    THIS WILL SET A LIMIT ON THE CONVERGENCE OF THE HARTREE POTENTIAL!!!\n");
        if(pct.imgpe==0) fprintf(ct.logfile, "    THIS CAN USUALLY BE CORRECTED BY INCREASING THE RADII IN THE PSEUDOPOTENTIAL FILES.\n");
        if(pct.imgpe==0) fprintf(ct.logfile, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    }



}                               /* end write_header */


/********************************/









static void init_write_pos (void)
{

    if(pct.imgpe==0) fprintf(ct.logfile, "\n\nInitial Ionic Positions And Displacements (Bohr)\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "Species      X           Y           Z           dX          dY          dZ");

    for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
    {

        ION &Atom = Atoms[ion];
        SPECIES &Type = Species[Atom.species];

        double crds[3], icrds[3];
        Rmg_L.to_cartesian_input (Atom.xtal, crds);
        Rmg_L.to_cartesian_input (Atom.ixtal, icrds);

        if(pct.imgpe==0) fprintf(ct.logfile, "\n  %-2s   %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f",
		Type.atomic_symbol,
                crds[0], crds[1], crds[2],
                crds[0] - icrds[0],
                crds[1] - icrds[1], 
                crds[2] - icrds[2]);

    }                           /* end for */

    if(pct.imgpe==0) fprintf(ct.logfile, "\n");


    if(pct.imgpe==0) fprintf(ct.logfile, "\n\nInitial Ionic Positions And Displacements (Angstrom)\n");
    if(pct.imgpe==0) fprintf(ct.logfile, "Species      X           Y           Z           dX          dY          dZ");

    for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
    {

        ION &Atom = Atoms[ion];
        SPECIES &Type = Species[Atom.species];

        double crds[3], icrds[3];
        Rmg_L.to_cartesian_input (Atom.xtal, crds);
        Rmg_L.to_cartesian_input (Atom.ixtal, icrds);

        if(pct.imgpe==0) fprintf(ct.logfile, "\n  %-2s   %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f",
		Type.atomic_symbol,
                crds[0] *a0_A, crds[1] *a0_A, crds[2] *a0_A,
                (crds[0] - icrds[0]) *a0_A,
                (crds[1] - icrds[1]) *a0_A, 
                (crds[2] - icrds[2]) *a0_A);

    }                           /* end for */

    if(pct.imgpe==0) fprintf(ct.logfile, "\n");

}                               /* end write_pos */

