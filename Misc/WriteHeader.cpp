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


    printf ("\n\n");
    printf ("                     * * * * * * * * * *\n");
    printf ("                     *    R   M   G    *\n");
    printf ("                     * * * * * * * * * *\n");
    printf ("\n");
    printf ("     -- A Real Space Multigrid Electronic structure code --\n");
    printf ("     --      More information at www.rmgdft.org          --\n");
        

    printf ("\nCode Revision:     %s", RMG_REVISION);
    printf ("\nBuild Date:        %s  %s", __DATE__, __TIME__);
    printf ("\nStart time:        %s", timeptr);

    printf("\nNOTICE: RMG internal pseudopotentials have switched to");
    printf("\nONCVP from Ultrasoft. You can revert to Ultrasoft by");
    printf("\nadding the input tag internal_pseudo_type=\"ultrasoft\" to");
    printf("\nyour input files.\n\n");

    printf ("\n");
    printf ("Files\n");
    printf ("   Control input file:        %s\n", ct.cfile);
    if (ct.runflag == 1)
	printf ("   Data input file:           %s\n", ct.infile);
    printf ("   Data output file:          %s\n", ct.outfile);
    
    printf ("\n");
    printf ("Run Setup\n");
    printf ("    Calculation type:         ");
    switch (ct.forceflag)
    {
    case MD_QUENCH:
        printf ("Quench electrons - Fixed ions SCF calculation\n");
        break;

    case MD_FASTRLX:
        printf ("Structure Optimization.\n");
        break;

    case BAND_STRUCTURE:
        printf ("Band structure calculation.\n");
        break;

    case MD_CVE:
        printf ("Molecular dynamics - CVE\n");
        break;

    case MD_CVT:
        printf ("Molecular dynamics - CVT\n");
        break;

    case MD_CPT:
        printf ("Molecular dynamics - CPT\n");
        break;

    case PLOT:
        printf ("Plot density in DX form.\n");
        break;

    case PSIPLOT:
        printf ("Plot Psi^2 in DX form.\n");
        break;

    case NEB_RELAX:
        printf ("Molecular dynamics using Nudged Elastic Band.\n");
        break;

    case TDDFT:
        printf ("Time dependent DFT (TDDFT) calculation \n");
        break;
    case Exx_only:
        printf ("calculate Exx integral's from saveed wave functions \n");
        break;
    case STM:
        printf ("calculate STM charge density \n");
        break;


    default:
        error_handler ("Unknown molecular dynamics method.");
    }
    printf ("    Description:              %s\n", ct.description.c_str());
    printf ("    Orbital Initialization:   ");
    switch (ct.runflag)
    {

	case RESTART:
	    printf ("Read from %s \n", ct.infile);
	    break;

	case RANDOM_START:
	    printf ("Random\n");
	    break;
	
	case LCAO_START:
	    printf ("LCAO (%d LCAO and %d random orbitals)\n",  ct.total_atomic_orbitals, ct.extra_random_lcao_states);
	    break;

	case MODIFIED_LCAO_START:
	    printf ("LCAO (%d MODIFIED LCAO and %d random orbitals)\n",  ct.init_states, ct.extra_random_lcao_states);
	    break;
	
	case INIT_FIREBALL:
	    printf ("Fireball\n");
	    break;
	
	case INIT_GAUSSIAN:
	    printf ("Gaussian\n");
	    break;
	
	case Start_TDDFT:
	    printf ("Initial start TD DFT calculation\n");
	    break;
	
	case Restart_TDDFT:
	    printf ("Restart TD DFT calculation\n");
	    break;
        
	default:
            printf ("Unknown start mode\n");
    }
    const std::string tstr = Functional::get_dft_name_rmg();
    printf("    Exchange Correlation:     %s\n", tstr.c_str());
    printf("    Spin Polarization:        ");
    if (ct.spin_flag)
	printf("ON\n");
    else
	printf("OFF\n");
    if (fabs(ct.background_charge) > 1e-6)
	printf("    System charge:            %6.2f\n", ct.background_charge);
    else
	printf("    System charge:            Neutral\n");


    
    printf ("\n");
    printf ("Processor Topology\n");  
    printf ("   Total PEs:                 %d\n", (get_PE_X() * get_PE_Y() * get_PE_Z()));
    printf ("   X direction:               %d\n", get_PE_X());
    printf ("   Y direction:               %d\n", get_PE_Y());
    printf ("   Z direction:               %d\n", get_PE_Z());
    printf ("   MG Threads/PE:             %d\n", ct.MG_THREADS_PER_NODE);
    printf ("   OMP Threads/PE:            %d\n", ct.OMP_THREADS_PER_NODE);
    
    
    
    printf ("\n");
    printf ("Grid Points");
    if (fabs (1.0 - get_anisotropy()) > 0.005)
	printf("  (Anisotropy: %5.3f)", get_anisotropy());
    printf ("\n");
    
    printf ("    X:  Total: %d   Per PE: %d   Spacing:%5.3f a0  \n", get_NX_GRID(), get_PX0_GRID(),   get_hxgrid() * get_xside());
    printf ("    Y:  Total: %d   Per PE: %d   Spacing:%5.3f a0  \n", get_NY_GRID(), get_PY0_GRID(),   get_hygrid() * get_yside());
    printf ("    Z:  Total: %d   Per PE: %d   Spacing:%5.3f a0  \n", get_NZ_GRID(), get_PZ0_GRID(),   get_hzgrid() * get_zside());
    printf ("\n");

    if(ct.coalesce_states)
    {
        printf ("    Coalescing states in X with factor %d\n", pct.coalesce_factor);
    }

    /* We compute the equivalent energy cutoff using the density of grid
     * points in the cell with a correction for the grid anisotropy.
     */
    double tpiba2 = 4.0 * PI * PI / (Rmg_L.celldm[0] * Rmg_L.celldm[0]);
    double t1 = coarse_pwaves->gcut * tpiba2;
    double t2 = fine_pwaves->gcut * tpiba2;
    int ibrav = get_ibrav_type();
    if(ibrav < 0) ibrav *= -1;
    printf ("    Equivalent energy cutoffs (psi,rho):  %8.3f   %8.3f Ry\n", t1, t2);
    printf ("\n");
    printf ("    Charge density grid:         %d times finer\n", get_FG_RATIO());


    printf ("\n");
    printf ("\n");
    printf ("Lattice (Cell) Setup\n");
    printf ("    Type:                       %s\n", lattice_type[ibrav]);
    printf ("    Volume (a0^3):              %8.2f\n", get_vel() * ct.psi_nbasis);
    printf ("    Boundary conditions:        ");
    switch (ct.boundaryflag)
    {

    case PERIODIC:
        printf ("Periodic\n");
        break;

    case CLUSTER:
        printf ("Cluster\n");
        break;

    case SURFACE:
        printf ("Surface\n");
        break;
	
    default:
	printf ("Unknown boundary conditions\n");

    }                           /* end switch */
    printf ("\n");

    // Print out original lattice vectors
    double *a0 = Rmg_L.a0i;
    double *a1 = Rmg_L.a1i;
    double *a2 = Rmg_L.a2i;
    printf ("    X Basis Vector:  %10.3f  %10.3f  %10.3f a0\n", a0[0], a0[1], a0[2]);
    printf ("    Y Basis Vector:  %10.3f  %10.3f  %10.3f a0\n", a1[0], a1[1], a1[2]);
    printf ("    Z Basis Vector:  %10.3f  %10.3f  %10.3f a0\n", a2[0], a2[1], a2[2]);
    
    double *b0 = Rmg_L.b0;
    double *b1 = Rmg_L.b1;
    double *b2 = Rmg_L.b2;
    printf ("    X Reci Vector:  %10.3f  %10.3f  %10.3f a0\n", b0[0], b0[1], b0[2]);
    printf ("    Y Reci Vector:  %10.3f  %10.3f  %10.3f a0\n", b1[0], b1[1], b1[2]);
    printf ("    Z Reci Vector:  %10.3f  %10.3f  %10.3f a0\n", b2[0], b2[1], b2[2]);
    
    printf ("\n");
    printf ("K-points\n");
    if(ct.is_gamma)
    {
	printf ("    Gamma Point Only (real orbitals)\n");
    }
    else
    {
	printf ("    Brillouin Zone sampling with %d K-points (orbitals are complex)\n", ct.num_kpts);
	printf ("\n");
	printf ("         Kx      Ky        Kz     Weight in crystal unit\n");
	for (kpt = 0; kpt < ct.num_kpts; kpt++)
	{
	    printf ("    %8.4f   %8.4f   %8.4f   %5.3f\n",
		    ct.kp[kpt].kpt[0], ct.kp[kpt].kpt[1], ct.kp[kpt].kpt[2], ct.kp[kpt].kweight);

	}

	printf ("\n");
	printf ("         Kx      Ky        Kz     Weight in 2PI/a \n");
        double kunit = twoPI /Rmg_L.celldm[0];
	for (kpt = 0; kpt < ct.num_kpts; kpt++)
	{
	    printf ("    %8.4f   %8.4f   %8.4f   %5.3f\n",
		    ct.kp[kpt].kvec[0]/kunit, ct.kp[kpt].kvec[1]/kunit, ct.kp[kpt].kvec[2]/kunit, ct.kp[kpt].kweight);

	}
    }

    
    printf ("\n");
    printf ("Atoms and States\n");
    printf ("    Number of ions:                          %lu\n", Atoms.size());
    printf ("    Number of species:                       %lu\n", Species.size());
    if (ct.spin_flag)
    {
	printf ("    Number of spin up states:                %d\n", ct.run_states);
	printf ("    Number of spin down states:              %d\n", ct.run_states);
    }
    else
    {
	printf ("    Number of states:                        %d\n", ct.run_states);
    }	
    
    
    



    printf ("\n");
    printf ("Run Parameters\n");
    printf ("    SCF Convergence criterion (potential):   %4.2e\n", ct.thr_rms);
    printf ("    SCF Convergence criterion (energy):      %4.2e\n", ct.thr_energy);
    printf ("    Max SCF steps:                           %d\n", ct.max_scf_steps);
    if (ct.forceflag == MD_FASTRLX)
	printf ("    Structural optimization force criterion: %d\n", ct.max_md_steps);
    if (ct.forceflag != MD_QUENCH)
    {
	printf ("    Max MD steps                             %d\n", ct.max_md_steps);
        printf ("    Timestep for molecular dynamics:         %12.8f\n", ct.iondt);
	printf ("    Restart data write frequency:            %d MD steps\n", ct.checkpoint);
    }
    
    
    
    printf ("\n");
    printf ("SCF Cycle Settings\n");
    printf ("    Charge density mixing:                   ");
    switch(ct.charge_mixing_type) {
        case 0:
            printf ("Linear  (Mixing constant %4.2f)\n", ct.mix);
            break;
	case 1:
	    printf ("Pulay  (Order:%d Scale:%4.2f Refresh:%d)\n", ct.charge_pulay_order, ct.charge_pulay_scale, ct.charge_pulay_refresh);
	    break;
	case 2:
	    printf ("Broyden\n");
	    break;
        default:
            printf ("Unknown charge mixing method\n");
    }
    
    printf ("    Hartree Solver:                          ");
    switch(ct.poisson_solver) {
        case POISSON_PFFT_SOLVER:
            printf ("PFFT\n");
            break;
	case MULTIGRID_SOLVER:
	    printf ("Multigrid\n");
	    break;
        default:
            printf ("Unknown Hartree solver\n");
    }
    
    printf ("    Interpolation type:                      ");
    switch(ct.interp_flag) {
        case BSPLINE_INTERPOLATION:
            printf ("B-spline   (Order %d  using trade_image %d)\n",
                ct.interp_order, ct.interp_trade);
            break;
        case PROLONG_INTERPOLATION:
            printf ("Prolong\n");
            break;
        case FFT_INTERPOLATION:
            printf ("FFT\n");
            break;
        default:
            printf ("Cubic\n");
    }

    
    printf ("\n");
    printf ("Subspace Diagonalization Options\n");
    
    printf ("    Frequency:                               every %d SCF step(s)\n", ct.diag);
    
    printf ("    Driver:                                  ");
    switch(ct.subdiag_driver) {
        case SUBDIAG_SCALAPACK:
            printf ("ScaLapack\n");
            break;
        case SUBDIAG_LAPACK:
#if CUDA_ENABLED
            printf ("Lapack changed to MAGMA\n");
#else
            printf ("Lapack\n");
#endif
            break;
        case SUBDIAG_MAGMA:
            printf ("MAGMA\n");
            break;
        case SUBDIAG_CUSOLVER:
            printf ("Cusolver\n");
            break;
        case SUBDIAG_ROCSOLVER:
            printf ("Rocsolver\n");
            break;
        default:
            printf ("Unknown diagonalization method");
    }
    
    printf ("    Initial diagonalization:                 ");
    if (ct.initdiag)
	printf("ON\n");
    else
	printf("OFF\n");
    
    printf ("    Folded spectrum:                         ");
    if (ct.use_folded_spectrum)
	printf("ON\n");
    else
	printf("OFF\n");
    

    printf ("\n");
    printf ("Filtering Cutoff  Parameters\n");  
    printf ("    Wavefunction grid (cparm):               %5.3f\n", ct.cparm);
    printf ("    Charge density grid (rhocparm):          %5.3f\n", ct.rhocparm);
    printf ("    Filter factor:                           %5.3f\n", ct.filter_factor);


    printf ("\n");
    printf ("Multigrid (MG) Parameters\n");
    if (ct.poisson_solver == MULTIGRID_SOLVER) {
	printf ("\n");
	printf ("    Poisson MG levels:                   %d\n", ct.poi_parm.levels);
	printf ("    Poisson global step:                 %-6.3f\n", ct.poi_parm.gl_step);
	printf ("    Poisson pre:                         %d\n", ct.poi_parm.gl_pre);
	printf ("    Poisson post:                        %d\n", ct.poi_parm.gl_pst);
    }

    printf ("    Psi MG levels:                           %d\n", ct.eig_parm.levels);
    printf ("    Psi global step:                         %-6.3f\n", ct.eig_parm.gl_step);
    printf ("    Psi pre:                                 %d\n", ct.eig_parm.gl_pre);
    printf ("    Psi post:                                %d\n", ct.eig_parm.gl_pst);

    printf ("\n");
    if (ct.kohn_sham_solver == DAVIDSON_SOLVER) {
        printf ("Davidson Parameters\n");
	printf ("    Davidson multiplier:                     %d\n", ct.davidx);
	printf ("    Davidson max step:                       %d\n", ct.david_max_steps);
	printf ("    Davidson unocc tol factor:               %-6.3f\n", ct.unoccupied_tol_factor);
    }

    printf ("\n");
    printf ("Blas Libraries\n");
    printf ("    CPU support with %s\n", RMG_BLAS_MSG);
#if HIP_ENABLED
    printf ("    GPU support with hipblas\n");
#endif
#if CUDA_ENABLED
    printf ("    GPU support with cublas\n");
#endif
    std::string serial("serial"), unknown("unknown"), msg(RMG_BLAS_MSG);
    size_t found = msg.find(serial);
    if (found!=std::string::npos)
    {
        printf("    WARNING!!!  serial blas library detected. Performance in hybrid mode may be suboptimal\n    with a CPU only build.");
        std::cout << "    WARNING!!!  serial blas library detected. Performance in hybrid mode may be suboptimal\n    with a CPU only build." << std::endl;

    }
    found = msg.find(unknown);
    if (found!=std::string::npos)
    {
        printf("    WARNING!!!  unknown blas library detected. Performance may be suboptimal\n    with a CPU only build.");
        std::cout << "WARNING!!!  unknown blas library detected. Performance may be suboptimal\n    with a CPU only build." << std::endl;
    }

#if 0
    /* Forces are updated under normalized constraint field */
    if (verify ("atom_constraints", NULL))
    {
        printf ("    Constrained per atom dynamics vector field.\n");
        for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
        {
            printf ("       % 10f % 10f % 10f % 10f\n",
					Atoms[ion].constraint.setA_coord[0],
					Atoms[ion].constraint.setA_coord[1],
					Atoms[ion].constraint.setA_coord[2],
                    Atoms[ion].constraint.setA_weight);
        }
    }
#endif
    
    /**********  Begin Species Table  ************/
    printf ("\n");


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

    


    printf ("\n");
    printf ("Atomic Species Information\n(PP = Pseudopotential, US = Ultrasoft, NC = Norm Conserving)\n");
    
    /*Begin table printout, vertical line first*/
    printf ("---------------------------------------------------------------");
    for (i=0; i < funct_spacing; i++)
	printf ("-");
    printf("\n");
    
    /*Table legend, line 1*/
    printf ("|Index|Symbol| Mass|Valence| PP |  Comp  |Local| Local|Nlocal|");
     
    for (i=0; i < funct_padding_left + 4; i++)
	printf(" ");
    printf("PP");
    for (i=0; i < funct_padding_right + 4; i++)
	printf(" ");
    printf("|\n");
    
    
    /*Table legend, line 2*/
    printf ("|     |      |     |       |Type|Gaussian|  l  |Radius|Radius|");
    
    for (i=0; i < funct_padding_left; i++)
	printf(" ");
    printf("Functional");
    for (i=0; i < funct_padding_right; i++)
	printf(" ");
    printf("|\n");

    printf ("---------------------------------------------------------------");
    for (i=0; i < funct_spacing; i++)
	printf ("-");
    printf("\n");


    for (auto& sp : Species)
    {
	printf("|%5d",sp.index + 1);
	printf("|%6.6s", sp.atomic_symbol);
	printf("|%5.1lf", sp.atomic_mass);
	printf("|%7.2lf", sp.zvalence);
	if (sp.is_norm_conserving)
	    printf("|  NC");     
	else
	    printf("|  US");     
	printf("|%8.2lf", sp.rc);
	printf("|%5d", sp.local);
	printf("|%6.2lf", sp.lradius); 
	printf("|%6.2lf", sp.nlradius); 
	printf("|%*s", funct_spacing,sp.functional.c_str());
	printf ("|\n");
    }
    
    printf ("---------------------------------------------------------------");
    for (i=0; i < funct_spacing; i++)
	printf ("-");
    printf("\n");
    /**********  End Species Table  ************/
    
    printf("\n");
    printf("Pseudopotential generation information:\n");
    for (auto& sp : Species)
    {
        printf("  %2s pseudopotential file: %s\n", sp.atomic_symbol, sp.pseudo_filename.c_str());
        printf("      Generation info     : %s\n", sp.generated.c_str());
        printf("      Author info         : %s\n", sp.author.c_str());
    }

    printf("\n\n");
    printf("Memory usage (Mbytes):     Min        Max       Total\n");
    printf("    wave functions      %8.2f   %8.2f   %8.2f\n",
                            (double)ct.psi_alloc[1] / 1000.0 / 1000.0,
                            (double)ct.psi_alloc[2] / 1000.0 / 1000.0,
                            (double)ct.psi_alloc[0] / 1000.0 / 1000.0);
    printf("    beta-functions      %8.2f   %8.2f   %8.2f\n",
                            (double)ct.beta_alloc[1] / 1000.0 / 1000.0,
                            (double)ct.beta_alloc[2] / 1000.0 / 1000.0,
                            (double)ct.beta_alloc[0] / 1000.0 / 1000.0);
    if(!ct.norm_conserving_pp)
    {
        printf("    q-functions         %8.2f   %8.2f   %8.2f\n",
                                (double)ct.q_alloc[1] / 1000.0 / 1000.0,
                                (double)ct.q_alloc[2] / 1000.0 / 1000.0,
                                (double)ct.q_alloc[0] / 1000.0 / 1000.0);
    }
    if(ct.xc_is_hybrid)
    {
        printf("    vexx                %8.2f   %8.2f   %8.2f\n",
                                (double)ct.vexx_alloc[1] / 1000.0 / 1000.0,
                                (double)ct.vexx_alloc[2] / 1000.0 / 1000.0,
                                (double)ct.vexx_alloc[0] / 1000.0 / 1000.0);
    }

    printf("\n");

    /* Write out the ionic positions and displacements */
    init_write_pos ();

#if 0
    if ((pct.imgpe == 0) && (verify ("pdb_atoms", NULL)))
        write_pdb ();
#endif

    crho_fract = ct.crho - ct.ionic_charge;
    if((fabs(crho_fract) > 1.0e-08) && ct.localize_localpp) {
        printf ("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
        printf ("    crho %e  %e", crho_fract, ct.crho - ct.nel);
         
        printf ("    WARNING: FRACTIONAL PART OF COMPENSATING CHARGES IS LARGER THAN TOLERANCE!!!\n");
        printf ("    THIS WILL SET A LIMIT ON THE CONVERGENCE OF THE HARTREE POTENTIAL!!!\n");
        printf ("    THIS CAN USUALLY BE CORRECTED BY INCREASING THE RADII IN THE PSEUDOPOTENTIAL FILES.\n");
        printf ("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    }



}                               /* end write_header */


/********************************/









static void init_write_pos (void)
{

    printf ("\n\nInitial Ionic Positions And Displacements (Bohr)\n");
    printf ("Species      X           Y           Z           dX          dY          dZ");

    for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
    {

        ION &Atom = Atoms[ion];
        SPECIES &Type = Species[Atom.species];

        double crds[3], icrds[3];
        Rmg_L.to_cartesian_input (Atom.xtal, crds);
        Rmg_L.to_cartesian_input (Atom.ixtal, icrds);

        printf ("\n  %-2s   %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f",
		Type.atomic_symbol,
                crds[0], crds[1], crds[2],
                crds[0] - icrds[0],
                crds[1] - icrds[1], 
                crds[2] - icrds[2]);

    }                           /* end for */

    printf ("\n");


    printf ("\n\nInitial Ionic Positions And Displacements (Angstrom)\n");
    printf ("Species      X           Y           Z           dX          dY          dZ");

    for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
    {

        ION &Atom = Atoms[ion];
        SPECIES &Type = Species[Atom.species];

        double crds[3], icrds[3];
        Rmg_L.to_cartesian_input (Atom.xtal, crds);
        Rmg_L.to_cartesian_input (Atom.ixtal, icrds);

        printf ("\n  %-2s   %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f",
		Type.atomic_symbol,
                crds[0] *a0_A, crds[1] *a0_A, crds[2] *a0_A,
                (crds[0] - icrds[0]) *a0_A,
                (crds[1] - icrds[1]) *a0_A, 
                (crds[2] - icrds[2]) *a0_A);

    }                           /* end for */

    printf ("\n");

}                               /* end write_pos */

