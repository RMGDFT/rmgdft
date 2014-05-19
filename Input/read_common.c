/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

/****f* QMD-MGDFT/read_control.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 * COPYRIGHT
 *   Copyright (C) 2005  Frisco Rose
 *                       Jerzy Bernholc
 *   
 * FUNCTION
 *	void read_common(CONTROL *c)
 * 		read information shared by RMG, ON, and NEGF.
 *    
 * INPUTS
 *   main control structure
 * OUTPUT
 *   variables in structure CONTROL c are updated
 *   in most of other file, the name is ct.... 
 *   see main.h for structure CONTROL
 * PARENTS
 *   main.c
 * CHILDREN
 * 
 * SOURCE */

#include <string.h>
#include <ctype.h>
#include <math.h>
#include "grid.h"
#include "main.h"
#include "common_prototypes.h"


void read_common ()
{
    int tmp, is, ibrav, flag;
    char *tbuf, *tptr;
    rmg_double_t celldm[6], a0[3], a1[3], a2[3], omega;
    static int run_count = - 1;
    int NX_GRID, NY_GRID, NZ_GRID;
    int FNX_GRID, FNY_GRID, FNZ_GRID;
    int PE_X, PE_Y, PE_Z;
    int FG_RATIO;

    run_count ++;
    

    my_malloc (tptr, MAX_PATH, char);
    tbuf = tptr;

    get_data ("processor_grid", tbuf, STR, NULL);
    PE_X = strtol(tbuf, &tbuf, 10);
    PE_Y = strtol(tbuf, &tbuf, 10);
    PE_Z = strtol(tbuf, &tbuf, 10);
    pct.pe_x = PE_X;
    pct.pe_y = PE_Y;
    pct.pe_z = PE_Z;

    /* read coarse grid info */
    get_data("coarse_grid", tbuf, STR, NULL);
    NX_GRID = strtol(tbuf, &tbuf, 10);
    NY_GRID = strtol(tbuf, &tbuf, 10);
    NZ_GRID = strtol(tbuf, &tbuf, 10);


    get_data("potential_grid_refinement",  &FG_RATIO, INT, "2");

    // Set up grid and node data
    set_grids(NX_GRID, NY_GRID, NZ_GRID, PE_X, PE_Y, PE_Z, FG_RATIO);

    get_data("threads_per_node",  &ct.THREADS_PER_NODE, INT, "1");

    FNX_GRID = NX_GRID * FG_RATIO;
    FNY_GRID = NY_GRID * FG_RATIO;
    FNZ_GRID = NZ_GRID * FG_RATIO;

   if(NPES != (PE_X*PE_Y*PE_Z))
        error_handler ("NPES %d  not equal to PE_X*PE_Y*PE_Z!i %d %d %d", NPES, PE_X, PE_Y, PE_Z);

#if GPU_ENABLED
    // Does the system support gpu direct collective operations
    get_data ("gpu_direct_collectives", &ct.gpu_direct_collectives, BOOL, "false");
#endif

    /* Read in the description */
    get_data ("description", &ct.description, STR, "QMD run");

    /* Read in the starting wavefunction file name */
    get_data ("input_wave_function_file", &ct.infile, STR, "Waves/wave.out");

    /* Read in the output wavefunction file name */
    get_data ("output_wave_function_file", &ct.outfile, STR, "Waves/wave.out");


    /* Set up and validate input options */
    char boundary_condition_type_opts[] = "Periodic\n"
                                          "Cluster\n"
                                          "Surface";
    get_data ("boundary_condition_type", NULL, INIT | OPT, boundary_condition_type_opts);

    /* Read in the boundary condition flag */
    get_data ("boundary_condition_type", &ct.boundaryflag, OPT, "Periodic");
    
    
    /*Order of Pulay mixing for charge density*/
    get_data ("charge_pulay_order", &ct.charge_pulay_order, INT, "5");
    
    /*Scale parameter for residuals in Pulay mixing*/
    get_data ("charge_pulay_scale", &ct.charge_pulay_scale, DBL, "0.50");

    /*How often to refresh Pulay history*/
    get_data ("charge_pulay_refresh", &ct.charge_pulay_refresh, INT, "0");
    
    /*Flag to test whether or not the modified metrics should be used in Pulay mixing*/
    get_data ("charge_pulay_special_metrics", &ct.charge_pulay_special_metrics, BOOL, "false");
    
    /*Weight to use for Pulay special metrics*/ 
    get_data ("charge_pulay_special_metrics_weight", &ct.charge_pulay_special_metrics_weight, DBL, "100.0");

    /* Read mixing parameter for linear mixing */
    get_data ("charge_density_mixing", &ct.mix, DBL, "0.5");

    /* Set up and validate input options */
    char exchange_correlation_type_opts[] = "LDA\n"
                                            "GGA BLYP\n"
                                            "GGA XB CP\n"
                                            "GGA XP CP\n"
                                            "GGA PBE\n"
                                            "MGGA TB09";
    get_data ("exchange_correlation_type", NULL, INIT | OPT, exchange_correlation_type_opts);

    /* Exchange correlation potential type flag */
    get_data ("exchange_correlation_type", &ct.xctype, OPT, "LDA");

    /* Number of scf steps */
    get_data ("max_scf_steps", &ct.max_scf_steps, INT, "500");

    /* RMS convergence criterion */
    get_data ("rms_convergence_criterion", &ct.thr_rms, DBL, "1.0E-7");

    /* force convergence criterion */
    get_data ("relax_max_force", &ct.thr_frc, DBL, "2.5E-3");
    
    /* Write pseudopotential plots */
    get_data ("write_pseudopotential_plots", NULL, BOOL, "false");

    /* Write wavefunctions into output file, every main.count of steps */
    get_data ("write_data_period", &ct.checkpoint, INT, "5");


    /* Set up and validate input options */
    char occupations_type_opts[] = "Fixed\n"
                                   "Fermi Dirac\n"
                                   "Gaussian\n"
                                   "Error Function";
    get_data ("occupations_type", NULL, INIT | OPT, occupations_type_opts);

    /* Fermi occupation flag */
    require (get_data ("occupations_type", &ct.occ_flag, OPT, NULL));

    /* Occupation width */
    get_data ("occupation_electron_temperature_eV", &ct.occ_width, DBL, "0.04");
    ct.occ_width *= eV_Ha;

    /* Occupation mixing */
    get_data ("occupation_number_mixing", &ct.occ_mix, DBL, "0.3");


    /*How many steps between writeout of eigenvalues */
    get_data ("write_eigvals_period", &ct.write_eigvals_period, INT, "5");

    /* Set up and validate input options */
    char calculation_mode_opts[] = "Quench Electrons\n"
                                   "Relax Structure\n"
                                   "Constant Volume And Energy\n"
                                   "Constant Temperature And Energy\n"
                                   "Constant Pressure And Energy\n"
                                   "Plot\n"
                                   "Psi Plot\n"
                                   "Band Structure Only\n"
                                   "NEB Relax\n"
                                   "Dimer Relax";
    get_data ("calculation_mode", NULL, INIT | OPT, calculation_mode_opts);

    /* Force flag */
    get_data ("calculation_mode", &ct.forceflag, OPT, "Quench Electrons");



/* do spin polarized calculation? */
    get_data ("spin_polarization", &ct.spin_flag, BOOL, "false");

    get_data ("equal_initial_density", &ct.init_equal_density_flag, BOOL, "false");
    /* Initialized spin up and down charge density equally? */

    /*Flag to write partial density of states*/
    get_data ("write_pdos", &ct.pdos_flag, BOOL, "false");

    /*maximum number of md steps */
    get_data ("max_md_steps", &ct.max_md_steps, INT, "100");


    /* read the electric field vector */
    get_data ("electric_field_vector", tbuf, STR, "0.0 0.0 1.0");
    ct.x_field_0 = strtod (tbuf, &tbuf);
    ct.y_field_0 = strtod (tbuf, &tbuf);
    ct.z_field_0 = strtod (tbuf, &tbuf);
    get_data ("electric_field_magnitude", &ct.e_field, DBL, "0.0");
    /* ct.e_field = sqrt(ct.x_field_0*ct.x_field_0 + ct.y_field_0*ct.y_field_0 + ct.z_field_0*ct.z_field_0); */

    /* Ionic timestep */
    get_data ("ionic_time_step", &ct.iondt, DBL, "50");
    
    /* Max number of sweeps in get_vh*/
    get_data ("hartree_max_sweeps", &ct.hartree_max_sweeps, INT, "100");
    
    /* Min number of sweeps in get_vh*/
    get_data ("hartree_min_sweeps", &ct.hartree_min_sweeps, INT, "5");

    /* Ratio between target RMS for get_vh and RMS total potential
     * This determines target RMS passed to get_vh*/
    get_data ("hartree_rms_ratio", &ct.hartree_rms_ratio, DBL, "50.0");

    /*Whether to use mask mask function for filtering PPs*/
    get_data ("mask_function_filtering", &ct.mask_function, BOOL, "false");


    /* check whether do spin polarized calculation or not*/
    if(ct.spin_flag)
    {
        get_data ("states_count_and_occupation_spin_up", ct.occupation_str_spin_up, STR, NULL);
        get_data ("states_count_and_occupation_spin_down", ct.occupation_str_spin_down, STR, NULL);
    } 
    else
    {
        get_data ("states_count_and_occupation", ct.occupation_str, STR, NULL);

    }
		



    /* Set up and validate input options */
    char crds_units_opts[] = "Bohr\n"
                               "Angstrom";
    get_data ("crds_units", NULL, INIT | OPT, crds_units_opts);

    /*This is not read into any variable */
    get_data ("crds_units", &tmp, OPT, "Bohr");

    /* Set up and validate input options */
    char bravais_lattice_type_opts[] = "None\n"
                                       "Cubic Primitive\n"
                                       "Cubic Face Centered\n"
                                       "Cubic Body Centered\n"
                                       "Hexagonal Primitive\n"
                                       "Hexagonal Rhombohedral (Trigonal)\n"
                                       "Tetragonal Primitive\n"
                                       "Tetragonal Body Centered\n"
                                       "Orthorhombic Primitive\n"
                                       "Orthorhombic Base Centered\n"
                                       "Orthorhombic Body Centered\n"
                                       "Orthorhombic Face Centered\n"
                                       "Monoclinic Primitive\n"
                                       "Monoclinic Base Centered\n"
                                       "Triclinic Primitive";
    get_data ("bravais_lattice_type", NULL, INIT | OPT, bravais_lattice_type_opts);

    /* lattice type */
    require (get_data ("bravais_lattice_type", &ibrav, OPT, NULL));
    set_ibrav_type(ibrav);

    /* Lattice constants */
    /* These should be conditionally read depending on ibrav etc. */
    require (get_data ("a_length", &celldm[0], DBL, NULL));
    require (get_data ("b_length", &celldm[1], DBL, NULL));
    require (get_data ("c_length", &celldm[2], DBL, NULL));
    require (get_data ("alpha", &celldm[3], DBL, NULL));
    require (get_data ("beta", &celldm[4], DBL, NULL));
    require (get_data ("gamma", &celldm[5], DBL, NULL));

    /*Transform to atomic units, which are used internally if input is in angstrom */
    if (verify ("crds_units", "Angstrom"))
    {
        celldm[0] *= A_a0;
        celldm[1] *= A_a0;
        celldm[2] *= A_a0;
    }

    /* Here we read celldm as a,b,c but for most lattice types code uses a, b/a, c/a */
    /* Every lattice type uses a, b/a, c/a except CUBIC_PRIMITIVE, CUBIC_FC and CUBIC_BC */
    if (!verify ("bravais_lattice_type", "Cubic Primitive") &&
        !verify ("bravais_lattice_type", "Cubic Face Centered") &&
        !verify ("bravais_lattice_type", "Cubic Body Centered"))
    {
        celldm[1] /= celldm[0];
        celldm[2] /= celldm[0];
    }
    

    /* initialize the lattice basis vectors */
    flag = 0;

    latgen (celldm, a0, a1, a2, &omega, &flag);

    /* Mehrstellen smoothings pre, post, step, vcycles */
    get_data ("kohn_sham_pre_smoothing", &ct.eig_parm.gl_pre, INT, "2");
    get_data ("kohn_sham_post_smoothing", &ct.eig_parm.gl_pst, INT, "1");
    get_data ("kohn_sham_time_step", &ct.eig_parm.gl_step, DBL, "0.3");
    get_data ("kohn_sham_mucycles", &ct.eig_parm.mucycles, INT, "1");
    get_data ("kohn_sham_fd_order", &ct.kohn_sham_fd_order, INT, "6");
    if((ct.kohn_sham_fd_order != 4) && (ct.kohn_sham_fd_order != 6)) {
        error_handler ("Kohn sham finite difference order must be 4 or 6.");
    }


    /* Poisson smoothings pre, post, step */
    get_data ("poisson_pre_smoothing", &ct.poi_parm.gl_pre, INT, "3");
    get_data ("poisson_post_smoothing", &ct.poi_parm.gl_pst, INT, "3");
    get_data ("poisson_finest_time_step", &ct.poi_parm.gl_step, DBL, "0.6");
    get_data ("poisson_coarse_time_step", &ct.poi_parm.sb_step, DBL, "0.6");
    get_data ("poisson_mucycles", &ct.poi_parm.mucycles, INT, "1");
    get_data ("poisson_coarsest_steps", &ct.poi_parm.coarsest_steps, INT, "80");

    /* Multigrid levels */
    get_data ("kohn_sham_mg_levels", &ct.eig_parm.levels, INT, "1");
    get_data ("poisson_mg_levels", &ct.poi_parm.levels, INT, "2");

    /*Fine grid for non-local pseudopotential */
    get_data ("fine_grid_non_local_pp", &ct.nxfgrid, INT, "4");

    /*Currently, fine grid has to be the same in each direction */
    ct.nzfgrid = ct.nyfgrid = ct.nxfgrid;

    /* Cutoff parameter */
    get_data ("energy_cutoff_parameter", &ct.cparm, DBL, "1.75");
    ct.qcparm = ct.cparm / (rmg_double_t) FG_RATIO;
    ct.betacparm = ct.cparm / (rmg_double_t) ct.nxfgrid;
    ct.cparm /= (rmg_double_t) FG_RATIO; 
    
    /*Blocking factor for scalapack*/
    get_data ("scalapack_block_factor", &ct.scalapack_block_factor, INT, "32");

    /* Whether to write full memory usage report at the end of calculation */
    get_data ("write_memory_report", &ct.write_memory_report, BOOL, "false");

    /* Potential acceleration constant step, default=1.0 */
    get_data ("potential_acceleration_constant_step", &ct.potential_acceleration_constant_step, DBL, "0.0");

    /* Potential acceleration poisson step, default=0.0 */
    get_data ("potential_acceleration_poisson_step", &ct.potential_acceleration_poisson_step, DBL, "0.0");

    /*Count number of species */
    require (get_data ("pseudopotential", &tmp, INIT | LIST, NULL));

    if (!run_count)
    {
	ct.num_species = tmp; 
	my_malloc (ct.sp, ct.num_species, SPECIES);
    }
    else 
    {
	if (tmp != ct.num_species)
	    error_handler("Inconsistency in number of species: %d was specified initially, but %d is given now", ct.num_species, tmp);
    }

    tbuf = tptr;
    is = 0;
    while (get_data ("pseudopotential", tbuf, ITEM | STR, NULL))
    {

        if (sscanf (tbuf, "%s %s", ct.sp[is].pseudo_symbol, ct.sp[is].pseudo_filename) != 2)
        {

            printf ("pseudo_symbol: %s pseudo_filename: %s", ct.sp[is].pseudo_symbol,
                    ct.sp[is].pseudo_filename);

            printf ("\n pseudopotential data: %s is malformed", tbuf);
            error_handler ("Malformed pseudopotential entry in the input file");
        }

        is++;

    }


    /* Set up and validate input options */
    char atomic_coordinate_type_opts[] = "Cell Relative\n"
                                         "Absolute";
    get_data ("atomic_coordinate_type", NULL, INIT | OPT, atomic_coordinate_type_opts);

    /* Absolute or cell relative coordinates */
    get_data ("atomic_coordinate_type", &ct.crd_flag, OPT, "Absolute");



    if ((FNX_GRID / (1 << ct.poi_parm.levels)) < 3)
        error_handler ("NX_GRID: too many hartree MG levels");
    if ((FNY_GRID / (1 << ct.poi_parm.levels)) < 3)
        error_handler ("NY_GRID: too many hartree MG levels");
    if ((FNZ_GRID / (1 << ct.poi_parm.levels)) < 3)
        error_handler ("NZ_GRID: too many hartree MG levels");
    if ((FNX_GRID % (1 << ct.poi_parm.levels)) != 0)
        error_handler ("NX_GRID not evenly divisible by 2^(poi_parm.levels)");
    if ((FNY_GRID % (1 << ct.poi_parm.levels)) != 0)
        error_handler ("NY_GRID not evenly divisible by 2^(poi_parm.levels)");
    if ((FNZ_GRID % (1 << ct.poi_parm.levels)) != 0)
        error_handler ("NZ_GRID not evenly divisible by 2^(poi_parm.levels)");

    /* Clean up malloc'ed memory */
    my_free (tptr);
    

}                               /* end read_control */
