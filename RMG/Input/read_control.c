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
 *	void read_control(CONTROL *c)
 * 		read all information from control input file
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
#include "main.h"

void read_control (char *file)
{
    int tmp, is;
    char *tbuf, *tptr;
    float ftmp;
    REAL time1;

    time1 = my_crtc ();
    
    get_data (ct.cfile, NULL, INIT | TAGS, NULL);


    my_malloc (tptr, MAX_PATH, char);
    tbuf = tptr;

    /* Read in the description */
    get_data ("description", &ct.description, STR, "QMD run");

    /* Read in the starting wavefunction file name */
    get_data ("input_wave_function_file", &ct.infile, STR, "Waves/wave.out");

    /* Read in the output wavefunction file name */
    get_data ("output_wave_function_file", &ct.outfile, STR, "Waves/wave.out");

    /* Set up and validate input options */
    char start_mode_opts[] = "Random Start\n"
                             "Restart From File\n"
                             "LCAO Start";
    get_data ("start_mode", NULL, INIT | OPT, start_mode_opts);

    /* Read in the initial run flag */
    get_data ("start_mode", &ct.runflag, OPT, "Random Start");

    /* Set up and validate input options */
    char z_average_output_mode_opts[] = "None\n"
                                        "potential and charge density\n"
                                        "wave functions";
    get_data ("z_average_output_mode", NULL, INIT | OPT, z_average_output_mode_opts);

    /* Read in zaverage output toggle */
    get_data ("z_average_output_mode", &ct.zaverage, OPT, "None");


    /* Set up and validate input options */
    char boundary_condition_type_opts[] = "Periodic\n"
                                          "Cluster\n"
                                          "Surface";
    get_data ("boundary_condition_type", NULL, INIT | OPT, boundary_condition_type_opts);

    /* Read in the boundary condition flag */
    get_data ("boundary_condition_type", &ct.boundaryflag, OPT, "Periodic");
    
    
    /* Set up and validate input options */
    char charge_mixing_type_opts[] = "Linear\n"
                                     "Pulay";
    get_data ("charge_mixing_type", NULL, INIT | OPT, charge_mixing_type_opts);
    
    /* Read type of charge density mixing */
    get_data ("charge_mixing_type", NULL, OPT, "Pulay");

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

    /* Read in the prjmixing parameter */
    get_data ("projector_mixing", &ct.prjmix, DBL, "0.5");

    /* Set up and validate input options */
    char exchange_correlation_type_opts[] = "LDA\n"
                                            "GGA BLYP\n"
                                            "GGA XB CP\n"
                                            "GGA XP CP\n"
                                            "GGA PBE";
    get_data ("exchange_correlation_type", NULL, INIT | OPT, exchange_correlation_type_opts);

    /* Exchange correlation potential type flag */
    get_data ("exchange_correlation_type", &ct.xctype, OPT, "LDA");

    /* Emin when do get_dos */
    get_data ("Emin", &ct.Emin, DBL, "-6.0");

    /* Emax when do get_dos */
    get_data ("Emax", &ct.Emax, DBL, "0.0");

    /* Number Energy points when do get_dos */
    get_data ("E_POINTS", &ct.E_POINTS, INT, "201");

    /* Number of scf steps */
    get_data ("max_scf_steps", &ct.max_scf_steps, INT, "500");

    /* RMS convergence criterion */
    get_data ("rms_convergence_criterion", &ct.thr_rms, DBL, "1.0E-7");

    /* force convergence criterion */
    get_data ("relax_max_force", &ct.thr_frc, DBL, "2.5E-3");
    
    char relax_mass_opts[] = "Atomic\n"
                                 "Equal";
    get_data ("relax_mass", NULL, INIT | OPT, relax_mass_opts);

    /* Mass to use for structural relaxation, either their atomic masses, or use the mass of carbon for all atoms*/
    get_data ("relax_mass", &ct.relax_mass, OPT, "Atomic");

    /* Write pseudopotential plots */
    get_data ("write_pseudopotential_plots", NULL, BOOL, "false");

    /* Write wavefunctions into output file, every main.count of steps */
    get_data ("write_data_period", &ct.checkpoint, INT, "5");


    /* sorting flag */
    get_data ("sort_wavefunctions", &ct.sortflag, BOOL, "true");
    

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

    /* Diagonalization info: initial diagonalization flag */
    get_data ("initial_diagnolization", &ct.initdiag, BOOL, "false");

    /* Diagonalization info: diagonalization period */
    get_data ("period_of_diagonalization", &ct.diag, INT, "10");

    /* stop diagonalizing after end_diag steps */
    get_data ("end_diagonalization_step", &ct.end_diag, INT, "100");

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

    /* Set up and validate input options */
    char relax_method_opts[] = "Fast Relax\n"
				"FIRE\n"
				"Quick Min\n"
				"MD Min\n"
                               "LBFGS";
    get_data ("relax_method", NULL, INIT | OPT, relax_method_opts);

    get_data ("relax_method", &ct.relax_method, OPT, "Fast Relax");
    
    /*Whether or not to use dynamic timestep in relax mode*/
    get_data ("relax_dynamic_timestep", NULL, BOOL, "false");



    if ( verify( "calculation_mode", "NEB Relax" ) )
    {
        /* set constraint type for switch in Common/constrain.c */
        ct.constrainforces = 5;
    }

    /* NEB spring constant */
    get_data ("neb_spring_constant", &ct.neb_spring_constant, DBL, "0.5");

    /* do spin polarized calculation? */
    get_data ("spin_polarization", &ct.spin_flag, BOOL, "false");

    get_data ("equal_initial_density", &ct.init_equal_density_flag, BOOL, "false");
    /* Initialized spin up and down charge density equally? */

    /*Flag to write partial density of states*/
    get_data ("write_pdos", &ct.pdos_flag, BOOL, "false");

    /*maximum number of md steps */
    get_data ("max_md_steps", &ct.max_md_steps, INT, "100");

    /* Retrieve number of rmg "restart like" (NEB/exchange/ARTs) steps to perform */
    get_data ("max_rmg_steps", &ct.max_rmg_steps, INT, "1" );

    /* read the electric field vector */
    get_data ("electric_field_vector", tbuf, STR, "0.0 0.0 1.0");
    ct.x_field_0 = strtod (tbuf, &tbuf);
    ct.y_field_0 = strtod (tbuf, &tbuf);
    ct.z_field_0 = strtod (tbuf, &tbuf);
    get_data ("electric_field_magnitude", &ct.e_field, DBL, "0.0");
    /* ct.e_field = sqrt(ct.x_field_0*ct.x_field_0 + ct.y_field_0*ct.y_field_0 + ct.z_field_0*ct.z_field_0); */
    tbuf = tptr;

    /* Set up and validate input options */
    char md_temperature_control_opts[] = "Nose Hoover Chains\n"
                                          "Anderson Rescaling";
    get_data ("md_temperature_control", NULL, INIT | OPT, md_temperature_control_opts);

    /* Set up and validate input options */
    /* Temperature Control Info */
    get_data ("md_temperature_control", &ct.tcontrol, OPT, "Nose Hoover Chains");

    /* Target MD Temperature */
    get_data ("md_temperature", &ct.nose.temp, DBL, "300");

    /* Nose Thermostats */
    get_data ("md_number_of_nose_thermostats", &ct.nose.m, INT, "5");

    /* Nose oscillation frequency */
    get_data ("md_nose_oscillation_frequency_THz", &ct.nose.fNose, DBL, "15.59");


    /* Nose randomize velocity flag */
    /* Default value depend on whether we do restart or not
     * If restarting, default is false, velocities from previous run should be good enough
     * if not restarting, defualt value will be true*/
    if (verify ("start_mode", "Restart From File"))
        get_data ("md_randomize_velocity", &ct.nose.randomvel, BOOL, "false");
    else
        get_data ("md_randomize_velocity", &ct.nose.randomvel, BOOL, "true");


    /* Set up and validate input options */
    char md_integration_order_opts[] = "2nd Velocity Verlet\n"
                                        "3rd Beeman-Velocity Verlet\n"
                                        "5th Beeman-Velocity Verlet";
    get_data ("md_integration_order", NULL, INIT | OPT, md_integration_order_opts);
    /* MD Integration flag */
    get_data ("md_integration_order", &ct.mdorder, OPT, "2nd Velocity Verlet");

    /* Ionic timestep */
    get_data ("ionic_time_step", &ct.iondt, DBL, "50");
    
    /* Maximum ionic timestep */
    get_data ("max_ionic_time_step", &ct.iondt_max, DBL, "150");

    if (ct.iondt_max < ct.iondt)
        error_handler("max_ionic_time_step (%f) has to be >= than ionic_time_step (%f)", ct.iondt_max, ct.iondt); 


    /*** Factor by which iondt is increased */
    get_data ("ionic_time_step_increase", &ct.iondt_inc, DBL, "1.1");
    
    /*Estimate default decrease factor, so that it takes about 3 decrease steps 
     * to go from max to starting time step*/
    {
        REAL ttt;
        char s_ttt[12] = { '\0' };

        ttt = pow (ct.iondt_max / ct.iondt, 0.3);
        ttt = 1.0/ ttt;

        if (ttt > 0.75) ttt = 0.75;
        if (ttt < 0.4) ttt = 0.4;

        snprintf(s_ttt, sizeof(s_ttt) - 1, "%f", ttt);
        /*** Factor by which iondt is decreased */
        get_data ("ionic_time_step_decrease", &ct.iondt_dec, DBL, s_ttt);
    }
    
    /*Number of steps after which iondt is increased (if F * v is positive for all those steps) */
    get_data ("dynamic_time_delay", &ct.relax_steps_delay, INT, "5");
    
    /*Current value of step counter, when this value reaches ct.relax_steps_delay timestep can be increased
     * This normally starts from zero and is increased if F *v > 0 or set to 0 otherwise*/
    get_data ("dynamic_time_counter", &ct.relax_steps_counter, INT, "0");



    /* Milliken population flag */
    get_data ("milliken_population", &ct.domilliken, BOOL, "false");

    /* Max number of sweeps in get_vh*/
    get_data ("hartree_max_sweeps", &ct.hartree_max_sweeps, INT, "100");
    
    /* Min number of sweeps in get_vh*/
    get_data ("hartree_min_sweeps", &ct.hartree_min_sweeps, INT, "5");

    /* Ratio between target RMS for get_vh and RMS total potential
     * This determines target RMS passed to get_vh*/
    get_data ("hartree_rms_ratio", &ct.hartree_rms_ratio, DBL, "50.0");

    /*Whether to use mask mask function for filtering PPs*/
    get_data ("mask_function_filtering", &ct.mask_function, BOOL, "false");


    /* Number of states */
    if(ct.spin_flag)
    {
    	get_data ("states_per_kpoint_up", &ct.num_states_up, INT, "0");
    	get_data ("states_per_kpoint_down", &ct.num_states_down, INT, "0");
    } 
    else
    {
    	get_data ("states_per_kpoint", &ct.num_states, INT, "0");

    }
    get_data ("states_per_kpoint", &ct.num_states, INT, "0");

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
		


    if (strcmp (ct.occupation_str, "") == 0)
    {

        /* number of unoccupied states above the Fermi level */
        get_data ("unoccupied_states_per_kpoint", &ct.num_unocc_states, INT, "10");

        /* number of excess electrons in the system (useful for doped systems) */
        get_data ("system_charge", &ct.background_charge, DBL, "0");

        /*Background charge is defined to be the opposite of system charge */
        ct.background_charge *= -1.0;

    }




    /* Set up and validate input options */
    char interpolation_type_opts[] = "Cubic Polynomial\n"
                                     "B-spline\n"
                                     "prolong";
    get_data ("interpolation_type", NULL, INIT | OPT, interpolation_type_opts);

    /*Interpolation type */
    get_data ("interpolation_type", &ct.interp_flag, OPT, "Cubic Polynomial");

    /*B-spline interpolation order */
    get_data ("b_spline_order", &ct.interp_order, INT, "5");

    /*B-spline order must be 8 or less */
    Dprintf ("Interpolation order was set to %d", ct.interp_order);
    if (ct.interp_order > 8)
	error_handler ("Interpolation order too high");
    if (ct.interp_order <= 0)
	error_handler ("Invalid interpolation order");

    get_data ("b_spline_trade_order", &ct.interp_trade, INT, "3");


    /* Get k-points and weights */
    require (get_data ("kpoints", &ct.num_kpts, INIT | LIST, NULL));
    my_malloc (ct.kp, ct.num_kpts, KPOINT);

    /* now we can read the kpoint data */
    tmp = 0;
    while (get_data ("kpoints", tbuf, ITEM | STR, NULL))
    {
        sscanf (tbuf, "%lf %lf %lf %lf",
                &ct.kp[tmp].kpt[0], &ct.kp[tmp].kpt[1], &ct.kp[tmp].kpt[2], &ct.kp[tmp].kweight);
        tmp++;
    }
    if (ct.num_kpts != tmp)
    {
        error_handler ("Count mismatch while reading k-point information\n");
    }

    /*  normalize the kweight to 1.0 if sum(kweights) > 0.0 */
    ftmp = 0.0;
    for (tmp = 0; tmp < ct.num_kpts; tmp++)
        ftmp += ct.kp[tmp].kweight;
    if (ftmp > 0.0)
        for (tmp = 0; tmp < ct.num_kpts; tmp++)
            ct.kp[tmp].kweight /= ftmp;


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
    require (get_data ("bravais_lattice_type", &ct.ibrav, OPT, NULL));

    /* Lattice constants */
    /* These should be conditionally read depending on ibrav etc. */
    require (get_data ("a_length", &ct.celldm[0], DBL, NULL));
    require (get_data ("b_length", &ct.celldm[1], DBL, NULL));
    require (get_data ("c_length", &ct.celldm[2], DBL, NULL));
    require (get_data ("alpha", &ct.celldm[3], DBL, NULL));
    require (get_data ("beta", &ct.celldm[4], DBL, NULL));
    require (get_data ("gamma", &ct.celldm[5], DBL, NULL));

    /*Transform to atomic units, which are used internally if input is in angstrom */
    if (verify ("crds_units", "Angstrom"))
    {
        ct.celldm[0] *= A_a0;
        ct.celldm[1] *= A_a0;
        ct.celldm[2] *= A_a0;
    }

    /* Here we read celldm as a,b,c but for most lattice types code uses a, b/a, c/a */
    /* Every lattice type uses a, b/a, c/a except CUBIC_PRIMITIVE, CUBIC_FC and CUBIC_BC */
    if (!verify ("bravais_lattice_type", "Cubic Primitive") &&
        !verify ("bravais_lattice_type", "Cubic Face Centered") &&
        !verify ("bravais_lattice_type", "Cubic Body Centered"))
    {
        ct.celldm[1] /= ct.celldm[0];
        ct.celldm[2] /= ct.celldm[0];
    }



    /* Mehrstellen smoothings pre, post, step, vcycles */
    get_data ("kohn_sham_pre_smoothing", &ct.eig_parm.gl_pre, INT, "2");
    get_data ("kohn_sham_post_smoothing", &ct.eig_parm.gl_pst, INT, "1");
    get_data ("kohn_sham_time_step", &ct.eig_parm.gl_step, DBL, "0.3");
    get_data ("kohn_sham_mucycles", &ct.eig_parm.mucycles, INT, "1");

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

    if ((PX0_GRID / (1 << ct.eig_parm.levels)) < 3)
        error_handler ("PX0_GRID: too many eigenvalue MG levels");
    if ((PY0_GRID / (1 << ct.eig_parm.levels)) < 3)
        error_handler ("PY0_GRID: too many eigenvalue MG levels");
    if ((PZ0_GRID / (1 << ct.eig_parm.levels)) < 3)
        error_handler ("PZ0_GRID: too many eigenvalue MG levels");
    if ((PX0_GRID % (1 << ct.eig_parm.levels)) != 0)
        error_handler ("PX0_GRID not evenly divisible by 2^(eig_parm.levels)");
    if ((PY0_GRID % (1 << ct.eig_parm.levels)) != 0)
        error_handler ("PY0_GRID not evenly divisible by 2^(eig_parm.levels)");
    if ((PZ0_GRID % (1 << ct.eig_parm.levels)) != 0)
        error_handler ("PZ0_GRID not evenly divisible by 2^(eig_parm.levels)");

    if ((FPX0_GRID / (1 << ct.poi_parm.levels)) < 3)
        error_handler ("PX0_GRID: too many hartree MG levels");
    if ((FPY0_GRID / (1 << ct.poi_parm.levels)) < 3)
        error_handler ("PY0_GRID: too many hartree MG levels");
    if ((FPZ0_GRID / (1 << ct.poi_parm.levels)) < 3)
        error_handler ("PZ0_GRID: too many hartree MG levels");
    if ((FPX0_GRID % (1 << ct.poi_parm.levels)) != 0)
        error_handler ("PX0_GRID not evenly divisible by 2^(poi_parm.levels)");
    if ((FPY0_GRID % (1 << ct.poi_parm.levels)) != 0)
        error_handler ("PY0_GRID not evenly divisible by 2^(poi_parm.levels)");
    if ((FPZ0_GRID % (1 << ct.poi_parm.levels)) != 0)
        error_handler ("PZ0_GRID not evenly divisible by 2^(poi_parm.levels)");


    /*Fine grid for non-local pseudopotential */
    get_data ("fine_grid_non_local_pp", &ct.nxfgrid, INT, "4");

    /*Currently, fine grid has to be the same in each direction */
    ct.nzfgrid = ct.nyfgrid = ct.nxfgrid;

    /* Cutoff parameter */
    get_data ("energy_cutoff_parameter", &ct.cparm, DBL, "1.75");
    ct.qcparm = ct.cparm / (REAL) FG_NX;
    ct.betacparm = ct.cparm / (REAL) ct.nxfgrid;
    ct.cparm /= (REAL) FG_NX; 

    /* Whether to write full memory usage report at the end of calculation */
    get_data ("write_memory_report", &ct.write_memory_report, BOOL, "false");

    /*Count number of species */
    require (get_data ("pseudopotential", &ct.num_species, INIT | LIST, NULL));

    my_malloc (ct.sp, ct.num_species, SPECIES);

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

    /* read the atomic positions and species */
    tbuf = tptr;


    /*Test for presence of pdb_atoms or atoms
     * if none, error_handler*/
    if (verify ("pdb_atoms", NULL) || verify ("atoms", NULL))
    {
        /* Test for presence of pdb_atom tag */
        if (get_data ("pdb_atoms", &tmp, INIT | LIST, NULL))
        {

            read_pdb ();

        }
        if (verify ("atoms", NULL))
            read_atoms ();

    }
    else
    {
        error_handler ("No ion coordinates available");
    }


    /*Read or initialize ionic velocities*/
    if (verify("ionic_velocities", NULL))
    {
	get_data ("ionic_velocities", NULL, INIT | LIST, NULL);
	
	tmp = 0;
	while (get_data ("ionic_velocities", tbuf, ITEM | STR, NULL))
	{
	    sscanf (tbuf, "%lf %lf %lf",
		    &ct.ions[tmp].velocity[0], &ct.ions[tmp].velocity[1], &ct.ions[tmp].velocity[2]);
	    tmp++;

	    if (tmp >= ct.num_ions) error_handler("Velocities specified for too many atoms");
	}

	if (tmp !=  ct.num_ions) 
	    error_handler("Mismatch between number of velocities (%d) and number of ions (%d)", tmp, ct.num_ions); 
    }

    else
    {
    
	for (tmp = 0; tmp < ct.num_ions; tmp++)
	{
	    ct.ions[tmp].velocity[0] = 0.0;
	    ct.ions[tmp].velocity[1] = 0.0;
	    ct.ions[tmp].velocity[2] = 0.0;
	}

    }
    
    
    /*Read or initialize ionic forces*/
    if (verify("ionic_forces", NULL))
    {
	get_data ("ionic_forces", NULL, INIT | LIST, NULL);
	tmp = 0;
	while (get_data ("ionic_forces", tbuf, ITEM | STR, NULL))
	{
	    sscanf (tbuf, "%lf %lf %lf",
		    &ct.ions[tmp].force[0][0], &ct.ions[tmp].force[0][1], &ct.ions[tmp].force[0][2]);
	    tmp++;

	    if (tmp >= ct.num_ions) error_handler("Forces specified for too many atoms");
	}

	if (tmp !=  ct.num_ions) 
	    error_handler("Mismatch between number of forces (%d) and number of ions (%d)", tmp, ct.num_ions); 
    }

    else
    {
    
	for (tmp = 0; tmp < ct.num_ions; tmp++)
	{
	    ct.ions[tmp].force[0][0] = 0.0;
	    ct.ions[tmp].force[0][1] = 0.0;
	    ct.ions[tmp].force[0][2] = 0.0;
	}

    }

    /*Read or intialize ionic force pointer*/
    if (verify("ionic_force_pointer", NULL))
    {
	get_data ("ionic_force_pointer", tbuf, ITEM | STR, NULL);
	sscanf (tbuf, "%d %d %d", &ct.fpt[0], &ct.fpt[1], &ct.fpt[2], &ct.fpt[3]);
    }
    else
    {
        ct.fpt[0] = 0;
        ct.fpt[1] = 1;
        ct.fpt[2] = 2;
        ct.fpt[3] = 3;
    }
    
    get_data ("ionic_force_pointer", tbuf, STR, "0 1 2 3");
    ct.fpt[0] = strtol (tbuf, &tbuf, 10);
    ct.fpt[1] = strtol (tbuf, &tbuf, 10);
    ct.fpt[2] = strtol (tbuf, &tbuf, 10);
    ct.fpt[3] = strtol (tbuf, &tbuf, 10);
    tbuf = tptr;

    get_data ("nose_positions", tbuf, STR, "0 0 0 0 0 0 0 0 0 0");
    for (tmp = 0; tmp < 10; tmp++)
	ct.nose.xx[tmp] = strtol (tbuf, &tbuf, 10);
    tbuf = tptr;
    
    get_data ("nose_velocities", tbuf, STR, "0 0 0 0 0 0 0 0 0 0");
    for (tmp = 0; tmp < 10; tmp++)
	ct.nose.xv[tmp] = strtol (tbuf, &tbuf, 10);
    tbuf = tptr;
    
    get_data ("nose_masses", tbuf, STR, "0 0 0 0 0 0 0 0 0 0");
    for (tmp = 0; tmp < 10; tmp++)
	ct.nose.xq[tmp] = strtol (tbuf, &tbuf, 10);
    tbuf = tptr;
    
    
	
    tmp = 0;
    get_data ("nose_forces", tbuf, INIT | LIST, NULL);
    while (get_data ("nose_forces", tbuf, ITEM | STR, NULL))
    {
	sscanf (tbuf, " %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&ct.nose.xf[tmp][0], &ct.nose.xf[tmp][1], &ct.nose.xf[tmp][2], &ct.nose.xf[tmp][3], &ct.nose.xf[tmp][4], 
		&ct.nose.xf[tmp][5], &ct.nose.xf[tmp][6], &ct.nose.xf[tmp][7], &ct.nose.xf[tmp][8], &ct.nose.xf[tmp][9]);
	tmp++;

	if (tmp >= 4) error_handler("Too many nose_forces, only 4 expected");
    }

    if ((tmp > 0) && (tmp !=  4)) 
	error_handler("Wrong number of Nose force (%d), 4 are expected", tmp); 





    /* get? multi-image topology if present in input */
    if (get_data( "image_communicator_topology", &tmp, INIT|LIST, NULL))



    /* get? multi-image topology if present in input */
    if (get_data( "image_communicator_topology", &tmp, INIT|LIST, NULL))



    /* get? multi-image topology if present in input */
    if (get_data( "image_communicator_topology", &tmp, INIT|LIST, NULL))
    {
        init_img_topo( tmp );
    }
    else
    {
        init_img_topo( 0 );
    }

    /* get? the vectors for constrained forces */
    if (get_data ("atom_constraints", &tmp, INIT | LIST, NULL))
    {
        if (ct.constrainforces == 0)
            ct.constrainforces = 1;

        if (tmp != ct.num_ions)
            error_handler ("Mismatch in number of atom constraints, check input.");

        /* Parse and normalize all constraint inputs */
        tbuf = tptr;
        tmp = 0;
        while (get_data ("atom_constraints", tbuf, ITEM | STR, NULL))
        {
            if (sscanf (tbuf, "%lf %lf %lf %lf",
                        &ct.ions[tmp].constraint.setA_coord[0],
                        &ct.ions[tmp].constraint.setA_coord[1],
                        &ct.ions[tmp].constraint.setA_coord[2],
                        &ct.ions[tmp].constraint.setA_weight) != 4)
            {
                printf ("\n Atom constraint data: %s is malformed", tbuf);
                error_handler ("Malformed constraint entry in the input file");
            }

            ftmp = 0.0;
            ftmp += ct.ions[tmp].constraint.setA_coord[0] * ct.ions[tmp].constraint.setA_coord[0];
            ftmp += ct.ions[tmp].constraint.setA_coord[1] * ct.ions[tmp].constraint.setA_coord[1];
            ftmp += ct.ions[tmp].constraint.setA_coord[2] * ct.ions[tmp].constraint.setA_coord[2];
            if (ftmp != ZERO)
            {
                ftmp = sqrt (ftmp);
                ct.ions[tmp].constraint.setA_coord[0] /= ftmp;
                ct.ions[tmp].constraint.setA_coord[1] /= ftmp;
                ct.ions[tmp].constraint.setA_coord[2] /= ftmp;
            }
            tmp++;
        }

        /* Make sure no errors in conversion */
        if (tmp != ct.num_ions)
        {
            printf ("Number of constraints is = %d != %d = ct.num_ions", is, ct.num_ions);
            error_handler ("Mismatch in number of ions and number of read constraints");
        }
    }

    /* Clean up malloc'ed memory */
    my_free (tptr);
    
    rmg_timings (READ_CONTROL_TIME, (my_crtc () - time1));

}                               /* end read_control */
