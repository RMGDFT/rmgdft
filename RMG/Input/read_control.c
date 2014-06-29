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
#include "grid.h"
#include "common_prototypes.h"
#include "main.h"


void read_control (char *file)
{
    int tmp, is;
    char *tbuf, *tptr;
    float ftmp;
    static int run_count = - 1;

    run_count ++;
    
    get_data (file, NULL, INIT | TAGS, NULL);

    read_common();
    my_malloc (tptr, MAX_PATH, char);
    tbuf = tptr;


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
    char charge_mixing_type_opts[] = "Linear\n"
                                     "Pulay";
    get_data ("charge_mixing_type", NULL, INIT | OPT, charge_mixing_type_opts);
    
    /* Read type of charge density mixing */
    get_data ("charge_mixing_type", NULL, OPT, "Pulay");

    /* Read in the prjmixing parameter */
    get_data ("projector_mixing", &ct.prjmix, DBL, "0.5");

    /* Emin when do get_dos */
    get_data ("Emin", &ct.Emin, DBL, "-6.0");

    /* Emax when do get_dos */
    get_data ("Emax", &ct.Emax, DBL, "0.0");

    /* Number Energy points when do get_dos */
    get_data ("E_POINTS", &ct.E_POINTS, INT, "201");


    
    char relax_mass_opts[] = "Atomic\n"
                                 "Equal";
    get_data ("relax_mass", NULL, INIT | OPT, relax_mass_opts);

    /* Mass to use for structural relaxation, either their atomic masses, or use the mass of carbon for all atoms*/
    get_data ("relax_mass", &ct.relax_mass, OPT, "Atomic");


    /* sorting flag */
    get_data ("sort_wavefunctions", &ct.sortflag, BOOL, "true");
    
    /* Diaonalization opts */
    char diagonalization_driver_opts[] = "lapack\n"
                                         "scalapack\n"
                                         "magma\n"
                                         "magmafs\n"
                                         "lapackfs\n";
//                                       "elpa\n"
    get_data ("subdiag_driver", NULL, INIT | OPT, diagonalization_driver_opts);
    get_data ("subdiag_driver", &ct.subdiag_driver, OPT, "scalapack");
#if !MAGMA_LIBS
    if(verify( "subdiag_driver", "magma" ))
          error_handler("    This version of RMG was not built with Magma.\n");
    if(verify( "subdiag_driver", "magmafs" ))
          error_handler("    This version of RMG was not built with Magma.\n");
#endif

    /* Diagonalization info: initial diagonalization flag */
    get_data ("initial_diagonalization", &ct.initdiag, BOOL, "false");

    /* Diagonalization info: diagonalization period */
    get_data ("period_of_diagonalization", &ct.diag, INT, "10");

    /* stop diagonalizing after end_diag steps */
    get_data ("end_diagonalization_step", &ct.end_diag, INT, "100");

    /* Read in the folded spectrum width */
    get_data ("folded_spectrum_width", &ct.folded_spectrum_width, DBL, "0.3");


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

    /*Flag to write partial density of states*/
    get_data ("write_pdos", &ct.pdos_flag, BOOL, "false");

    /* Retrieve number of rmg "restart like" (NEB/exchange/ARTs) steps to perform */
    get_data ("max_rmg_steps", &ct.max_rmg_steps, INT, "1" );

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
        rmg_double_t ttt;
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

    get_data("scf_steps_offset", &tmp, INT, "0"); 
    ct.scf_steps = tmp;
    
    get_data("total_scf_steps_offset", &tmp, INT, "0"); 
    ct.total_scf_steps = tmp;
    
    get_data("md_steps_offset", &tmp, INT, "0"); 
    ct.md_steps = tmp;




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
    require (get_data ("kpoints", &tmp, INIT | LIST, NULL));
    if (!run_count)
    {
	ct.num_kpts = tmp;
	my_malloc (ct.kp, ct.num_kpts, KPOINT);
    }
    else
    {
	if (tmp != ct.num_kpts)
	    error_handler("Inconsistency in number of kpoints: %d was specified initially, but %d is given now", ct.num_kpts, tmp);
    }

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
    




    /*Fine grid for non-local pseudopotential */
    get_data ("fine_grid_non_local_pp", &ct.nxfgrid, INT, "4");

    /*Currently, fine grid has to be the same in each direction */
    ct.nzfgrid = ct.nyfgrid = ct.nxfgrid;

    /* Cutoff parameter */
    get_data ("energy_cutoff_parameter", &ct.cparm, DBL, "1.75");
    ct.qcparm = ct.cparm / (rmg_double_t) get_FG_RATIO();
    ct.betacparm = ct.cparm / (rmg_double_t) ct.nxfgrid;
    ct.cparm /= (rmg_double_t) get_FG_RATIO(); 
    
    // Norm conserving pseudo potential flag
    get_data ("norm_conserving_pp", &ct.norm_conserving_pp, BOOL, "false");
    
    /* scalapack allreduce or point to point for dist matrices */
    get_data ("scalapack_global_sums", &ct.scalapack_global_sums, BOOL, "true");





    /* read the atomic positions and species */
    tbuf = tptr;


    /*Test for presence of pdb_atoms or atoms
     * if none, error_handler*/
    if (verify ("pdb_atoms", NULL) || verify ("atoms", NULL))
    {
	int num_ions = 0;
	
	get_data ("atoms", &num_ions, INIT | LIST, NULL);
	tmp = 0;
	get_data ("pdb_atoms", &tmp, INIT | LIST, NULL);

	/*Figure out number of atoms, number of pdb_atoms should equal to number of "atoms" unless one of them is zero*/
	if (num_ions)
	{
	    if (tmp) 
		if (tmp != num_ions)
		    error_handler("pdb_atoms and atoms specify different numbers of atoms: %d and %d, respectively", tmp, num_ions);
	}
	else
	{
	    if (!tmp) error_handler("Both pdb_atoms and atoms specify zero ions");
	    num_ions = tmp;
	}

	
	if (!run_count)
	{
	    ct.num_ions = num_ions;
	    my_calloc (ct.ions, ct.num_ions, ION);
	}
	else
	    if (num_ions != ct.num_ions)
		error_handler("Inconsistency in number of ions: %d was specified initially, but %d is given now", ct.num_ions, num_ions);

	/* Test for presence of pdb_atom tag */
        if (get_data ("pdb_atoms", &tmp, INIT | LIST, NULL))
	    read_pdb ();

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
	    if (tmp >= ct.num_ions) error_handler("# of velocities(%d) >= # of atoms(%d)", tmp, ct.num_ions);
	    
	    sscanf (tbuf, "%lf %lf %lf",
		    &ct.ions[tmp].velocity[0], &ct.ions[tmp].velocity[1], &ct.ions[tmp].velocity[2]);
	    tmp++;

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
	    if (tmp >= ct.num_ions) error_handler("Forces specified for too many atoms");
	    
	    sscanf (tbuf, "%lf %lf %lf",
		    &ct.ions[tmp].force[0][0], &ct.ions[tmp].force[0][1], &ct.ions[tmp].force[0][2]);
	    tmp++;

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
	sscanf (tbuf, "%d %d %d %d", &ct.fpt[0], &ct.fpt[1], &ct.fpt[2], &ct.fpt[3]);
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
	if (tmp >= 4) error_handler("Too many nose_forces, only 4 expected");
	
	sscanf (tbuf, " %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&ct.nose.xf[tmp][0], &ct.nose.xf[tmp][1], &ct.nose.xf[tmp][2], &ct.nose.xf[tmp][3], &ct.nose.xf[tmp][4], 
		&ct.nose.xf[tmp][5], &ct.nose.xf[tmp][6], &ct.nose.xf[tmp][7], &ct.nose.xf[tmp][8], &ct.nose.xf[tmp][9]);
	tmp++;

    }

    if ((tmp > 0) && (tmp !=  4)) 
	error_handler("Wrong number of Nose force (%d), 4 are expected", tmp); 


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


    // We check multigrid levels last since we have to be sure that the grid dims have already 
    // been read in.
    if ((get_NX_GRID() / (1 << ct.eig_parm.levels)) < 3)
        error_handler ("NX_GRID: too many eigenvalue MG levels");
    if ((get_NY_GRID() / (1 << ct.eig_parm.levels)) < 3)
        error_handler ("NY_GRID: too many eigenvalue MG levels");
    if ((get_NZ_GRID() / (1 << ct.eig_parm.levels)) < 3)
        error_handler ("NZ_GRID: too many eigenvalue MG levels");
    if ((get_NX_GRID() % (1 << ct.eig_parm.levels)) != 0)
        error_handler ("NX_GRID not evenly divisible by 2^(eig_parm.levels)");
    if ((get_NY_GRID() % (1 << ct.eig_parm.levels)) != 0)
        error_handler ("NY_GRID not evenly divisible by 2^(eig_parm.levels)");
    if ((get_NZ_GRID() % (1 << ct.eig_parm.levels)) != 0)
        error_handler ("NZ_GRID not evenly divisible by 2^(eig_parm.levels)");


    /* Clean up malloc'ed memory */
    my_free (tptr);
    

}                               /* end read_control */
