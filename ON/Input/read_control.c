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
 *   see md.h for structure CONTROL
 * PARENTS
 *   md.c
 * CHILDREN
 * 
 * SOURCE */

/* !!!Determine unique list of tags!!! The following command line should be run 
 * and the search terms should be verified as unique if this file is modified
 * 
 * grep get_data read_control.c | cut -d"(" -f2 | cut -d"," -f1 | sort | more
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include "main.h"

static void read_orbitals ();
static void read_kpoints ();


void read_control (void)
{
    int tmp, mpi_nprocs;
    char *tbuf, *tptr;
    FILE *fhand;
    int is,ns;


    /* Open the input file for reading */

    get_data (ct.cfile, NULL, INIT | TAGS, NULL);

    my_malloc (tptr, MAX_PATH, char);
    tbuf = tptr;

    /* Read in the description */
    get_data ("description", &ct.description, STR, "QMD run");

    /* Read in the starting wavefunction file name */
    get_data ("input_wave_function_file", &ct.infile, STR, "wave.out");

    /* Read in the output wavefunction file name */
    get_data ("output_wave_function_file", &ct.outfile, STR, "wave.out");

   char start_mode_opts[] = "Random Start\n"
                            "Restart From File\n"
                            "FIREBALL Start\n"
                            "Gaussian Start";
    get_data ("start_mode", NULL, INIT | OPT, start_mode_opts);


    /* Read in the initial run flag */
    get_data ("start_mode", &ct.runflag, OPT, "Random Start");


    /* Set up and validate input options */
    char boundary_condition_type_opts[] = "Periodic\n"
                                          "Cluster\n"
                                          "Surface";
    get_data ("boundary_condition_type", NULL, INIT | OPT, boundary_condition_type_opts);

    /* Read in the boundary condition flag */
    get_data ("boundary_condition_type", &ct.boundaryflag, OPT, "Periodic");

    /* Read mixing parameter */
    get_data ("charge_density_mixing", &ct.mix, DBL, "0.5");


  /*Order of Pulay mixing for charge density*/
    get_data ("charge_pulay_order", &ct.charge_pulay_order, INT, "1");


    /*How often to refresh Pulay history*/
    get_data ("charge_pulay_refresh", &ct.charge_pulay_refresh, INT, "100");

    /*Flag to test whether or not the modified metrics should be used in * Pulay mixing*/
    get_data ("charge_pulay_special_metrics", &ct.charge_pulay_special_metrics, BOOL, "false");

    /*Weight to use for Pulay special metrics*/
    get_data ("charge_pulay_special_metrics_weight", &ct.charge_pulay_special_metrics_weight, DBL, "100.0");


  /* Set up and validate input options */
    char exchange_correlation_type_opts[] = "LDA\n"
                                            "GGA BLYP\n"
                                            "GGA XB CP\n"
                                            "GGA XP CP\n"
                                            "GGA PBE";
    get_data ("exchange_correlation_type", NULL, INIT | OPT, exchange_correlation_type_opts);


    /* Exchange correlation potential type flag */
    get_data ("exchange_correlation_type", &ct.xctype, OPT, "LDA");

    /* Number of scf steps */
    get_data ("max_scf_steps", &ct.max_scf_steps, INT, "100");
    get_data ("freeze_orbital_step", &ct.freeze_orbital_step, INT, "90");

    /* RMS convergence criterion */
    get_data ("rms_convergence_criterion", &ct.thr_rms, DBL, "1.0E-7");

    /* force convergence criterion */
    get_data ("fast_relax_max_force", &ct.thr_frc, DBL, "2.5E-3");


  /* Write pseudopotential plots */
    get_data ("write_pseudopotential_plots", NULL, BOOL, "false");

    /* Write wavefunctions into output file, every md count of steps */
    get_data ("write_waves_to_file", &ct.checkpoint, BOOL, "true");
    if (ct.checkpoint)
    {
        get_data ("md_steps_til_write_waves", &ct.checkpoint, INT, "10");
        /*If steps_til_write_waves is 0 this is most likely an error
         * 1 was probably intended (i.e.) writing out every step*/
        if (ct.checkpoint == 0)
            ct.checkpoint = 1;
    }


    /* How often calculate energy, print out eigenvalues and occupancies  */
    get_data ("print_energy_and_eigenvalues", &ct.outcount, BOOL, "true");
    if (ct.outcount)
    {
        get_data ("scf_steps_til_energy_and_eigenv_print", &ct.outcount, INT, "1");
        /*If steps_til_energy_and_eigenv_print is 0 this is most likely an error
         * 1 was probably intended (i.e.) writing out every step*/
        if (ct.outcount == 0)
            ct.outcount = 1;
    }

    /* Set up and validate input options */
    char occupations_type_opts[] = "Fixed\n"
                                   "Fermi Dirac\n"
                                   "Gaussian\n"
                                   "Error Function";
    get_data ("occupations_type", NULL, INIT | OPT, occupations_type_opts);


    /* Fermi occupation flag */
    get_data ("occupations_type", &ct.occ_flag, OPT, NULL);

    /* Occupation width */
    get_data ("occupation_electron_temperature_eV", &ct.occ_width, DBL, "0.04");
    ct.occ_width *= eV_Ha;

    /* Occupation mixing */
    get_data ("occupation_number_mixing", &ct.occ_mix, DBL, "0.3");

    /* states occupancy count; overrides background charge */
    get_data("states_count_and_occupation", ct.occupation_str, STR, "");
    get_data("bg_begin", &ct.bg_begin, DBL, "0.0");
    get_data("bg_end", &ct.bg_end, DBL, "10.0");
    get_data("BT", &ct.BT, DBL, "4.0");


   /* Set up and validate input options */
    char calculation_mode_opts[] = "Quench Electrons\n"
                                   "Fast Relax\n"
                                   "Constant Volume And Energy\n"
                                   "Constant Temperature And Energy\n"
                                   "Constant Pressure And Energy\n"
                                   "Plot\n"
                                   "Psi Plot\n"
                                   "Band Structure Only\n"
                                   "NEB Relax";
    get_data ("calculation_mode", NULL, INIT | OPT, calculation_mode_opts);

    /* Force flag */
    get_data ("calculation_mode", &ct.forceflag, OPT, "Quench Electrons");

    /* do spin polarized calculation? */
    get_data ("enable_spin_polarized_calculation", &ct.spin_flag, BOOL, "false");


    /*maximum number of md steps */
    get_data ("max_md_steps", &ct.max_md_steps, INT, "10");

    /* Ionic timestep */
    get_data ("ionic_time_step", &ct.iondt, DBL, "50");

   /* read the electric field vector */
    get_data ("electric_field_vector", tbuf, STR, "0.0 0.0 1.0");
    ct.x_field_0 = strtod (tbuf, &tbuf);
    ct.y_field_0 = strtod (tbuf, &tbuf);
    ct.z_field_0 = strtod (tbuf, &tbuf);
    get_data ("electric_field_magnitude", &ct.e_field, DBL, "0.0");
    ct.constrainforces = 0;

    /*Whether to use mask mask function for filtering PPs*/
    get_data ("mask_function_filtering", &ct.mask_function, BOOL, "false");

    /* -------------------------------- */
    /*                                  */
    /*   order-n code specific inputs   */
    /*                                  */
    /* -------------------------------- */

    /* Set up and validate input options */
    char mixing_opts[] = "Steepest Descent\n"
                         "Pulay\n"
                         "KAIN";
    get_data ("mg_method", NULL, INIT | OPT, mixing_opts);
    /* read mg_eig method */
    get_data("mg_method", &ct.mg_method, OPT, "Pulay");

    /* read mg steps */
    get_data("mg_steps", &ct.mg_steps, INT, "2");

    /* read orbital movable centers option */
    get_data("do_movable_orbital_centers", &ct.movingCenter, BOOL, "false");

    /* read orbital movable centers steps */
    get_data("movable_orbital_centers_steps", &ct.movingSteps, INT, "40");
    if(ct.movingSteps == 0) 
    {
        if(pct.gridpe == 0) 
        {
            printf("\n  *************WARNING: ct.movingSteps reset to 100 \n");
        }
        ct.movingSteps = 100;

    }




    /* Number of states */
    get_data ("number_of_orbitals", &ct.num_states, INT, "0");
    if (ct.num_states > MAX_STATES)
    {
        printf("\n increase MAX_STATES in params.h %d ", ct.num_states);
        error_handler("Too many states specified in input file");
    }


    /* Get k-points and weights */
    read_kpoints ();

    /* Set up and validate input options */
    char length_units_opts[] = "Bohr\n"
                               "Angstrom";
    get_data ("length_units", NULL, INIT | OPT, length_units_opts);

    /*This is not read into any variable */
    get_data ("length_units", &tmp, OPT, "Bohr");


    /* Bravais lattice type */

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
    if (verify ("length_units", "Angstrom"))
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



    /* number of excess electrons in the system (useful for doped systems) */
    get_data("system_charge", &ct.background_charge, DBL, "0");

    /*Background charge is defined to be the opposite of system charge */
    ct.background_charge *= -1.0;



    /* read coarse grid info */
    get_data("coarse_grid", tbuf, STR, NULL);
    NX_GRID = strtol(tbuf, &tbuf, 10);
    NY_GRID = strtol(tbuf, &tbuf, 10);
    NZ_GRID = strtol(tbuf, &tbuf, 10);

    get_data("potential_grid_refinement",  &FG_NX, INT, NULL);

    FG_NY = FG_NX;
    FG_NZ = FG_NX; 

    FNX_GRID = NX_GRID * FG_NX;
    FNY_GRID = NY_GRID * FG_NY;
    FNZ_GRID = NZ_GRID * FG_NZ;

    /*Currently, fine grid has to be the same in each direction */
    get_data("beta_grid_refinement",  &ct.nxfgrid, INT, "4");

    ct.nzfgrid = ct.nyfgrid = ct.nxfgrid;
    BETA_NX = ct.nxfgrid;
    BETA_NY = ct.nyfgrid;
    BETA_NZ = ct.nzfgrid;

    /* read in the processor grid info */

    get_data ("num_processor", &NPES, INT, NULL);
    get_data ("processor_grid", tbuf, STR, NULL);
    pct.pe_x = strtol(tbuf, &tbuf, 10);
    pct.pe_y = strtol(tbuf, &tbuf, 10);
    pct.pe_z = strtol(tbuf, &tbuf, 10);

    PE_X = pct.pe_x;
    PE_Y = pct.pe_y;
    PE_Z = pct.pe_z;
    
    get_data ("kpoints_per_processor", &pct.pe_kpoint, INT, "1");
    get_data ("Hamiltonia_processor_grid", tbuf, STR, "1 1");
    pct.nprow = strtol(tbuf, &tbuf, 10);
    pct.npcol = strtol(tbuf, &tbuf, 10);

    if(NPES != pct.pe_x * pct.pe_y * pct.pe_z ) 
    {
        printf("\n NPES = %d", NPES);
        printf("\n pct.pe_x, y,z = %d %d %d",pct.pe_x, pct.pe_y, pct.pe_z);
        error_handler("bad decomposion of processor grid");
    }

    if(NPES < pct.pe_kpoint * pct.nprow * pct.npcol)
    {
        printf("\n NPES = %d", NPES);
      printf("\n pct.pe_kpoint, pct.nprow, pct.npcol = %d %d %d",pct.pe_kpoint, pct.nprow, pct.npcol);
        error_handler("bad decomposion of processor grid");
    }

    if(pct.nprow > pct.npcol ) 
    {
        printf("\n pct.nprow, pct.npcol = %d %d ", pct.nprow, pct.npcol);
        error_handler("pct.nprow should be smaller than pct.npcol");
    }


    MPI_Comm_size (MPI_COMM_WORLD, &mpi_nprocs);
    if(NPES != mpi_nprocs) 
    {
        printf("\n NPES, mpi_nproc %d %d", NPES, mpi_nprocs);
        error_handler("bad NPES: job and input not match");
    }



    /* Mehrstellen smoothings pre, post, step */
    get_data ("kohn_sham_pre_smoothing", &ct.eig_parm.gl_pre, INT, "2");
    get_data ("kohn_sham_post_smoothing", &ct.eig_parm.gl_pst, INT, "1");
    get_data ("kohn_sham_time_step", &ct.eig_parm.gl_step, DBL, "0.3");

    /* Poisson smoothings pre, post, step */
    get_data ("poisson_pre_smoothing", &ct.poi_parm.gl_pre, INT, "2");
    get_data ("poisson_post_smoothing", &ct.poi_parm.gl_pst, INT, "1");
    get_data ("poisson_time_step", &ct.poi_parm.gl_step, DBL, "0.5");

    /* Multigrid levels */
    get_data ("kohn_sham_mg_levels", &ct.eig_parm.levels, INT, "1");
    get_data ("poisson_mg_levels", &ct.poi_parm.levels, INT, "2");

    PX0_GRID = NX_GRID/pct.pe_x;
    PY0_GRID = NY_GRID/pct.pe_y;
    PZ0_GRID = NZ_GRID/pct.pe_z;
    P0_BASIS = PX0_GRID * PY0_GRID * PZ0_GRID;
    S0_BASIS = (PX0_GRID+2) * (PY0_GRID+2) * (PZ0_GRID+2);


    FPX0_GRID = PX0_GRID * FG_NX;
    FPY0_GRID = PY0_GRID * FG_NY;
    FPZ0_GRID = PZ0_GRID * FG_NZ;
    FP0_BASIS = FPX0_GRID * FPY0_GRID * FPZ0_GRID;

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

    /* spin or not */
    get_data ("do_spin_polarized", &ct.spin, BOOL, "false");

    /* Cutoff parameter */
    get_data ("energy_cutoff_parameter", &ct.cparm, DBL, "1.75");
    ct.qcparm = ct.cparm / (REAL) FG_NX;
    ct.betacparm = ct.cparm / (REAL) ct.nxfgrid;
    ct.cparm /= (REAL) FG_NX;

    /* Output some information for GW calculation. */
    get_data ("output_information_for_GW", &ct.flag_gw, INT, "0");

    /*Count number of species */
    require (get_data ("pseudopotential", &ct.num_species, INIT | LIST, NULL));

    my_malloc (ct.sp, ct.num_species, SPECIES);

    is = 0;
    while (get_data ("pseudopotential", tbuf, ITEM | STR, NULL))
    {

        if (sscanf (tbuf, "%s %s", ct.sp[is].pseudo_symbol,
                    ct.sp[is].pseudo_filename) != 2)
        {

            printf ("pseudo_symbol: %s pseudo_filename: %s",
                    ct.sp[is].pseudo_symbol,
                    ct.sp[is].pseudo_filename);

            printf ("\n pseudopotential data: %s is malformed", tbuf);
            error_handler ("Malformed pseudopotential entry in the input file");
        }

        is++;

    }


    /* read info about atomic orbital filenames for each species for
     * LCAO start */
    char s[32], fn[MAX_PATH];
    require (get_data ("atomic_orbital_files", &ns, INIT | LIST, NULL));
    if(ns != ct.num_species) printf(" \n number of species %d is not equal to number of atomic orbital filesi %d", ct.num_species, ns);

    while (get_data ("atomic_orbital_files", tbuf, ITEM | STR,  NULL)) 
    {
        if (sscanf (tbuf, " %s %s ", s, fn))
        {
            /* search for the species among the known species */
            is = 0;
            while (is < ns && strcmp (s, ct.sp[is].pseudo_symbol))
                is++;

            if (is < ns)        /* we've found it */
                strcpy (ct.file_atomic_orbit[is], fn);
        }
    }


    /* Set up and validate input options */
    char atomic_coordinate_type_opts[] = "Cell Relative\n"
                                         "Absolute";
    get_data ("atomic_coordinate_type", NULL, INIT | OPT, atomic_coordinate_type_opts);

    /* Absolute or cell relative coordinates */
    get_data ("atomic_coordinate_type", &ct.crd_flag, OPT, "Absolute");


    /* read the atomic positions and species */
    read_atoms ();

    read_orbitals ();


    /* Close input file */
    my_free (tptr);



}                               /* end read_control */


static void read_kpoints ()
{
    int ik, nk;
    double w;
    char tbuf [MAX_PATH];

    /* first let's count the number of k-points */
    nk = 0;


    require (get_data ("kpoints", &nk, INIT | LIST, NULL));
    my_malloc (ct.kp, nk, KPOINT);

    /* now we can read the kpoint data */
    ik = 0;
    while (get_data ("kpoints", tbuf, ITEM | STR, NULL))
    {
        sscanf (tbuf, "%lf %lf %lf %lf",
                &ct.kp[ik].kpt[0], &ct.kp[ik].kpt[1], &ct.kp[ik].kpt[2], &ct.kp[ik].kweight);
        ik++;
    }
    if (nk != ik)
    {
        printf ("nk = %d != %d = ik\n", nk, ik);
        error_handler ("Mismatch while reading k-point information\n");
    }

    /*  normalize the kweight to 1.0 if sum(kweights) > 0.0 */
    w = 0.0;
    for (ik = 0; ik < nk; ik++)
        w += ct.kp[ik].kweight;
    if (w > 0.0)
        for (ik = 0; ik < nk; ik++)
            ct.kp[ik].kweight /= w;

    ct.num_kpts = nk;


    /* optional kpoints consistency check */
    if (get_data ("number_of_kpoints", &nk, INT, "0"))
        if (nk != ct.num_kpts)
        {
            printf ("number_of_kpoints = %d != %d = kpoint count\n", nk, ct.num_kpts);
            error_handler ("Mismatch between number_of_kpoints and kpoint count\n");
        }

}




static void read_orbitals ()
{
    int ni, is, ns;
    int i, ist, st1;
    int num_lines;
    double crds[3], radius;
    double bohr;
    int num_tem, movable, frozen;
    char tbuf [MAX_PATH];


    get_data("orbitals", &num_lines, INIT | LIST, NULL);



    bohr = 1.0;
    if (verify ("length_units", "Angstrom"))
    {
        bohr= A_a0;
    }

    /* read and count coordinates for each ion */
    ni = 0;
    ns = 0;
    ist = 0;
    while (get_data ("orbitals", tbuf, ITEM | STR, NULL))
    {
        int args = 0;
        args = sscanf (tbuf, "%d %lf %lf %lf %lf %d %d",
                &num_tem, &crds[0], &crds[1], &crds[2], &radius,
                &movable, &frozen
                );

        if (args < 7)
        {
            printf ("Error reading orbital info %d args %d\n", ni, args);
            error_handler ("Not enough arguments in orbital info list!");
        }

        for(i = 0; i < num_tem; i++)
        {
            states[ni + i].crds[0] = crds[0] * bohr;
            states[ni + i].crds[1] = crds[1] * bohr;
            states[ni + i].crds[2] = crds[2] * bohr;
            states[ni + i].radius = radius;
            states[ni + i].movable = movable;
            states[ni + i].frozen = frozen;
            states[ni + i].n_orbital_same_center = num_tem;
            states[ni + i].gaussian_orbital_index = i;
        }


        ni += num_tem;

    }                         

    ct.num_states = ni;

    /* optional number of atoms consistency check */
    if (get_data ("number_of_orbitals", &ni, INT, "0"))
    {
        if (ni != ct.num_states)
        {
            printf ("number_of_orbitals = %d != %d = orbital info count\n", ni, ct.num_states);
            error_handler ("Mismatch between number_of_orbitals and orbitals count");
        }
    }


}


