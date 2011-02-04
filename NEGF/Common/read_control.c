/************************** SVN Revision Information **************************
 **    $Id: read_control.c 1242 2011-02-02 18:55:23Z luw $    **
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
#include "md.h"

#define get_data(__PAR0__, __PAR1__, __PAR2__, __PAR4__) get_input( fhand, __PAR0__, __PAR1__, __PAR2__, __PAR4__)

static void read_atoms (FILE * fhand, char *tbuf0);
static void read_orbitals (FILE * fhand, char *tbuf0);
static void read_kpoints (FILE * fhand, char *tbuf);


void read_control (void)
{
    int tmp, mpi_nprocs;
    char *tbuf, *tptr;
    FILE *fhand;
    int i;

    if(pct.thispe == 0) 
    {
        printf("\n ********************************");
        printf("\n  control file read informations");
        printf("\n ********************************\n");
    }

    /* Open the input file for reading */
    my_fopen (fhand, ct.cfile, "r");
    fputc (' ', fhand);

    my_malloc (tptr, MAX_PATH, char);
    tbuf = tptr;

    /* Read in the description */
    get_data ("description", &ct.description, STR, "QMD run");

    /*   special for transport code */
    get_data ("metalic", &ct.metal, BOOL, "true");
    get_data ("num_blocks", &ct.num_blocks, INT, "3");

    /* read coarse grid info */
    get_data("blocks_dim", tbuf, STR, NULL);

    for(i=0; i <ct.num_blocks; i++)
    {
        ct.block_dim[i] = strtol(tbuf, &tbuf, 10);
    }

    get_data ("potential_compass", tbuf, STR, "NULL");
    potentialCompass.type = strtol(tbuf, &tbuf, 10);
    potentialCompass.box1.x1 = strtol(tbuf, &tbuf, 10);
    potentialCompass.box1.x2 = strtol(tbuf, &tbuf, 10);
    potentialCompass.box1.y1 = strtol(tbuf, &tbuf, 10);
    potentialCompass.box1.y2 = strtol(tbuf, &tbuf, 10);
    potentialCompass.box1.z1 = strtol(tbuf, &tbuf, 10);
    potentialCompass.box1.z2 = strtol(tbuf, &tbuf, 10);

    get_data ("chargedensity_compass", tbuf, STR, "NULL");
    chargeDensityCompass.type = strtol(tbuf, &tbuf, 10);
    chargeDensityCompass.box1.x1 = strtol(tbuf, &tbuf, 10);
    chargeDensityCompass.box1.x2 = strtol(tbuf, &tbuf, 10);
    chargeDensityCompass.box1.y1 = strtol(tbuf, &tbuf, 10);
    chargeDensityCompass.box1.y2 = strtol(tbuf, &tbuf, 10);
    chargeDensityCompass.box1.z1 = strtol(tbuf, &tbuf, 10);
    chargeDensityCompass.box1.z2 = strtol(tbuf, &tbuf, 10);

    /* read Leads potential window */
/*
    get_data("probe_potential_window", tbuf, STR, NULL);

    for(i=0; i <cei.num_probe; i++)
    {
        cei.probe_window_start[i] = strtol(tbuf, &tbuf, 10);
        cei.probe_window_end[i] = strtol(tbuf, &tbuf, 10);

          printf (" hello %d %d %d \n", cei.num_probe, 
          cei.probe_window_start[i], cei.probe_window_end[i]);
    }
*/

    /* Read in the starting wavefunction file name */
    get_data ("input_wave_function_file", &ct.infile, STR, "wave.out");

    /* Read in the output wavefunction file name */
    get_data ("output_wave_function_file", &ct.outfile, STR, "wave.out");

    /* Read in the initial run flag */
    get_data ("start_mode_NEGF", &ct.runflag, INT, "112");
    get_data ("average_plane_rho", tbuf, STR, "1  0  0  0  1");
    for(i=0; i < 5; i++)
    {
        ct.plane[i] = strtol(tbuf, &tbuf, 10);
    }
    if(pct.thispe == 0) 
    {
        printf(" \n average plane of rho ");
        for(i=0; i < 5; i++) printf(" %d ", ct.plane[i]);
        printf("\n");
    }


    /* Read in the boundary condition flag */
    get_data ("boundary_condition_type", &ct.boundaryflag, OPT, "Periodic");

    /* Read mixing parameter */
    get_data ("charge_density_mixing", &ct.mix, DBL, "0.5");

    /* Exchange correlation potential type flag */
    get_data ("exchange_correlation_type", &ct.xctype, OPT, "LDA");

    /* Number of scf steps */
    get_data ("max_scf_steps", &ct.max_scf_steps, INT, "100");
    ct.scfpermd = ct.max_scf_steps;

    /* RMS convergence criterion */
    get_data ("rms_convergence_criterion", &ct.thr_rms, DBL, "1.0E-7");

    /* force convergence criterion */
    get_data ("fast_relax_max_force", &ct.thr_frc, DBL, "2.5E-3");

    /* Write wavefunctions into output file, every md count of steps */
    get_data ("do_write_waves_to_file", &ct.checkpoint, BOOL, "true");
    if (ct.checkpoint)
    {
        get_data ("md_steps_til_write_waves", &ct.checkpoint, INT, "10");
        /*If steps_til_write_waves is 0 this is most likely an error
         * 1 was probably intended (i.e.) writing out every step*/
        if (ct.checkpoint == 0)
            ct.checkpoint = 1;
    }


    /* How often calculate energy, print out eigenvalues and occupancies  */
    get_data ("do_print_energy_and_eigenvalues", &ct.outcount, BOOL, "true");
    if (ct.outcount)
    {
        get_data ("scf_steps_til_energy_and_eigenv_print", &ct.outcount, INT, "1");
        /*If steps_til_energy_and_eigenv_print is 0 this is most likely an error
         * 1 was probably intended (i.e.) writing out every step*/
        if (ct.outcount == 0)
            ct.outcount = 1;
    }

    /* Fermi occupation flag */
    get_data ("occupations_type", &ct.occ_flag, OPT, NULL);

    /* Occupation width */
    get_data ("occupation_electron_temperature_eV", &ct.occ_width, DBL, "0.04");
    ct.occ_width *= eV_Ha;

    /* Occupation mixing */
    get_data ("occupation_number_mixing", &ct.occ_mix, DBL, "0.3");

    /* states occupancy count; overrides background charge */
    get_data("states_count_and_occupation", ct.occupation_str, STR, "");


    /* Force flag */
    get_data ("calculation_mode", &ct.forceflag, OPT, "Quench Electrons");
    if(pct.thispe == 0) printf("\n forceflag %d", ct.forceflag);

    /*maximum number of md steps */
    get_data ("max_md_steps", &ct.max_md_steps, INT, "10");

    /* Ionic timestep */
    get_data ("ionic_time_step", &ct.iondt, DBL, "50");


    /* -------------------------------- */
    /*                                  */
    /*   order-n code specific inputs   */
    /*                                  */
    /* -------------------------------- */

    /* read mg_eig method */
    get_data("mg_method", &ct.mg_method, OPT, "Pulay");

    /* read mg steps */
    get_data("mg_steps", &ct.mg_steps, INT, "2");

    /* read orbital movable centers option */
    get_data("do_movable_orbital_centers", &ct.movingCenter, BOOL, "false");

    /* read orbital movable centers steps */
    get_data("movable_orbital_centers_steps", &ct.movingSteps, INT, "40");






    /* Number of states */
    get_data ("number_of_orbitals", &ct.num_states, INT, "0");
    if (ct.num_states > MAX_STATES)
    {
        printf("\n increase MAX_STATES in params.h %d ", ct.num_states);
        error_handler("Too many states specified in input file");
    }


    /* Get k-points and weights */
    read_kpoints (fhand, tbuf);


    /*This is not read into any variable */
    get_data ("length_units", &tmp, OPT, "Bohr");

    /* Bravais lattice type */
    get_data ("bravais_lattice_type", &ct.ibrav, OPT, NULL);

    /* Lattice constants */
    /* These should be conditionally read depending on ibrav etc. */
    get_data ("a_length", &ct.celldm[0], DBL, NULL);
    get_data ("b_length", &ct.celldm[1], DBL, NULL);
    get_data ("c_length", &ct.celldm[2], DBL, NULL);
    get_data ("alpha", &ct.celldm[3], DBL, NULL);
    get_data ("beta", &ct.celldm[4], DBL, NULL);
    get_data ("gamma", &ct.celldm[5], DBL, NULL);

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

    get_data("potential_grid_refinement",  &RHO_NX, INT, NULL);

    RHO_NY = RHO_NX;
    RHO_NZ = RHO_NX; 

    FNX_GRID = NX_GRID * RHO_NX;
    FNY_GRID = NY_GRID * RHO_NY;
    FNZ_GRID = NZ_GRID * RHO_NZ;


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

    get_data ("kpoints_per_processor", &pct.pe_kpoint, INT, "1");
    get_data ("Hamiltonia_processor_grid", tbuf, STR, "1 1");
    NPROW = strtol(tbuf, &tbuf, 10);
    NPCOL = strtol(tbuf, &tbuf, 10);

    if(NPES != pct.pe_x * pct.pe_y * pct.pe_z ) 
    {
        printf("\n NPES = %d", NPES);
        printf("\n pct.pe_x, y,z = %d %d %d",pct.pe_x, pct.pe_y, pct.pe_z);
        error_handler("bad decomposion of processor grid");
    }

    if(NPES < pct.pe_kpoint * NPROW * NPCOL)
    {
        printf("\n NPES = %d", NPES);
        printf("\n pct.pe_kpoint, NPROW, NPCOL = %d %d %d",pct.pe_kpoint, NPROW, NPCOL);
        error_handler("bad decomposion of processor grid");
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


    FPX0_GRID = PX0_GRID * RHO_NX;
    FPY0_GRID = PY0_GRID * RHO_NY;
    FPZ0_GRID = PZ0_GRID * RHO_NZ;
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
    ct.qcparm = ct.cparm / (REAL) RHO_NX;
    ct.betacparm = ct.cparm / (REAL) ct.nxfgrid;
    ct.cparm /= (REAL) RHO_NX;


    /* read the atomic positions and species */
    read_atoms (fhand, tbuf);
    read_orbitals (fhand, tbuf);


    /* Close input file */
    my_free (tptr);
    fclose (fhand);

    potentialCompass.box1.x1 *= RHO_NX;
    potentialCompass.box1.x2 *= RHO_NX;
    potentialCompass.box1.y1 *= RHO_NY;
    potentialCompass.box1.y2 *= RHO_NY;
    potentialCompass.box1.z1 *= RHO_NZ;
    potentialCompass.box1.z2 *= RHO_NZ;

    chargeDensityCompass.box1.x1 *= RHO_NX;
    chargeDensityCompass.box1.x2 *= RHO_NX;
    chargeDensityCompass.box1.y1 *= RHO_NY;
    chargeDensityCompass.box1.y2 *= RHO_NY;
    chargeDensityCompass.box1.z1 *= RHO_NZ;
    chargeDensityCompass.box1.z2 *= RHO_NZ;


}                               /* end read_control */



static void read_kpoints (FILE * fhand, char *tbuf)
{
    int ik, nk;
    double w;

    /* first let's count the number of k-points */
    nk = 0;
    while (get_data ("kpoints", tbuf, LIST, NULL))
        nk++;

    my_malloc (ct.kp, nk, KPOINT);

    /* now we can read the kpoint data */
    ik = 0;
    while (get_data ("kpoints", tbuf, LIST, NULL))
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


static void read_atoms (FILE * fhand, char *tbuf0)
{
    int ni, is, ns;
    int ist, st1;
    char *tbuf = tbuf0;
    char pseudo_symbol[MAX_SPECIES][32];

    /* Absolute or cell relative coordinates */
    get_data ("atomic_coordinate_type", &ct.crd_flag, OPT, "Absolute");


    /*Count number of ions in input file*/
#if 1
    ct.num_ions = listlen(fhand, "atoms");


#else
    ni = 0;
    while (get_data ("atoms", tbuf, LIST, NULL))
        ni++;
    ct.num_ions = ni;
#endif


    /*Allocate memory for ions*/
    my_calloc(ct.ions, ct.num_ions, ION);



    /* read and count coordinates for each ion */
    ni = 0;
    ns = 0;
    ist = 0;
    while (get_data ("atoms", tbuf, LIST, NULL))
    {
        char s[32];
        int args = 0;
        args = sscanf (tbuf, "%s %lf %lf %lf %d %d %d",
                s,
                &ct.ions[ni].crds[0], &ct.ions[ni].crds[1], &ct.ions[ni].crds[2], 
                &ct.ions[ni].movable, &ct.ions[ni].frozen, &ct.ions[ni].n_loc_states
                );

        if (args < 7)
        {
            printf ("Error reading ion %d args %d\n", ni, args);
            error_handler ("Not enough arguments in atoms list!");
        }

        /* search for the species among existent pseudo_symbols */
        is = 0;
        while (is < ns && strcmp (s, pseudo_symbol[is]))
            is++;

        if (is == ns)           /* it's a new species */
            strcpy (pseudo_symbol[ns++], s);

        ct.ions[ni].species = is;


        /* Handle lattice and angstrom units */
        if (verify ("atomic_coordinate_type", "Cell Relative"))
        {
            ct.ions[ni].xtal[0] = ct.ions[ni].crds[0];
            ct.ions[ni].xtal[1] = ct.ions[ni].crds[1];
            ct.ions[ni].xtal[2] = ct.ions[ni].crds[2];
        }
        else if (verify ("length_units", "Angstrom"))
        {
            ct.ions[ni].crds[0] *= A_a0;
            ct.ions[ni].crds[1] *= A_a0;
            ct.ions[ni].crds[2] *= A_a0;
        }

        for(st1 = ist; st1 < ist+ct.ions[ni].n_loc_states; st1++)
        {
            state_to_ion[st1] = ni;
            states[st1].atomic_orbital_index = st1-ist;
        }

        ist += ct.ions[ni].n_loc_states;


        ni++;

    }                           /* end while (get_data( "atoms", tbuf, LIST, NULL) != 0) */


    ct.num_species = ns;


    /* allocate data structures for the species and setup the sp[j].pseudo_symbol */
    my_malloc (ct.sp, ns, SPECIES);

    for (is = 0; is < ns; is++)
    {
        strcpy (ct.sp[is].pseudo_symbol, pseudo_symbol[is]);
        sprintf (ct.sp[is].pseudo_filename, "%s.pp", pseudo_symbol[is]);
    }




    /* read info about pseudopotential filenames for each species */
    tbuf = tbuf0;
    while ((get_data ("pseudopotential", tbuf, STR | SEQ, NULL)) != 0)
    {
        char s[32], fn[MAX_PATH];
        if (sscanf (tbuf, " %s %s ", s, fn))
        {
            /* search for the species among the known species */
            is = 0;
            while (is < ns && strcmp (s, ct.sp[is].pseudo_symbol))
                is++;

            if (is < ns)        /* we've found it */
                strcpy (ct.sp[is].pseudo_filename, fn);
        }
    }



    /* read info about atomic orbital filenames for each species for LCAO start */
    tbuf = tbuf0;
    while ((get_data ("atomic_orbital_files", tbuf, STR | SEQ, NULL)) != 0)
    {
        char s[32], fn[MAX_PATH];
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







    /* optional number of atoms consistency check */
    if (get_data ("number_of_atoms", &ni, INT, "0"))
        if (ni != ct.num_ions)
        {
            printf ("number_of_atoms = %d != %d = atoms count\n", ni, ct.num_ions);
            error_handler ("Mismatch between number_of_atoms and atoms count");
        }

    /* optional number of species consistency check */
    if (get_data ("number_of_species", &ns, INT, "0"))
        if (ns != ct.num_species)
        {
            printf ("number_of_species = %d != %d = species count\n", ns, ct.num_species);
            error_handler ("Mismatch between number_of_species and species count");
        }




}


static void read_orbitals (FILE * fhand, char *tbuf0)
{
    int ni, is, ns;
    int i, ist, st1;
    int num_lines;
    double crds[3], radius;
    int num_tem, movable, frozen;
    char *tbuf = tbuf0;


    /*Count number of ions in input file*/

    num_lines = listlen(fhand, "orbitals");

    /* read and count coordinates for each ion */
    ni = 0;
    ns = 0;
    ist = 0;
    while (get_data ("orbitals", tbuf, LIST, NULL))
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
            states[ni + i].crds[0] = crds[0];
            states[ni + i].crds[1] = crds[1];
            states[ni + i].crds[2] = crds[2];
            states[ni + i].radius = radius;
            states[ni + i].movable = movable;
            states[ni + i].frozen = frozen;
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


