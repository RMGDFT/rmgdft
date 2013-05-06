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

    read_common();
    my_malloc (tptr, MAX_PATH, char);
    tbuf = tptr;

   char start_mode_opts[] = "Random Start\n"
                            "Restart From File\n"
                            "FIREBALL Start\n"
                            "Gaussian Start";
    get_data ("start_mode", NULL, INIT | OPT, start_mode_opts);


    /* Read in the initial run flag */
    get_data ("start_mode", &ct.runflag, OPT, "Random Start");



    get_data ("freeze_orbital_step", &ct.freeze_orbital_step, INT, "90");



    get_data("bg_begin", &ct.bg_begin, DBL, "0.0");
    get_data("bg_end", &ct.bg_end, DBL, "10.0");
    get_data("BT", &ct.BT, DBL, "4.0");


    ct.constrainforces = 0;


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


    /* number of excess electrons in the system (useful for doped systems) */
    get_data("system_charge", &ct.background_charge, DBL, "0");

    /*Background charge is defined to be the opposite of system charge */
    ct.background_charge *= -1.0;



    get_data ("kpoints_per_processor", &pct.pe_kpoint, INT, "1");
    get_data ("Hamiltonia_processor_grid", tbuf, STR, "1 1");
    pct.nprow = strtol(tbuf, &tbuf, 10);
    pct.npcol = strtol(tbuf, &tbuf, 10);


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



    PX0_GRID = NX_GRID/pct.pe_x;
    PY0_GRID = NY_GRID/pct.pe_y;
    PZ0_GRID = NZ_GRID/pct.pe_z;
    P0_BASIS = PX0_GRID * PY0_GRID * PZ0_GRID;

    pct.PX0_GRID = PX0_GRID;
    pct.PY0_GRID = PY0_GRID;
    pct.PZ0_GRID = PZ0_GRID;
    pct.P0_BASIS = P0_BASIS;


    S0_BASIS = (PX0_GRID+2) * (PY0_GRID+2) * (PZ0_GRID+2);


    FPX0_GRID = PX0_GRID * FG_NX;
    FPY0_GRID = PY0_GRID * FG_NY;
    FPZ0_GRID = PZ0_GRID * FG_NZ;
    FP0_BASIS = FPX0_GRID * FPY0_GRID * FPZ0_GRID;

    pct.FPX0_GRID = FPX0_GRID;
    pct.FPY0_GRID = FPY0_GRID;
    pct.FPZ0_GRID = FPZ0_GRID;
    pct.FP0_BASIS = FP0_BASIS;


    /* Output some information for GW calculation. */
    get_data ("output_information_for_GW", &ct.flag_gw, INT, "0");


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


