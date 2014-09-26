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
#include "prototypes_on.h"
#include "init_var.h"
static void read_orbitals ();


void read_control (char *file)
{
    char *tbuf, *tptr;
    int is,ns;


    /* Open the input file for reading */

    get_data (file, NULL, INIT | TAGS, NULL);

    read_common();
    my_malloc (tptr, MAX_PATH, char);
    tbuf = tptr;

   char start_mode_opts[] = "Random Start\n"
                            "Restart From File\n"
                            "FIREBALL Start\n"
                            "Gaussian Start\n"
                            "Restart TDDFT";
    get_data ("start_mode", NULL, INIT | OPT, start_mode_opts);


    /* Read in the initial run flag */
    get_data ("start_mode", &ct.runflag, OPT, "Random Start");



    get_data ("freeze_orbital_step", &ct.freeze_orbital_step, INT, "90");



    get_data("bg_begin", &ct.bg_begin, DBL, "-100");
    get_data("bg_end", &ct.bg_end, DBL, "622");
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

    /* number of waves to plot, by default plot zero waves */
    get_data("num_waves_to_plot", &ct.num_waves, INT, "0");

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
    get_data ("number_of_atoms", &ct.num_ions, INT, "0");


    allocate_states();
    get_state_to_proc(states);




    /* number of excess electrons in the system (useful for doped systems) */
    get_data("system_charge", &ct.background_charge, DBL, "0");

    /*Background charge is defined to be the opposite of system charge */
    ct.background_charge *= -1.0;



    get_data ("kpoints_per_processor", &pct.pe_kpoint, INT, "1");
    get_data ("Hamiltonia_processor_grid", tbuf, STR, "1 1");
    pct.scalapack_nprow = strtol(tbuf, &tbuf, 10);
    pct.scalapack_npcol = strtol(tbuf, &tbuf, 10);




    if(NPES < pct.pe_kpoint * pct.scalapack_nprow * pct.scalapack_npcol)
    {
        printf("\n NPES = %d", NPES);
      printf("\n pct.pe_kpoint, pct.scalapack_nprow, pct.scalapack_npcol = %d %d %d",pct.pe_kpoint, pct.scalapack_nprow, pct.scalapack_npcol);
        error_handler("bad decomposion of processor grid");
    }

    if(pct.scalapack_nprow > pct.scalapack_npcol ) 
    {
        printf("\n pct.scalapack_nprow, pct.scalapack_npcol = %d %d ", pct.scalapack_nprow, pct.scalapack_npcol);
        error_handler("pct.scalapack_nprow should be smaller than pct.scalapack_npcol");
    }



//    get_PX0_GRID() = get_NX_GRID()/pct.pe_x;
//    get_PY0_GRID() = get_NY_GRID()/pct.pe_y;
//    get_PZ0_GRID() = get_NZ_GRID()/pct.pe_z;
//    get_P0_BASIS() = get_PX0_GRID() * get_PY0_GRID() * get_PZ0_GRID();


    S0_BASIS = (get_PX0_GRID()+2) * (get_PY0_GRID()+2) * (get_PZ0_GRID()+2);


//    get_FPX0_GRID() = get_PX0_GRID() * get_FG_NX();
//    get_FPY0_GRID() = get_PY0_GRID() * get_FG_NY();
//    get_FPZ0_GRID() = get_PZ0_GRID() * get_FG_NZ();
//    get_FP0_BASIS() = get_FPX0_GRID() * get_FPY0_GRID() * get_FPZ0_GRID();


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



static void read_orbitals ()
{
    int ni;
    int i;
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


