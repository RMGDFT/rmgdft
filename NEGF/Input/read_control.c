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
#include "init_var.h"
#include "LCR.h"
#include "twoParts.h"

static void read_orbitals ();
static void read_kpoints ();


void read_control (char *file)
{
    int tmp, mpi_nprocs;
    char *tbuf, *tptr;
    FILE *fhand;
    int is,ns;
    int i;

    /* Open the input file for reading */

    get_data (file, NULL, INIT | TAGS, NULL);

    read_common();
    my_malloc (tptr, MAX_PATH, char);
    tbuf = tptr;

   char start_mode_opts[] = "Random Start\n"
                            "Restart From File\n"
                            "LCAO Start";
    get_data ("start_mode", NULL, INIT | OPT, start_mode_opts);


    /* Read in the initial run flag */
    get_data ("start_mode", &ct.runflag, OPT, "Random Start");





    get_data ("Simpson_depth", &ct.simpson_depth, INT, "0");
    get_data ("Simpson_tol", &ct.simpson_tol, DBL, "0.001");

    get_data ("bg_begin", &ct.bg_begin, DBL, "-100");
    get_data ("bg_end",   &ct.bg_end, DBL, "622");

    get_data ("gbias_begin", &ct.gbias_begin, DBL, "-100");
    get_data ("gbias_end",   &ct.gbias_end, DBL, "622");
    get_data ("BT", &ct.BT, DBL, "4.0");
    get_data ("gate_bias", &ct.gate_bias, DBL, "0.0");

//  vcomp_Lbegin vcomp_Lend controls the left buffer region in center part which must be aligned with left lead
    get_data ("vcomp_Lbegin", &ct.vcomp_Lbegin, INT, "-10");
    get_data ("vcomp_Lend",   &ct.vcomp_Lend, INT, "-9");

//  vcomp_Rbegin vcomp_Rend controls the right buffer region in center part which must be aligned with right lead
    get_data ("vcomp_Rbegin", &ct.vcomp_Rbegin, INT, "-8");
    get_data ("vcomp_Rend", &ct.vcomp_Rend, INT, "-7");
    
//  by default, we will not print 3Ddos for transmission peaks 
    get_data ("auto_3Ddos", &ct.auto_3Ddos, INT, "0");

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


    MPI_Comm_size (MPI_COMM_WORLD, &mpi_nprocs);
    if(NPES != mpi_nprocs) 
    {
        printf("\n NPES, mpi_nproc %d %d", NPES, mpi_nprocs);
        error_handler("bad NPES: job and input not match");
    }



    /* read the atomic positions and species */
    read_atoms ();

    read_orbitals ();


/**************************************************
 ******   NEGF specific input  ********************
***************************************************/
    get_data ("metalic", &ct.metal, BOOL, "false");
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

    /* Read in the initial run flag */
    get_data ("start_mode_NEGF", &ct.runflag, INT, "112");
    get_data ("average_plane_rho", tbuf, STR, "1  0  0  0  1");
    for(i=0; i < 5; i++)
    {
        ct.plane[i] = strtol(tbuf, &tbuf, 10);
    }
    if(pct.gridpe == 0)
    {
        printf(" \n average plane of rho ");
        for(i=0; i < 5; i++) printf(" %d ", ct.plane[i]);
        printf("\n");
    }



    potentialCompass.box1.x1 *= get_FG_RATIO();
    potentialCompass.box1.x2 *= get_FG_RATIO();
    potentialCompass.box1.y1 *= get_FG_RATIO();
    potentialCompass.box1.y2 *= get_FG_RATIO();
    potentialCompass.box1.z1 *= get_FG_RATIO();
    potentialCompass.box1.z2 *= get_FG_RATIO();

    chargeDensityCompass.box1.x1 *= get_FG_RATIO();
    chargeDensityCompass.box1.x2 *= get_FG_RATIO();
    chargeDensityCompass.box1.y1 *= get_FG_RATIO();
    chargeDensityCompass.box1.y2 *= get_FG_RATIO();
    chargeDensityCompass.box1.z1 *= get_FG_RATIO();
    chargeDensityCompass.box1.z2 *= get_FG_RATIO();



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
    int num_tem, movable, frozen;
    char tbuf [MAX_PATH];


    get_data("orbitals", &num_lines, INIT | LIST, NULL);

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


