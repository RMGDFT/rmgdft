/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/read_LCR.c *****
 * NAME
 *   Ab initio O(n) real space code with 
 *   localized orbitals and multigrid acceleration
 *   Version: 3.0.0
 * COPYRIGHT
 *   Copyright (C) 2001  Wenchang Lu,
 *                       Jerzy Bernholc
 * FUNCTION
 *   void read_LCR()
 *   read all informations from control input file
 * INPUTS
 *   no explicit input
 * OUTPUT
 *   read the control parameters for the non-linear green function calculations.
 * PARENTS
 *   md.c
 * CHILDREN
 *   nothing
 * SOURCE
 */


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include "main.h"
#include "init_var.h"
#include "LCR.h"



/* Reads a line from the input file into tbuf skips blank lines and
   comment lines. Calls error and terminates if eof encountered. */
char *get_line (char *buf, FILE * fh);
char *get_num (char *str);


/* Reads and parses the input control file */
void read_LCR ()
{

    int iprobe;
    char *tptr;
    FILE *fhand;
    char tbuf[200];
    char newname[200];


    for (iprobe = 0; iprobe < cei.num_subsystem; iprobe++)
    {

        /* Open the input file for reading */
        sprintf (newname, "%s%s%d", pct.image_path[pct.thisimg], "LCR.dat", iprobe);
        my_fopen (fhand, newname, "r");

        /* Read in the starting wavefunction file name */
        strcpy (newname, get_line (tbuf, fhand));
        sprintf (lcr[iprobe].name, "%s%s", pct.image_path[pct.thisimg], newname);

        /* Read in the calculated lead potential/rho file name */
        strcpy (newname, get_line (tbuf, fhand));
        sprintf (lcr[iprobe].lead_name, "%s%s", pct.image_path[pct.thisimg], newname);

        lcr[iprobe].NX_GRID = atoi (get_line (tbuf, fhand));

        if (NULL == (tptr = get_num (tbuf)))
        {
            printf ("\n probe %d missing the get_NY_GRID()\n", iprobe);
            error_handler ("need get_NY_GRID()");
        }
        else
            lcr[iprobe].NY_GRID = atoi (tptr);

        if (NULL == (tptr = get_num (tptr)))
        {
            printf ("\n probe %d missing the get_NZ_GRID()\n", iprobe);
            error_handler ("need get_NZ_GRID()");
        }
        else
            lcr[iprobe].NZ_GRID = atoi (tptr);


        lcr[iprobe].x0 = atoi (get_line (tbuf, fhand));

        if (NULL == (tptr = get_num (tbuf)))
        {
            printf ("\n probe %d missing the y0\n", iprobe);
            error_handler ("need y0");
        }
        else
            lcr[iprobe].y0 = atoi (tptr);

        if (NULL == (tptr = get_num (tptr)))
        {
            printf ("\n probe %d missing the z0\n", iprobe);
            error_handler ("need z0 ");
        }
        else
            lcr[iprobe].z0 = atoi (tptr);

        lcr[iprobe].x1 = atoi (get_line (tbuf, fhand));

        if (NULL == (tptr = get_num (tbuf)))
        {
            printf ("\n probe %d missing the y1\n", iprobe);
            error_handler ("need y1 ");
        }
        else
            lcr[iprobe].y1 = atoi (tptr);

        if (NULL == (tptr = get_num (tptr)))
        {
            printf ("\n probe %d missing the z1\n", iprobe);
            error_handler ("need z1 ");
        }
        else
            lcr[iprobe].z1 = atoi (tptr);


        lcr[iprobe].x2 = atoi (get_line (tbuf, fhand));

        if (NULL == (tptr = get_num (tbuf)))
        {
            printf ("\n probe %d missing the y2\n", iprobe);
            error_handler ("need y2 ");
        }
        else
            lcr[iprobe].y2 = atoi (tptr);

        if (NULL == (tptr = get_num (tptr)))
        {
            printf ("\n probe %d missing the z2\n", iprobe);
            error_handler ("need z2 ");
        }
        else
            lcr[iprobe].z2 = atoi (tptr);


        lcr[iprobe].num_ions = atoi (get_line (tbuf, fhand));



        lcr[iprobe].state_begin = atoi (get_line (tbuf, fhand));

        if (NULL == (tptr = get_num (tbuf)))
        {
            printf ("\n probe %d missing the state_middle\n", iprobe);
            error_handler ("need state_middle ");
        }
        else
            lcr[iprobe].state_middle = atoi (tptr);

        if (NULL == (tptr = get_num (tptr)))
        {
            printf ("\n probe %d missing the state_end\n", iprobe);
            error_handler ("need state_end ");
        }
        else
            lcr[iprobe].state_end = atoi (tptr);

        lcr[iprobe].num_states = lcr[iprobe].state_end - lcr[iprobe].state_begin;


        lcr[iprobe].ion_begin = atoi (get_line (tbuf, fhand));

        lcr[iprobe].EF_new = atof (get_line (tbuf, fhand));

        if (NULL == (tptr = get_num (tbuf)))
        {
            printf ("\n probe %d missing the  EF_old\n", iprobe);
            error_handler ("need state_end ");
        }
        else
            lcr[iprobe].EF_old = atof (tptr);
        if (NULL == (tptr = get_num (tptr)))
        {
            printf ("\n probe %d missing the  bias\n", iprobe);
            error_handler ("need state_end ");
        }
        else
            lcr[iprobe].bias = atof (tptr);

        lcr[iprobe].xside = atof (get_line (tbuf, fhand));
        if (NULL == (tptr = get_num (tbuf)))
        {
            printf ("\n probe %d missing the  x_shift\n", iprobe);
            error_handler ("need x_shift ");
        }
        else
        {
            lcr[iprobe].x_shift = atof (tptr);
        }

        lcr[iprobe].yside = atof (get_line (tbuf, fhand));
        if (NULL == (tptr = get_num (tbuf)))
        {
            printf ("\n probe %d missing the  y_shift\n", iprobe);
            error_handler ("need y_shift ");
        }
        else
        {
            lcr[iprobe].y_shift = atof (tptr);
        }


    }


}

/********/
