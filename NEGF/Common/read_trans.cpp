#include "negf_prototypes.h"
/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include "pmo.h"

char *get_num (char *str);
char *get_line (char *buf, FILE * fh);

/* Read the parameters for NEGF calculations */

void read_trans (complex_energy_integral * cei)
{

    char *tptr;
    FILE *fhand;
    char tbuf[200];
    int idx;
	
    /* Open the input file for reading */
    sprintf (tbuf, "%s%s", pct.image_path[pct.thisimg], "trans.in");
    my_fopen (fhand, tbuf, "r");

    /* Read in the initial run flag */
    cei->num_probe = atoi (get_line (tbuf, fhand));
    assert(cei->num_probe <= NUM_PROBE_MAX);
    tptr = tbuf;
	my_malloc_init( cei->probe_in_block, cei->num_probe, int );
    for (idx = 0; idx < cei->num_probe; idx++)
    {
        tptr = get_num (tptr);
        cei->probe_in_block[idx] = atoi (tptr);
    }

    /* Read in number of subsystem and their order */
    cei->num_subsystem = atoi (get_line (tbuf, fhand));
    assert(cei->num_subsystem < NUM_SUBSYSTEM_MAX);
    tptr = tbuf;
	my_malloc_init( cei->subsystem_idx, cei->num_subsystem, int );
    for (idx = 0; idx < cei->num_subsystem; idx++)
    {
        tptr = get_num (tptr);
        cei->subsystem_idx[idx] = atoi (tptr);
    }


    /* read Leads potential window */

    cei->num_probe_window = atoi (get_line (tbuf, fhand));
    tptr = tbuf;
    my_malloc_init( cei->probe_window_start, cei->num_probe_window, int );
    my_malloc_init( cei->probe_window_end, cei->num_probe_window, int );

    for (idx = 0; idx < cei->num_probe_window; idx++)
    {
        tptr = get_num (tptr);
        cei->probe_window_start[idx] = atof (tptr);
        tptr = get_num (tptr);
        cei->probe_window_end[idx] = atof (tptr);
    }


    /* read dos_window for integration */

    cei->num_dos_window = atoi (get_line (tbuf, fhand));
    tptr = tbuf;
    my_malloc_init( cei->dos_window_start, cei->num_dos_window, int );
    my_malloc_init( cei->dos_window_end, cei->num_dos_window, int );

    for (idx = 0; idx < cei->num_dos_window; idx++)
    {
        tptr = get_num (tptr);
        cei->dos_window_start[idx] = atof (tptr);
        tptr = get_num (tptr);
        cei->dos_window_end[idx] = atof (tptr);
    }


/*
    ct.num_cond_curve = atoi (get_line (tbuf, fhand));
    my_malloc_init( ct.cond_probe1, ct.num_cond_curve, int );
    my_malloc_init( ct.cond_probe2, ct.num_cond_curve, int );
    for (idx = 0; idx < ct.num_cond_curve; idx++)
    {
        ct.cond_probe1[idx] = atof (get_line (tbuf, fhand));
        ct.cond_probe2[idx] = atof (tptr = get_num (tbuf));
    }

*/

	cei->ncircle = atoi (get_line (tbuf, fhand));
    cei->nmax_gq1 = atoi (get_line (tbuf, fhand));
    cei->nmax_gq2 = atoi (get_line (tbuf, fhand));

    cei->EB = atof (get_line (tbuf, fhand));
    cei->KT = atof (get_line (tbuf, fhand));
    cei->GAMMA = atof (get_line (tbuf, fhand));
    cei->DELTA2 = atof (get_line (tbuf, fhand));
    cei->DELTA = atof (get_line (tbuf, fhand));
    cei->Npulaysave = atoi (get_line (tbuf, fhand));
    cei->Npulayrefresh = atoi (get_line (tbuf, fhand));
    cei->pulaymix = atof (get_line (tbuf, fhand));

    pmo.nrow = atoi (get_line (tbuf, fhand));
    pmo.ncol = atoi (get_line (tbuf, fhand));

#if CUDA_ENABLED
    pmo.nrow = 1;
    pmo.ncol = 1;
#endif


    if (pct.gridpe == 0)
    {
        printf ("\n  PARAMETERS control the complex energy integral  \n");

        printf (" number of probe:                     %d \n", cei->num_probe);
        for (idx = 0; idx < cei->num_probe; idx++)
        {
           printf (" probe in which block?:   %d \n", cei->probe_in_block[idx]);
        }
        printf (" number of subsystem:                     %d \n", cei->num_subsystem);
        for (idx = 0; idx < cei->num_subsystem; idx++)
        {
           printf (" Order of subsystem in the input file?:   %d \n", cei->subsystem_idx[idx]);
        }
        printf (" energy point in semicircle:          %d \n", cei->ncircle);
        printf (" energy point in strait line:         %d \n", cei->nmax_gq1);
        printf (" energy point in noneq calc:          %d \n", cei->nmax_gq2);
        printf (" bottom of valance band :             %f eV \n", cei->EB);
        printf (" Fermi Dirac temperature :            %f eV \n", cei->KT);
        printf (" parameter GAMMA         :            %f eV \n", cei->GAMMA);
        printf (" imaginary energy for strait line:    %f eV \n", cei->DELTA2);
        printf (" imaginary energy for noneq calc :    %f eV \n", cei->DELTA);
        printf (" pulay steps in rho:                  %d \n", cei->Npulaysave);
        printf (" pulay refresh steps in rho:          %d \n", cei->Npulayrefresh);
        printf (" pulay mixing parameter in rho:       %f \n", cei->pulaymix);

        printf (" \n parallel matrix grid nrow, ncol:  %d %d \n", pmo.nrow, pmo.ncol);



    }

    fclose (fhand);
}
