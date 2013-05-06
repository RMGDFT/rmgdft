/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/* routine to read the cond.input file */

#include <stdio.h>
#include <stdlib.h>
#include "main.h"


char *get_num (char *str);
char *get_line (char *buf, FILE * fh);


void read_cond_input (double *emin, double *emax, int *E_POINTS, double *E_imag, double *KT, int *kpoint)
{
    FILE *fhand;
    int idx;
    char tbuf[200], *tptr;

    /* Open the input file for reading */
    my_fopen (fhand, "cond.input", "r");

    ct.num_states = atoi (get_line (tbuf, fhand));

    if (NULL == (tptr = get_num (tbuf)))
    {
        printf (" missing the lcr[1].num_states \n");
        exit (0);
    }
    lcr[1].num_states = atoi (tptr);

    if (NULL == (tptr = get_num (tptr)))
    {
        printf (" missing the lcr[2].num_states \n");
        exit (0);
    }
    lcr[2].num_states = atoi (tptr);

    ct.num_blocks = atoi (get_line (tbuf, fhand));
    tptr = tbuf;
    for (idx = 0; idx < ct.num_blocks; idx++)
    {
        tptr = get_num (tptr);
        ct.block_dim[idx] = atoi (tptr);
    }

    *emin = atof (get_line (tbuf, fhand));
    *emax = atof (tptr = get_num (tbuf));
    *E_POINTS = atoi (get_num (tptr));

    *E_imag = atof (get_line (tbuf, fhand));
    *KT = atof (get_line (tbuf, fhand));
    kpoint[0] = atoi (get_line (tbuf, fhand));
    kpoint[1] = atoi (tptr = get_num (tbuf));
    kpoint[2] = atoi (tptr = get_num (tptr));

    ct.num_cond_curve = atoi (get_line (tbuf, fhand));
    my_malloc_init( ct.cond_probe1, ct.num_cond_curve, int ); 
    my_malloc_init( ct.cond_probe2, ct.num_cond_curve, int ); 
    for (idx = 0; idx < ct.num_cond_curve; idx++)
    {
        ct.cond_probe1[idx] = atof (get_line (tbuf, fhand));
        ct.cond_probe2[idx] = atof (tptr = get_num (tbuf));
    }
/*
    ct.cond_probe1 = atof (get_line (tbuf, fhand));
    ct.cond_probe2 = atof (tptr = get_num (tbuf));
*/

}
