/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*

init_nonlocal_comm.c

    Sets up the parameters for communication in calculation of
<psi_i|kb_ion> <kb_ion|psi_j>

each processor calculates the part of <psi|kb> in which psi are stored in
this processor.


This should be called after get_nlop.c
*/




#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "main_on.h"



void init_nonlocal_comm(void)
{
    int ion, idx, item;
    int size;


    my_malloc(num_nonlocal_ion, NPES, int);
    for (idx = 0; idx < NPES; idx++)
        num_nonlocal_ion[idx] = 0;
    num_nonlocal_ion[pct.gridpe] = pct.n_ion_center;
    global_sums_int(num_nonlocal_ion, &NPES);

    max_ion_nonlocal = 0;
    for (idx = 0; idx < NPES; idx++)
    {
        if (max_ion_nonlocal < num_nonlocal_ion[idx])
            max_ion_nonlocal = num_nonlocal_ion[idx];
    }

    /* ionidx_allproc is same as pct.ionidx[] except that this store those
       from all processors */


    if (ionidx_allproc != NULL)
        my_free(ionidx_allproc);
    my_malloc( ionidx_allproc, max_ion_nonlocal * NPES, int );
    for (idx = 0; idx < max_ion_nonlocal * NPES; idx++)
        ionidx_allproc[idx] = 0.0;

    for (ion = 0; ion < pct.n_ion_center; ion++)
        ionidx_allproc[pct.gridpe * max_ion_nonlocal + ion] = pct.ionidx[ion];

    item = max_ion_nonlocal * NPES;
    global_sums_int(ionidx_allproc, &item);

    size = ct.state_per_proc * max_ion_nonlocal * ct.max_nl;

    /* allocate mem for kbpsi <psi|kb>  */

    if (kbpsi != NULL)
    {
        my_free(kbpsi);
        my_free(kbpsi_comm);
        my_free(partial_kbpsi_x);
        my_free(partial_kbpsi_y);
        my_free(partial_kbpsi_z);
    }


    my_malloc_init( kbpsi, size, rmg_double_t );
    my_malloc_init( kbpsi_comm, size, rmg_double_t );
    my_malloc_init( partial_kbpsi_x, size, rmg_double_t );
    my_malloc_init( partial_kbpsi_y, size, rmg_double_t );
    my_malloc_init( partial_kbpsi_z, size, rmg_double_t );

}
