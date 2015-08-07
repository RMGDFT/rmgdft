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
#include "main.h"
#include "prototypes_on.h"
#include "init_var.h"



void init_nonlocal_comm(void)
{
    int ion, idx, item;
    int size;
    int pair_find;
    int ion1, ion2, ion_global1, ion_global2;


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
        my_free(kbpsi_res);
        my_free(partial_kbpsi_x);
        my_free(partial_kbpsi_y);
        my_free(partial_kbpsi_z);
    }


    my_malloc_init( kbpsi, size, double );
    my_malloc_init( kbpsi_comm, size, double );
    my_malloc_init( kbpsi_res, size, double );
    my_malloc_init( partial_kbpsi_x, size, double );
    my_malloc_init( partial_kbpsi_y, size, double );
    my_malloc_init( partial_kbpsi_z, size, double );


    int proc1, proc2, proc3;
    int loop;
    int *matrix_pairs, *proc_mark;

    my_calloc( matrix_pairs, NPES * NPES, int );
    my_calloc( proc_mark, NPES, int );

    my_calloc( kbpsi_comm_send, NPES, int );
    my_calloc( kbpsi_comm_recv, NPES, int );


    for (proc1 = 0; proc1 < NPES * NPES; proc1++)
        matrix_pairs[proc1] = 0;

    proc1 = pct.gridpe;
    for (proc2 = proc1 + 1; proc2 < NPES; proc2++)
    {
        for (ion1 = 0; ion1 <num_nonlocal_ion[proc1]; ion1++)
            for (ion2 = 0; ion2 <num_nonlocal_ion[proc2]; ion2++)
            {
                ion_global1 = ionidx_allproc[proc1 * max_ion_nonlocal + ion1] ;
                ion_global2 = ionidx_allproc[proc2 * max_ion_nonlocal + ion2] ;
                if(ion_global1 == ion_global2) 
                {
                    matrix_pairs[proc1 * NPES + proc2] = 1;
                    matrix_pairs[proc2 * NPES + proc1] = 1;
                    break;
                }
            }
    }

    item = NPES * NPES;
    global_sums_int(matrix_pairs, &item);


    //   if (pct.gridpe == 0)
    //   {
    //       printf("\n initial communication matrix ");
    //       for (i = 0; i < NPES; i++)
    //       {
    //           printf("\n");
    //           for (j = 0; j < NPES; j++)
    //               printf(" %d ", matrix_pairs[i * NPES + j]);
    //       }
    //   }

    kbpsi_num_loop = 0;
    for (loop = 0; loop < NPES; loop++)
    {
        kbpsi_comm_send[loop] = -1;
        kbpsi_comm_recv[loop] = -1;
        for (proc1 = 0; proc1 < NPES; proc1++)
            proc_mark[proc1] = 1;
        pair_find = 0;
        for (proc1 = 0; proc1 < NPES; proc1++)
        {
            for (proc3 = proc1+1; proc3 < NPES + proc1 + 1; proc3++)
            {
                proc2 = proc3%NPES;
                if (proc_mark[proc2] && matrix_pairs[proc1 * NPES + proc2]) 
                {
                    pair_find++;
                    if(pct.gridpe == proc1) kbpsi_comm_send[loop] = proc2;
                    if(pct.gridpe == proc2) kbpsi_comm_recv[loop] = proc1;
                    matrix_pairs[proc1 * NPES + proc2]=0;
                    proc_mark[proc2] = 0;
                    break;
                }
            } 
        }

        if (pair_find == 0)
        {
            kbpsi_num_loop = loop;
            break;
        }
    }


    dprintf("\n kbpsi_num_loop %d", kbpsi_num_loop);
    //    for (loop = 0; loop < NPES; loop++)
    //      dprintf("\n\n %d  %d  %d  %d   loooop\n\n", pct.gridpe, loop,
    //kbpsi_comm_send[loop], kbpsi_comm_recv[loop]);

    //    for (loop = 0; loop < num_loop_kbpsi; loop++)
    //    {
    //        dprintf("\nLoop: %d  PE:%d send %d  ", loop, pct.gridpe, comm_pair[loop]);
    //
    //    }

    my_free(matrix_pairs);
    my_free(proc_mark);

}
