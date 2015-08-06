/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "main.h"
#include "init_var.h"
#include "LCR.h"

#define         LDEBUG  0

void nlforce_par_D (STATE *states, double *forces)
{
    int st1, st2, ixyz;
    int nC, nL, i, ntot, ion;
    int idx1, idx2, size;
    int N_res, N_int, Num_ene_points;
    double *par_Hij_tri, *par_Sij_tri, *GHG_tri, *GHG_en_tri, *S_matrix;
    ION *iptr;


    /*allocate memory for par_Hij_tri, par_Sij_tri, GHG_tri and GHG_en_tri  */
    ntot = 0;
    for (i = 0; i < ct.num_blocks - 1; i++)
        ntot += (ct.block_dim[i] * ct.block_dim[i] + ct.block_dim[i] * ct.block_dim[i + 1]);
    ntot += ct.block_dim[ct.num_blocks - 1] * ct.block_dim[ct.num_blocks - 1];

    my_malloc_init( GHG_tri, ntot, double );
    my_malloc_init( GHG_en_tri, ntot, double );

    /* get the integration of multipliation of GHG and GHG*epslion over energy points */
    multi_GHG_munu (GHG_tri, GHG_en_tri);
    my_free (sigma_all);

    my_malloc_init( par_Hij_tri, ntot, double );
    my_malloc_init( par_Sij_tri, ntot, double );

    size = ct.num_states * ct.num_states;
    my_malloc_init( S_matrix, size, double );

    idx1 = ct.num_states - lcr[2].num_states;
    idx2 = ct.num_states - lcr[2].num_states / 2;
    for (ion = 0; ion < ct.num_ions; ion++)
    {
        if (ct.ions[ion].movable)
        {

            /* Making a loop for x , y, and z direction, corresponding to ixyz = 1, 2, 3 */
            for ( ixyz = 1; ixyz <= 3; ixyz++)
            {

                /* Calculating the part of partial_Hij/partial_x due to partial_beta and partial_Dnm*/
                ion_partial_Hij_and_Sij (ion, ixyz, work_matrix, S_matrix);
             
                for (st1 = 0; st1 < lcr[1].num_states / 2; st1++)
                    for (st2 = 0; st2 < lcr[1].num_states / 2; st2++)
                    {
                        work_matrix[st1 * ct.num_states + st2] = 0.0;
                        S_matrix[st1 * ct.num_states + st2] = 0.0;
                    }
                for (st1 = 0; st1 < lcr[1].num_states / 2; st1++)
                    for (st2 = 0; st2 < lcr[2].num_states / 2; st2++)
                    {
                        work_matrix[(st2 + idx2) * ct.num_states + st1] = 0.0;
                        work_matrix[st1 * ct.num_states + st2 + idx2] = 0.0;
                        S_matrix[(st2 + idx2) * ct.num_states + st1] = 0.0;
                        S_matrix[st1 * ct.num_states + st2 + idx2] = 0.0;
                    }
                for (st1 = lcr[2].num_states / 2; st1 < lcr[2].num_states; st1++)
                    for (st2 = lcr[2].num_states / 2; st2 < lcr[2].num_states; st2++)
                    {
             
                        work_matrix[(st1 + idx1) * ct.num_states + st2 + idx1] = 0.0;
                        S_matrix[(st1 + idx1) * ct.num_states + st2 + idx1] = 0.0;
                    }
             
                whole_to_tri_real (par_Hij_tri, work_matrix, ct.num_blocks, ct.block_dim);
                whole_to_tri_real (par_Sij_tri, S_matrix, ct.num_blocks, ct.block_dim);

                forces[ion*3 + ixyz -1] += 2.0 * trace_AB_tri( par_Hij_tri, GHG_tri, ct.num_blocks, ct.block_dim);
                forces[ion*3 + ixyz -1] -= 2.0 * trace_AB_tri( par_Sij_tri, GHG_en_tri, ct.num_blocks, ct.block_dim);
             
            }

        } /* end if */

    }   /* end for ion number */


    my_free(S_matrix);


    /* Calculating the force due to another part of partial Halmitonian (partial_Vnuc/partial_R) */
    tri_to_row(GHG_tri, work_matrix, ct.num_blocks, ct.block_dim);
    nlforce_partial_H_part2 (states, states1, work_matrix, forces);


    my_free(par_Hij_tri);
    my_free(par_Sij_tri);
    my_free(GHG_tri);
    my_free(GHG_en_tri);


}
