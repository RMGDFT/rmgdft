/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include "init_var.h"


void mat_dist_to_global (double *dist_mat, double *global_mat, int *desca) 
{

    int i, j,  k;
    int ictxt, mb, nprow, npcol, myrow, mycol;
    int jj, kk;


    ictxt = desca[1]; 
    mb = desca[4];

    Cblacs_gridinfo (ictxt, &nprow, &npcol, &myrow, &mycol);

    for(j = 0; j < ct.num_states * ct.num_states; j++)
        global_mat[j] = 0.0;


    for(k=0; k < MXLCOL; k++)
    {
        for(j =0; j < MXLLDA; j++)
        {


            jj = (j/mb ) * nprow * mb + myrow * mb + j - j/mb * mb; 
            kk = (k/mb ) * npcol * mb + mycol * mb + k - k/mb * mb; 

            if(jj < ct.num_states && kk < ct.num_states)
                global_mat[kk * ct.num_states + jj] = dist_mat[k * MXLLDA + j];

        }
    }

    int n2 = ct.num_states * ct.num_states;
    global_sums(global_mat, &n2, pct.grid_comm);

}

