/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include "init_var.h"


void mat_global_to_dist (double *global_mat, double *dist_mat, int *desca) 
{

    int i, j,  k;
    int ictxt, mb, nprow, npcol, myrow, mycol;
    int jj, kk;


    ictxt = desca[1]; 
    mb = desca[4];

    Cblacs_gridinfo (ictxt, &nprow, &npcol, &myrow, &mycol);

    for(j = 0; j < MXLLDA * MXLCOL; j++)
        dist_mat[j] = 0.0;


    for(k=0; k < MXLCOL; k++)
    {
        for(j =0; j < MXLLDA; j++)
        {


            jj = (j/mb ) * nprow * mb + myrow * mb + j - j/mb * mb; 
            kk = (k/mb ) * npcol * mb + mycol * mb + k - k/mb * mb; 

            if(jj < ct.num_states && kk < ct.num_states)
                dist_mat[k * MXLLDA + j] = global_mat[kk * ct.num_states + jj] ;

        }
    }


}

