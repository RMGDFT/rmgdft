#include "negf_prototypes.h"
/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <assert.h>

#include "main.h"
#include "init_var.h"
#include "LCR.h"

/*
 *
 *      calculate the trace of a distributed matrix in
 *     parallel matrix operations.
 *
 */

double pmo_trace(std::complex<double> *matrix, int *desca)
{

    int mycol, myrow, nprow, npcol;
    int ictxt = desca[1], mb = desca[4], nb = desca[5], mxllda = desca[8];
    int m = desca[2], n = desca[3];
    int i,j, ii,jj, icrow, iccol;
    int idx;
    double tem;


    tem = 0.0;

    Cblacs_gridinfo (ictxt, &nprow, &npcol, &myrow, &mycol);

    if(m !=n ) 
    {
        rmg_printf("\n m n %d %d \n", m, n);
        rmg_error_handler(__FILE__, __LINE__, "not a square matrix in pmo_trace m!=n");
    }

    for(i =0; i <m; i++)
    {

        j=i;

        icrow = (i/mb) % nprow;
        iccol = (j/nb) % npcol;
        /* global index (i,i) is in processor (icrow, iccol) */

        if(icrow == myrow && iccol == mycol)
        {

            ii= (i/mb) / nprow * mb + i%mb;
            jj= (j/nb) / npcol * nb + j%nb;

            idx = jj *mxllda + ii;

            tem +=std::real( matrix[idx]);
        }

    }

    return tem;
}


