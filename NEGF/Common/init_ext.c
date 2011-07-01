/************************** SVN Revision Information **************************
 **    $Id: init_nuc.c 1434 2011-05-19 15:39:10Z luw $    **
******************************************************************************/
 
/*


    init_ext.c

    Set up gate voltage


*/



#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "md.h"


void init_ext (double *vext, int gbias_begin, int gbias_end, int reach_begin, int reach_end, double gate_bias)
{

    int ii, jj, kk, ix, iy, iz;
    int x_proc;  
    double v_external;    
    double x_locate;

    if (pct.gridpe == 0)
    {

        printf(" Begin init_ext ...\n");

    }                           /* end if */

    fflush(NULL);

    my_barrier();


    /* Grab some memory for storage */
    /*Bikan Tan */
    my_calloc( vext, FP0_BASIS, double );



    pe2xyz(pct.gridpe, &ii, &jj, &kk);
    x_proc = ii * FPX0_GRID;

    //    printf("\n my ii = %d and FPX0_GRID = %d  gate_bias = %f \n", ii, FPX0_GRID, gate_bias);   
    //    printf("\n gbias_begin = %d and gbias_end = %d \n", gbias_begin, gbias_end);   

    for (ix = 0; ix < FPX0_GRID; ix++)
    {

        x_locate = (ix + x_proc)/2.0; // everything is based on coarse grid now!
        v_external = gate_bias * (  2 / (1 + exp( -2.5*abs(x_locate - gbias_begin)/abs(reach_begin - gbias_begin) ) ) + 2 / (1 + exp( -2.5*abs(gbias_end - x_locate)/abs(reach_end - gbias_end ) ) ) - 3 );

        if (x_locate > gbias_begin && x_locate < gbias_end)  
        {

            for (iy = 0; iy < FPY0_GRID; iy++)
            {

                for (iz = 0; iz < FPZ0_GRID; iz++)
                {

                    vext[ix * FPY0_GRID * FPZ0_GRID + iy * FPZ0_GRID + iz] = v_external;


                }                   /* end for */

            }                       /* end for */

        }

        if (jj == 0 && kk == 0)
            printf("x_locate = %5f      vext  = %6.3f \n ", x_locate, v_external );

    }                           /* end for */



    if (pct.gridpe == 0)
    {

        printf(" init_ext done\n");

    }                           /* end if */

    /* Wait until everyone gets here */
    fflush(NULL);
    my_barrier();

}                               /* end init_ext */


