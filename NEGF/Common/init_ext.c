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
#include "main.h"
#include "init_var_negf.h"
#include "LCR.h"


void init_ext (double *vext, double gbias_begin, double gbias_end,  double BT, double gate_bias)
{

    int ix, iy, iz;
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



    x_proc = get_FPX_OFFSET();

    printf("\n gate_bias = %f \n", gate_bias);   
    printf("\n gbias_begin = %f and gbias_end = %f  BT = %f \n", gbias_begin, gbias_end, BT);   

    for (ix = 0; ix < get_FPX0_GRID(); ix++)
    {

        x_locate = (ix + x_proc)/2.0; // everything is based on coarse grid now!
        v_external = gate_bias / (1.0 + exp((gbias_begin - x_locate)/BT)  + exp((x_locate - gbias_end)/BT) ); 


            for (iy = 0; iy < get_FPY0_GRID(); iy++)
            {

                for (iz = 0; iz < get_FPZ0_GRID(); iz++)
                {

                        vext[ix * get_FPY0_GRID() * get_FPZ0_GRID() + iy * get_FPZ0_GRID() + iz] = v_external;


                }                   /* end for */

            }                       /* end for */


//            printf("x_locate = %5f      vext  = %10.7f \n ", x_locate, v_external );

    }                           /* end for */



    if (pct.gridpe == 0)
    {

        printf(" init_ext done\n");

    }                           /* end if */

    /* Wait until everyone gets here */
    fflush(NULL);
    my_barrier();

}                               /* end init_ext */


