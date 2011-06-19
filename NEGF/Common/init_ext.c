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


void init_ext (double *vext, int x_begin, int x_end, double level_diff)
{

    int ii, jj, kk, ix, iy, iz;
    int x_proc, x_locate;  
    

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
    
    for (ix = 0; ix < FPX0_GRID; ix++)
    {

	    for (iy = 0; iy < FPY0_GRID; iy++)
	    {


		    for (iz = 0; iz < FPZ0_GRID; iz++)
		    {

			    x_locate = ix + x_proc;
			    if (x_locate > x_begin && x_locate < x_end)  
				    vext[ix * FPY0_GRID * FPZ0_GRID + iy * FPZ0_GRID + iz] = level_diff * (2 / (1 + exp(-100*(x_locate-x_begin)) ) + 2 / (1 + exp(-100*(x_end-x_locate)) ) - 3 );


		    }                   /* end for */

	    }                       /* end for */

	    printf("x_locate = %5d      vext  = %6.3f \n ", x_locate, vext[ix * FPY0_GRID * FPZ0_GRID + iy * FPZ0_GRID + iz] );

    }                           /* end for */



    if (pct.gridpe == 0)
    {

        printf(" init_ext done\n");

    }                           /* end if */

    /* Wait until everyone gets here */
    fflush(NULL);
    my_barrier();

}                               /* end init_ext */


