/************************** SVN Revision Information **************************
 **    $Id: init_nuc.c 1434 2011-05-19 15:39:10Z luw $    **
******************************************************************************/
 
/*
    init_comp.c

    Set up compensating voltage *vcomp to compensate the misalignment caused by
    different background charge in lead calculation and center part calculation

*/



#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"


void init_comp (double *vcomp)
{

    int ix, iy, iz;
    int x_proc;  
    double v_comp;    
    double x_locate;
    double delta_v1;
    double delta_v2;
    double delta_v3;
    double delta_v4;
    int lb = ct.vcomp_Lbegin * 2;  //lb,le,rb and re is now base on FINE grid
    int le = ct.vcomp_Lend * 2;
    int rb = ct.vcomp_Rbegin * 2;
    int re = ct.vcomp_Rend * 2;

    if (pct.gridpe == 0)
    {

        printf(" Begin init_comp ...\n");

    }                           /* end if */

    fflush(NULL);

    my_barrier();

    x_proc = pct.FPX_OFFSET;

    printf(" vcomp_Lbegin = %d  and vcomp_Lend = %d  on FINE GRID \n", lb, le);   
    printf(" vcomp_Rbegin = %d  and vcomp_Rend = %d  on FINE GRID \n", rb, re);   
    printf(" vcomp_Lbegin = %d  and vcomp_Lend = %d  on coarse GRID \n", ct.vcomp_Lbegin, ct.vcomp_Lend);   
    printf(" vcomp_Rbegin = %d  and vcomp_Rend = %d  on coarse GRID \n", ct.vcomp_Rbegin, ct.vcomp_Rend);   
   
    delta_v1 = zvec[lb - 1] - zvec[lb];
    delta_v2 = zvec[lb - 1] - zvec[le];
    delta_v3 = zvec[re] - zvec[rb];
    delta_v4 = zvec[re] - zvec[re - 1];
    printf(" delta_v1  = %f \n", delta_v1);   
    printf(" delta_v2  = %f \n", delta_v2);   
    printf(" delta_v3  = %f \n", delta_v3);   
    printf(" delta_v4  = %f \n", delta_v4);   

    for (ix = 0; ix < FPX0_GRID; ix++)
    {
        x_locate = ix + x_proc; // everything is based on FINE grid!

        if ( x_locate >= lb && x_locate <= le )
        
            v_comp = (delta_v2 - delta_v1) * (x_locate - lb) / (le - lb) + delta_v1; 
        
        else if (x_locate > rb && x_locate < re)

            v_comp = (delta_v4 - delta_v3) * (x_locate - rb)/ (re - rb) + delta_v3;

        else if (x_locate > le && x_locate <= rb)

            v_comp = (delta_v3 - delta_v2) * (x_locate - le)/ (rb - le) + delta_v2;

        else

            v_comp = 0;


        for (iy = 0; iy < FPY0_GRID; iy++)
        {

            for (iz = 0; iz < FPZ0_GRID; iz++)
            {

                vcomp[ix * FPY0_GRID * FPZ0_GRID + iy * FPZ0_GRID + iz] = v_comp;


            }                   /* end for */

        }                       /* end for */


	printf("x_locate = %5f      vcomp  = %10.7f \n ", x_locate, v_comp );

    }                           /* end for */



    if (pct.gridpe == 0)
    {

        printf(" init_comp done\n");

    }                           /* end if */

    /* Wait until everyone gets here */
    fflush(NULL);
    my_barrier();

}                               /* end init_comp */


