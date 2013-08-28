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


void init_comp (REAL *vh)
{

    int ix, iy, iz, i;
    int pxoff;  
    double v_comp, v_sum, v_reference;    
    int x_locate;
    double delta_v1;
    double delta_v2;
    double delta_v3;
    double delta_v4;
    int lb = ct.vcomp_Lbegin * 2;  //lb,le,rb and re is now base on FINE grid
    int le = ct.vcomp_Lend * 2;
    int rb = ct.vcomp_Rbegin * 2;
    int re = ct.vcomp_Rend * 2;
    
    REAL t1;
    REAL *zvec;
    my_malloc_init( zvec, FNX_GRID, REAL );

    if (pct.gridpe == 0)
    {

        printf(" Begin init_comp ...\n");

    }                           /* end if */

    fflush(NULL);

    my_barrier();

    /* Get this processors offset */
    pxoff = pct.FPX_OFFSET;
    /* Zero out result vector */
    for (ix = 0; ix < FNX_GRID; ix++)
        zvec[ix] = ZERO;


    /* Loop over this processor */
    for (ix = 0; ix < FPX0_GRID; ix++)
    {
        t1 = 0.0;
        for (iy = 0; iy < FPY0_GRID; iy++)
            for (iz = 0; iz < FPZ0_GRID; iz++)

		    t1 += vh[ix * FPY0_GRID * FPZ0_GRID + iy * FPZ0_GRID + iz];


	zvec[ix + pxoff] = t1/ FNY_GRID / FNZ_GRID;

    }                           /* end for */


    /* Now sum over all processors */
    ix = FNX_GRID;
    global_sums(zvec, &ix, pct.grid_comm);


    printf(" vcomp_Lbegin = %d  and vcomp_Lend = %d  on FINE GRID \n", lb, le);   
    printf(" vcomp_Rbegin = %d  and vcomp_Rend = %d  on FINE GRID \n", rb, re);   
    printf(" vcomp_Lbegin = %d  and vcomp_Lend = %d  on coarse GRID \n", ct.vcomp_Lbegin, ct.vcomp_Lend);   
    printf(" vcomp_Rbegin = %d  and vcomp_Rend = %d  on coarse GRID \n", ct.vcomp_Rbegin, ct.vcomp_Rend);   
   
    delta_v1 = zvec[lb - 2] - zvec[lb];
    delta_v2 = zvec[lb - 2] - zvec[le];
    delta_v3 = zvec[re + 1] - zvec[rb];
    delta_v4 = zvec[re + 1] - zvec[re - 1];
    printf(" delta_v1  = %f \n", delta_v1);   
    printf(" delta_v2  = %f \n", delta_v2);   
    printf(" delta_v3  = %f \n", delta_v3);   
    printf(" delta_v4  = %f \n", delta_v4);   

    v_reference = 0;
    for (i = 0; i < lb; i++)
    {
	    v_reference = v_reference + zvec[i];

    }
    v_reference = v_reference / lb;

    for (ix = 0; ix < FPX0_GRID; ix++)
    {
        x_locate = ix + pxoff; // everything is based on FINE grid!

        if ( x_locate >= lb && x_locate <= le )
        
            v_comp = (delta_v2 - delta_v1) * (x_locate - lb) / (le - lb) + delta_v1; 
        
        else if (x_locate > rb && x_locate < re)

            v_comp = (delta_v4 - delta_v3) * (x_locate - rb)/ (re - rb) + delta_v3;

        else if (x_locate > le && x_locate <= rb)
	{       
                v_sum = zvec[x_locate];
		for (i = 1; i < 32; i++)
		{
                    v_sum = v_sum + zvec[x_locate-i] + zvec[x_locate+i];
                    
		}
                v_sum = v_sum / 63.0;

		v_comp = v_reference - v_sum;
	}
        else

            v_comp = 0;


        for (iy = 0; iy < FPY0_GRID; iy++)
        {

            for (iz = 0; iz < FPZ0_GRID; iz++)
            {

                vcomp[ix * FPY0_GRID * FPZ0_GRID + iy * FPZ0_GRID + iz] = v_comp;


            }                   /* end for */

        }                       /* end for */


	printf("x_locate = %5d      vcomp  = %10.7f \n ", x_locate, v_comp );

    }                           /* end for */



    if (pct.gridpe == 0)
    {

        printf(" init_comp done\n");

    }                           /* end if */

    /* Wait until everyone gets here */
    fflush(NULL);
    my_barrier();
    my_free(zvec);

}                               /* end init_comp */


