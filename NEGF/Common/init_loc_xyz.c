/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 

/*

                      init_loc_xyz.c 


*/




#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "md.h"



void init_loc_xyz ()
{
    int ion;
    SPECIES *sp;
    ION *iptr;


    int L0_NLDIM;
    REAL hxgrid, hygrid, hzgrid;



    /* Loop over all the ions on this processor */

    for (ion = 0; ion < ct.num_ions; ion++)
    {
        /* Generate ion pointer */
        iptr = &ct.ions[ion];

        /* Get species type */
        sp = &ct.sp[iptr->species];


        L0_NLDIM = sp->ldim_coar;
        hxgrid = ct.hxgrid * ct.xside;
        hygrid = ct.hygrid * ct.yside;
        hzgrid = ct.hzgrid * ct.zside;

        get_start (L0_NLDIM, ct.ions[ion].crds[0], ct.xcstart, hxgrid,
                   &iptr->ixstart_loc, &iptr->xcstart_loc);
        get_start (L0_NLDIM, ct.ions[ion].crds[1], ct.ycstart, hygrid,
                   &iptr->iystart_loc, &iptr->ycstart_loc);
        get_start (L0_NLDIM, ct.ions[ion].crds[2], ct.zcstart, hzgrid,
                   &iptr->izstart_loc, &iptr->zcstart_loc);

        iptr->xcstart_loc /= ct.xside;
        iptr->ycstart_loc /= ct.yside;
        iptr->zcstart_loc /= ct.zside;
        iptr->ixend_loc = iptr->ixstart_loc + L0_NLDIM - 1;
        iptr->iyend_loc = iptr->iystart_loc + L0_NLDIM - 1;
        iptr->izend_loc = iptr->izstart_loc + L0_NLDIM - 1;


    }

    if (pct.thispe == 0)
    {

        printf ("\n init_loc_xyz.c  done\n");

    } 
    my_barrier ();
    fflush (NULL);


}
