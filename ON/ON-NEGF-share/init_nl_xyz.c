/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 

/*

                      init_nl_xyz.c 




*/




#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "md.h"



void init_nl_xyz (void)
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


        L0_NLDIM = sp->nldim;
        hxgrid = ct.hxgrid * ct.xside;
        hygrid = ct.hygrid * ct.yside;
        hzgrid = ct.hzgrid * ct.zside;

        get_start (L0_NLDIM, ct.ions[ion].crds[0], ct.xcstart, hxgrid,
                   &iptr->ixstart, &iptr->nlxcstart);
        get_start (L0_NLDIM, ct.ions[ion].crds[1], ct.ycstart, hygrid,
                   &iptr->iystart, &iptr->nlycstart);
        get_start (L0_NLDIM, ct.ions[ion].crds[2], ct.zcstart, hzgrid,
                   &iptr->izstart, &iptr->nlzcstart);

        iptr->nlxcstart /= ct.xside;
        iptr->nlycstart /= ct.yside;
        iptr->nlzcstart /= ct.zside;
        iptr->ixend = iptr->ixstart + L0_NLDIM - 1;
        iptr->iyend = iptr->iystart + L0_NLDIM - 1;
        iptr->izend = iptr->izstart + L0_NLDIM - 1;


    }                           /* end for ion */

    if (pct.gridpe == 0)
    {

        printf (" init_nl_xyz.c  done\n");

    }                           /* end if */
    my_barrier ();
    fflush (NULL);


}                               /* end init_nl_xyz.c */
