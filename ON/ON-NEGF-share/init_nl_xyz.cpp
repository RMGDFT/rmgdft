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
#include "main.h"
#include "prototypes_on.h"



void init_nl_xyz (void)
{
    int ion;
    SPECIES *sp;
    ION *iptr;


    int L0_NLDIM;
    double hxgrid, hygrid, hzgrid;



    /* Loop over all the ions on this processor */

    for (ion = 0; ion < ct.num_ions; ion++)
    {
        /* Generate ion pointer */
        iptr = &Atoms[ion];

        /* Get species type */
        sp = &ct.sp[iptr->species];


        L0_NLDIM = sp->nldim;
        hxgrid = get_hxgrid() * get_xside();
        hygrid = get_hygrid() * get_yside();
        hzgrid = get_hzgrid() * get_zside();

        get_start (L0_NLDIM, Atoms[ion].crds[0], ct.xcstart, hxgrid,
                   &iptr->ixstart, &iptr->nlxcstart);
        get_start (L0_NLDIM, Atoms[ion].crds[1], ct.ycstart, hygrid,
                   &iptr->iystart, &iptr->nlycstart);
        get_start (L0_NLDIM, Atoms[ion].crds[2], ct.zcstart, hzgrid,
                   &iptr->izstart, &iptr->nlzcstart);

        iptr->nlxcstart /= get_xside();
        iptr->nlycstart /= get_yside();
        iptr->nlzcstart /= get_zside();
        iptr->ixend = iptr->ixstart + L0_NLDIM - 1;
        iptr->iyend = iptr->iystart + L0_NLDIM - 1;
        iptr->izend = iptr->izstart + L0_NLDIM - 1;


    }                           /* end for ion */

    if (pct.gridpe == 0)
    {

        printf (" init_nl_xyz.c  done\n");

    }                           /* end if */
    MPI_Barrier(pct.img_comm);
    fflush (NULL);


}                               /* end init_nl_xyz.c */
