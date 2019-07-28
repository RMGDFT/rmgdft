#include "negf_prototypes.h"
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
#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include "prototypes_on.h"



void init_loc_xyz ()
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
        iptr = &ct.ions[ion];

        /* Get species type */
        sp = &ct.sp[iptr->species];


        L0_NLDIM = sp->ldim_coar;
        hxgrid = get_hxgrid() * get_xside();
        hygrid = get_hygrid() * get_yside();
        hzgrid = get_hzgrid() * get_zside();

        get_start (L0_NLDIM, ct.ions[ion].crds[0], ct.xcstart, hxgrid,
                   &iptr->ixstart_loc, &iptr->xcstart_loc);
        get_start (L0_NLDIM, ct.ions[ion].crds[1], ct.ycstart, hygrid,
                   &iptr->iystart_loc, &iptr->ycstart_loc);
        get_start (L0_NLDIM, ct.ions[ion].crds[2], ct.zcstart, hzgrid,
                   &iptr->izstart_loc, &iptr->zcstart_loc);

        iptr->xcstart_loc /= get_xside();
        iptr->ycstart_loc /= get_yside();
        iptr->zcstart_loc /= get_zside();
        iptr->ixend_loc = iptr->ixstart_loc + L0_NLDIM - 1;
        iptr->iyend_loc = iptr->iystart_loc + L0_NLDIM - 1;
        iptr->izend_loc = iptr->izstart_loc + L0_NLDIM - 1;


    }

    if (pct.gridpe == 0)
    {

        printf ("\n init_loc_xyz.c  done\n");

    } 
    MPI_Barrier(pct.img_comm);
    fflush (NULL);


}
