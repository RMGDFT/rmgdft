/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*

                        partial_vloc.c


    Sets up the ket part of the non-local operators.



*/




#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "method.h"
#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include "fftw3.h"


void partial_vloc ()
{
    int ion, idx, overlap, st1;
    int tot_prj, ion1, index;
    int PROJECTOR_SPACE, ix, iy, iz;
    double r, ax[3], bx[3], xc, yc, zc, t1, invdr;
    double *vloc_x, *vloc_y, *vloc_z;
    ION *iptr;
    SPECIES *sp;


    my_barrier ();

    /*  get total number of vnuc on this processor */
    /*  pct.n_ion_center: number of ions whose nl projector overlap
     *  with the states on this processor */

    pct.n_ion_center_loc = 0;
    tot_prj = 0;
    for (ion = 0; ion < ct.num_ions; ion++)
    {
        overlap = 0;
        for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
        {
            index = (st1 - ct.state_begin) * ct.num_ions + ion;
            if (ion_orbit_overlap_region_loc[index].flag == 1)
                overlap = 1;
        }
        if (overlap == 1)
        {
            pct.ionidx_loc[pct.n_ion_center_loc] = ion;
            pct.n_ion_center_loc += 1;
            tot_prj += 1;
        }
    }

    PROJECTOR_SPACE = ct.max_lpoints * tot_prj;

    /*allocate memorry for weight factor of vnuc_x */
    if (vnuc_x != NULL)
        my_free (vnuc_x);
    my_malloc_init( vnuc_x, PROJECTOR_SPACE, double );

    /*allocate memorry for weight factor of vnuc_y */
    if (vnuc_y != NULL)
        my_free (vnuc_y);
    my_malloc_init( vnuc_y, PROJECTOR_SPACE, double);

    /*allocate memorry for weight factor of vnuc_z */
    if (vnuc_z != NULL)
        my_free (vnuc_z);
    my_malloc_init( vnuc_z, PROJECTOR_SPACE, double);

    vloc_x = vnuc_x;
    vloc_y = vnuc_y;
    vloc_z = vnuc_z;
    for (ion1 = 0; ion1 < pct.n_ion_center_loc; ion1++)
    {
        ion = pct.ionidx_loc[ion1];

        /* Generate ion pointer */
        iptr = &ct.ions[ion];

        /* Get species type */
        sp = &ct.sp[iptr->species];

        invdr = 1. / sp->drlig;
 
        idx = 0;
        xc = iptr->xcstart_loc;
        for (ix = 0; ix < sp->ldim_coar; ix++)
        {
 
            yc = iptr->ycstart_loc;
            for (iy = 0; iy < sp->ldim_coar; iy++)
            {
 
                zc = iptr->zcstart_loc;
                for (iz = 0; iz < sp->ldim_coar; iz++)
                {
 
                    ax[0] = xc - iptr->xtal[0];
                    ax[1] = yc - iptr->xtal[1];
                    ax[2] = zc - iptr->xtal[2];
 
                    r = metric(ax);
                    t1 = linint (sp->drlocalig, r, invdr);
                    to_cartesian(ax, bx);
                    r += 1.0e-10;

                    vloc_x[idx] = -t1 * bx[0] / r;
                    vloc_y[idx] = -t1 * bx[1] / r;
                    vloc_z[idx] = -t1 * bx[2] / r;

 
                    idx++;
                    zc += get_hzgrid();
                }/* end for ix*/
 
                yc += get_hygrid();
            }/* end for iy*/
 
            xc += get_hxgrid();
        }/* end for ix*/

        vloc_x += ct.max_lpoints;
        vloc_y += ct.max_lpoints;
        vloc_z += ct.max_lpoints;
    } /* end for ion1 */


    if (pct.gridpe == 0)
    {

        printf (" partial_vloc.c  done\n");

    }
    my_barrier ();
    fflush (NULL);

}
