/************************** SVN Revision Information **************************
 **    $Id: get_start.c 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 
#include "md.h"
#include <math.h>



void get_start (int cdim, double crds, double cstart, double hgrid, int *istart, double *nlcstart)
{
    double t1, t2;
    int ic;

    /* t1 = nb of nodes between the boundary and crds */
    t1 = (crds - cstart) / hgrid;

    /* get the integral part of t1 in ic */
    t1 = modf (t1, &t2);
    ic = (int) t2;
    if (t1 > 0.5)
        ic++;

    *istart = ic - cdim / 2;
    *nlcstart = (cstart + hgrid * *istart);
}

/* Generate index array "pvec" and 
   an array "dvec" with value 1 if the node is in the region */
void get_index_array (int *pvec, int *dvec, int dim,
                      int *Aix, int *Aiy, int *Aiz, int *icount,
                      int index_low[3], int index_high[3])
{
    int idx, iz, iy, ix;

    *icount = idx = 0;

    for (ix = 0; ix < dim; ix++)
        for (iy = 0; iy < dim; iy++)
            for (iz = 0; iz < dim; iz++)
            {

                dvec[idx] = FALSE;
                if ((((Aix[ix] >= index_low[0]) && (Aix[ix] <= index_high[0])) &&
                     ((Aiy[iy] >= index_low[1]) && (Aiy[iy] <= index_high[1])) &&
                     ((Aiz[iz] >= index_low[2]) && (Aiz[iz] <= index_high[2]))))
                {


                    pvec[*icount] =
                        PY0_GRID * PZ0_GRID * (Aix[ix] % PX0_GRID) +
                        PZ0_GRID * (Aiy[iy] % PY0_GRID) + (Aiz[iz] % PZ0_GRID);

                    dvec[idx] = TRUE;

                    (*icount)++;


                }               /* end if */
                idx++;


            }                   /* end for */
}

/* Generate range of indices over which the projectors */
/* will be mapped onto the global grid.                */
void get_Ai (int cdim, int *Ai, int ngrid, int istart)
{
    int ix, idx;

    idx = istart;
    for (ix = 0; ix < cdim; ix++)
    {

        if (idx < 0)
        {

            Ai[ix] = ngrid + idx;

        }
        else if (idx > (ngrid - 1))
        {

            Ai[ix] = idx - ngrid;

        }
        else
        {

            Ai[ix] = idx;

        }                       /* end if */

        idx++;

    }                           /* end for */

}
