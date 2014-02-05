/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "grid.h"
#include "main.h"
#include "common_prototypes.h"

/* Sets Qnm function (part of ultrasfot pseudpotential*/


void get_QI (void)
{
    int is, idx, idx1, i, j, k, l, ion, size;
    int ix, iy, iz, species, ih, num;
    int lpx[9][9], lpl[9][9][9];
    int nh;
    int ilow, jlow, klow, ihi, jhi, khi, map, icount;
    int Aix[FNX_GRID], Aiy[FNY_GRID], Aiz[FNZ_GRID];
    int icut, itmp, icenter, alloc;
    int *pvec, *ivec, *dvec;
    rmg_double_t x[3], cx[3], r, invdr, ap[25][9][9], ylm[25];
    rmg_double_t xc, yc, zc;
    ION *iptr;
    rmg_double_t *qnmlig, *QI_tpr;
    SPECIES *sp;

    // If norm conserving pp just return
    if(ct.norm_conserving_pp) return;

    /* Initialize some coefficients (not sure what exactly)*/
    aainit (ct.max_l + 1, 2 * ct.max_l + 1, 2 * ct.max_l + 1, 4 * ct.max_l + 1, (ct.max_l + 1) * (ct.max_l + 1), ap, lpx,
            lpl);

    alloc = ct.max_Qpoints;
    my_malloc (pvec, 2 * alloc, int);
    dvec = pvec + alloc;

    for (ion = 0; ion < ct.num_ions; ion++)
    {

        /*Release memory first */
        if (pct.Qdvec[ion])
            my_free (pct.Qdvec[ion]);
        if (pct.Qindex[ion])
            my_free (pct.Qindex[ion]);
        if (pct.augfunc[ion])
            my_free (pct.augfunc[ion]);

        /*Let those empty pointers point to NULL */
        pct.Qdvec[ion] = NULL;
        pct.Qindex[ion] = NULL;
        pct.augfunc[ion] = NULL;

        /*Initialize this */
        pct.Qidxptrlen[ion] = 0;
        
        for (idx = 0; idx < alloc; idx++)
        {
            pvec[idx] = 0;
            dvec[idx] = 0;
        }

        icount = 0;


        iptr = &ct.ions[ion];
        species = iptr->species;
        sp = &ct.sp[species];

        nh = sp->nh;

        icenter = sp->qdim / 2;
        icut = (icenter + 1) * (icenter + 1);

        map = get_index (pct.gridpe, iptr, Aix, Aiy, Aiz, &ilow, &ihi, &jlow, &jhi, &klow, &khi,
                         sp->qdim, get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(),
                         ct.psi_fnxgrid, ct.psi_fnygrid, ct.psi_fnzgrid,
                         &iptr->Qxcstart, &iptr->Qycstart, &iptr->Qzcstart);

        if (map)
        {

            /* Generate index arrays */
            icount = idx = 0;
            for (ix = 0; ix < sp->qdim; ix++)
            {
                for (iy = 0; iy < sp->qdim; iy++)
                {
                    for (iz = 0; iz < sp->qdim; iz++)
                    {
                        dvec[idx] = FALSE;
                        if ((((Aix[ix] >= ilow) && (Aix[ix] <= ihi)) &&
                             ((Aiy[iy] >= jlow) && (Aiy[iy] <= jhi)) &&
                             ((Aiz[iz] >= klow) && (Aiz[iz] <= khi))))
                        {

                            /* Cut it off if required */
                            itmp = ((ix - icenter) * (ix - icenter) +
                                    (iy - icenter) * (iy - icenter) +
                                    (iz - icenter) * (iz - icenter));

                            if (icut >= itmp)
                            {
                                pvec[icount] =
                                                get_FPY0_GRID() * get_FPZ0_GRID() * ((Aix[ix]-get_FPX_OFFSET()) % get_FPX0_GRID()) +
                                                get_FPZ0_GRID() * ((Aiy[iy]-get_FPY_OFFSET()) % get_FPY0_GRID()) +
                                                ((Aiz[iz]-get_FPZ_OFFSET()) % get_FPZ0_GRID());


                                dvec[idx] = TRUE;

                                icount++;

                            }   /* end if icut */

                        }
                        idx++;
                    }
                }
            }


            if (icount)
            {


                my_calloc (pct.Qdvec[ion], idx, int);

                ivec = pct.Qdvec[ion];
                for (idx1 = 0; idx1 < idx; idx1++)
                    ivec[idx1] = (int) dvec[idx1];



                my_calloc (pct.Qindex[ion], icount + 128, int);

                ivec = pct.Qindex[ion];
                for (idx1 = 0; idx1 < icount; idx1++)
                    ivec[idx1] = (int) pvec[idx1];


                size = nh * (nh + 1) / 2;
                my_calloc (pct.augfunc[ion], size * icount + 128, rmg_double_t);

            }

            pct.Qidxptrlen[ion] = icount;


            invdr = 1.0 / sp->drqlig;
            QI_tpr = pct.augfunc[ion];
            qnmlig = sp->qnmlig;

            idx = 0;
            icount = 0;
            /* Generate index arrays */
            xc = iptr->Qxcstart;
            for (ix = 0; ix < sp->qdim; ix++)
            {

                yc = iptr->Qycstart;
                for (iy = 0; iy < sp->qdim; iy++)
                {

                    zc = iptr->Qzcstart;
                    for (iz = 0; iz < sp->qdim; iz++)
                    {

                        if (dvec[idx])
                        {
                            x[0] = xc - iptr->xtal[0];
                            x[1] = yc - iptr->xtal[1];
                            x[2] = zc - iptr->xtal[2];

                            r = metric (x);
                            to_cartesian (x, cx);
                            ylmr2 (cx, ylm);

                            num = 0;
                            for (i = 0; i < nh; i++)
                            {

                                for (j = i; j < nh; j++)
                                {

                                    idx1 = num * pct.Qidxptrlen[ion] + icount;
                                    QI_tpr[idx1] = qval (i, j, r, invdr, qnmlig,sp->nhtol,
                                                         sp->nhtom, sp->indv, ylm, ap,
                                                         lpx, lpl, sp);

                                    num++;
                                }
                            }
                            icount++;
                        }

                        idx++;

                        zc += ct.hzzgrid;

                    }           /* end for iz */
                    yc += ct.hyygrid;

                }               /* end for iy */
                xc += ct.hxxgrid;

            }                   /* end for ix */

        }                       /*end if map */

    }                           /*end for ion */

    my_free (pvec);

}
