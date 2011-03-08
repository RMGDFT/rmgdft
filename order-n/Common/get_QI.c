/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "md.h"


void get_QI(void)
{
    int isp, idx, idx1, i, j, k, l, ion, size;
    int ix, iy, iz, nty, ih, num;
    int lpx[9][9], lpl[9][9][9], nhtol[MAX_SPECIES][18];
    int nhtom[MAX_SPECIES][18], indv[MAX_SPECIES][18], nh[MAX_SPECIES];
    int ilow, jlow, klow, ihi, jhi, khi, map, icount;
    int *Aix, *Aiy, *Aiz;
    int icut, itmp, icenter, lmax, alloc;
    int *pvec, *ivec, *dvec;
    REAL x[3], cx[3], r, invdr, ap[25][9][9], ylm[25];
    REAL xc, yc, zc;
    ION *iptr;
    REAL *qnmlig_tpr, *QI_tpr;
    SPECIES *sp;

    my_malloc_init( Aix, FNX_GRID, int );
    my_malloc_init( Aiy, FNY_GRID, int );
    my_malloc_init( Aiz, FNZ_GRID, int );

    for (k = 0; k < 9; k++)
    {
        for (j = 0; j < 9; j++)
        {
            for (i = 0; i < 25; i++)
            {
                ap[i][j][k] = 0.0;
                if (i < 9)
                    lpl[i][j][k] = 0;
            }
            lpx[k][j] = 0;
        }
    }

    for (i = 0; i < MAX_SPECIES; i++)
    {
        for (j = 0; j < 18; j++)
        {
            nhtol[i][j] = 0;
            nhtom[i][j] = 0;
            indv[i][j] = 0;
        }
        nh[i] = 0;
    }

    lmax = 0;

    for (isp = 0; isp < ct.num_species; isp++)
    {

        sp = &ct.sp[isp];
        if (sp->llbeta[sp->nbeta - 1] > lmax)
            lmax = sp->llbeta[sp->nbeta - 1];

        ih = 0;
        for (i = 0; i < sp->nbeta; i++)
        {
            l = sp->llbeta[i];
            for (j = 0; j < 2 * l + 1; j++)
            {
                nhtol[isp][ih] = l;
                nhtom[isp][ih] = j;
                indv[isp][ih] = i;
                ih = ih + 1;
            }
        }
        nh[isp] = ih;

    }                           /*end for isp */

    aainit(lmax + 1, 2 * lmax + 1, 2 * lmax + 1, 4 * lmax + 1,
           (lmax + 1) * (lmax + 1), ap, lpx, lpl);

    alloc = ct.max_Qpoints;
    my_malloc_init( pvec, 2 * alloc, int );
    dvec = pvec + alloc;

    for (ion = 0; ion < ct.num_ions; ion++)
    {

        iptr = &ct.ions[ion];
        nty = iptr->species;
        sp = &ct.sp[nty];

        icenter = sp->qdim / 2;
        icut = (icenter + 1) * (icenter + 1);

        map =
            get_index(iptr, Aix, Aiy, Aiz, &ilow, &ihi, &jlow, &jhi, &klow,
                      &khi, sp->qdim, FPX0_GRID, FPY0_GRID, FPZ0_GRID,
                      ct.vh_nxgrid, ct.vh_nygrid, ct.vh_nzgrid,
                      &iptr->Qxcstart, &iptr->Qycstart, &iptr->Qzcstart);

        /* If there is a mapping for this ion then we have to generate */
        /* the projector.                                              */

        for (idx = 0; idx < alloc; idx++)
        {
            pvec[idx] = 0;
            dvec[idx] = 0;
        }

        icount = 0;
        if (map)
        {

            /* Generate index arrays */
            idx = 0;
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
                            itmp = (ix - icenter) * (ix - icenter) +
                                (iy - icenter) * (iy - icenter) + (iz - icenter) * (iz - icenter);

                            if (icut >= itmp)
                            {

                                pvec[icount] =
                                    FPY0_GRID * FPZ0_GRID * (Aix[ix] %
                                                             FPX0_GRID) +
                                    FPZ0_GRID * (Aiy[iy] % FPY0_GRID) + (Aiz[iz] % FPZ0_GRID);

                                dvec[idx] = TRUE;

                                icount++;

                            }   /* end if icut */

                        }       /*end if */
                        idx++;

                    }           /* end for iz */

                }               /* end for iy */

            }                   /* end for ix */

            if (pct.Qdvec[ion] != NULL)
                my_free(pct.Qdvec[ion]);

            if (icount)
            {
                my_calloc( pct.Qdvec[ion], idx, int );
                if (NULL == pct.Qdvec[ion])
                    error_handler("Unable to allocate Qdevc array space");

                ivec = pct.Qdvec[ion];
                for (idx1 = 0; idx1 < idx; idx1++)
                {

                    ivec[idx1] = (int) dvec[idx1];

                }
            }
            else
                pct.Qdvec[ion] = NULL;


            if (pct.Qindex[ion] != NULL)
                my_free(pct.Qindex[ion]);

            if (icount)
            {
                my_calloc( pct.Qindex[ion], icount + 128, int );
                if (NULL == pct.Qindex[ion])
                    error_handler("Unable to allocate Qindex array space");

                ivec = pct.Qindex[ion];
                for (idx1 = 0; idx1 < icount; idx1++)
                {

                    ivec[idx1] = (int) pvec[idx1];

                }
            }
            else
                pct.Qindex[ion] = NULL;

            if (pct.augfunc[ion] != NULL)
                my_free(pct.augfunc[ion]);

            if (icount)
            {

                size = nh[nty] * (nh[nty] + 1) / 2;
                my_calloc( pct.augfunc[ion], size * icount + 128, REAL );
                if (NULL == pct.augfunc[ion])
                    error_handler("Unable to allocate augfunc array space");


            }
            else
                pct.augfunc[ion] = NULL;

        }                       /*end for if map */

        pct.Qidxptrlen[ion] = icount;

        if (map)
        {
            invdr = ONE / sp->drqlig;
            QI_tpr = pct.augfunc[ion];
            qnmlig_tpr = sp->qnmlig;

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

                            r = metric(x);
                            to_cartesian(x, cx);
                            ylmr2(cx, ylm);

                            num = 0;
                            for (i = 0; i < nh[nty]; i++)
                            {

                                for (j = i; j < nh[nty]; j++)
                                {

                                    idx1 = num * pct.Qidxptrlen[ion] + icount;
                                    QI_tpr[idx1] =
                                        qval(i, j, r, invdr, qnmlig_tpr,
                                             &nhtol[nty][0], &nhtom[nty][0],
                                             &indv[nty][0], ylm, ap, lpx, lpl, sp);


                                    num++;
                                }       /*end for j */

                            }   /* end for i */

                            icount++;

                        }       /* end if */

                        idx++;

                        zc += ct.hzzgrid;

                    }           /* end for iz */
                    yc += ct.hyygrid;

                }               /* end for iy */
                xc += ct.hxxgrid;

            }                   /* end for ix */

        }                       /*end if map */

    }                           /*end for ion */

    if (pct.thispe == 0)
    {

        printf(" get_QI.c  done\n");

    }                           /* end if */
    my_barrier();
    fflush(NULL);

    my_free(pvec);
    my_free(Aix);
    my_free(Aiy);
    my_free(Aiz);

}
