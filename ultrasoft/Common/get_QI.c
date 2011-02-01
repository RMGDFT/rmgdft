/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"

/* Sets Qnm function (part of ultrasfot pseudpotential*/


void get_QI (void)
{
    int is, idx, idx1, i, j, k, l, ion, size;
    int ix, iy, iz, nty, ih, num;
    int lpx[9][9], lpl[9][9][9];
    int nhtol[MAX_SPECIES][18];
    int nhtom[MAX_SPECIES][18];
    int indv[MAX_SPECIES][18];
    int nh[MAX_SPECIES];
    int ilow, jlow, klow, ihi, jhi, khi, map, icount;
    int Aix[FNX_GRID], Aiy[FNY_GRID], Aiz[FNZ_GRID];
    int icut, itmp, icenter, lmax, alloc;
    int *pvec, *ivec, *dvec;
    REAL x[3], cx[3], r, invdr, ap[25][9][9], ylm[25];
    REAL xc, yc, zc;
    ION *iptr;
    REAL *qnmlig_tpr, *QI_tpr;
    SPECIES *sp;



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

    for (is = 0; is < ct.num_species; is++)
    {
        for (j = 0; j < 18; j++)
        {
            nhtol[is][j] = 0;
            nhtom[is][j] = 0;
            indv[is][j] = 0;
        }
        nh[is] = 0;
    }

    lmax = 0;

    for (is = 0; is < ct.num_species; is++)
    {

        sp = &ct.sp[is];
        if (sp->llbeta[sp->nbeta - 1] > lmax)
            lmax = sp->llbeta[sp->nbeta - 1];

        ih = 0;
        for (i = 0; i < sp->nbeta; i++)
        {
            l = sp->llbeta[i];
            for (j = 0; j < 2 * l + 1; j++)
            {
                nhtol[is][ih] = l;
                nhtom[is][ih] = j;
                indv[is][ih] = i;
                ih = ih + 1;
            }
        }
        nh[is] = ih;

    }

    aainit (lmax + 1, 2 * lmax + 1, 2 * lmax + 1, 4 * lmax + 1, (lmax + 1) * (lmax + 1), ap, lpx,
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



        iptr = &ct.ions[ion];
        nty = iptr->species;
        sp = &ct.sp[nty];

        icenter = sp->qdim / 2;
        icut = (icenter + 1) * (icenter + 1);

        map = get_index (iptr, Aix, Aiy, Aiz, &ilow, &ihi, &jlow, &jhi, &klow, &khi,
                         sp->qdim, FPX0_GRID, FPY0_GRID, FPZ0_GRID,
                         ct.psi_fnxgrid, ct.psi_fnygrid, ct.psi_fnzgrid,
                         &iptr->Qxcstart, &iptr->Qycstart, &iptr->Qzcstart);

        /*
           if(pct.thispe==0) printf("qdim=%d  Qxcstart=%f  Qycstart=%f  Qzcstart=%f\n",sp->qdim,iptr->Qxcstart,iptr->Qycstart,iptr->Qzcstart);
         */
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
                                pvec[icount] = (FPY0_GRID * FPZ0_GRID * (Aix[ix] % FPX0_GRID) +
                                                FPZ0_GRID * (Aiy[iy] % FPY0_GRID) +
                                                (Aiz[iz] % FPZ0_GRID));

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


                size = nh[nty] * (nh[nty] + 1) / 2;
                my_calloc (pct.augfunc[ion], size * icount + 128, REAL);

            }

            pct.Qidxptrlen[ion] = icount;


            invdr = 1.0 / sp->drqlig;
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

                            r = metric (x);
                            to_cartesian (x, cx);
                            ylmr2 (cx, ylm);

                            num = 0;
                            for (i = 0; i < nh[nty]; i++)
                            {

                                for (j = i; j < nh[nty]; j++)
                                {

                                    idx1 = num * pct.Qidxptrlen[ion] + icount;
                                    QI_tpr[idx1] = qval (i, j, r, invdr, qnmlig_tpr, &nhtol[nty][0],
                                                         &nhtom[nty][0], &indv[nty][0], ylm, ap,
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
