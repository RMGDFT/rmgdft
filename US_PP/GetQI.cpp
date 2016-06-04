/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "grid.h"

#include "make_conf.h"
#include "const.h"
#include "params.h"
#include "grid.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
//#include <complex>
//#include "Kpoint.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "AtomicInterpolate.h"
#include "transition.h"


extern double Atomic_inv_a;
extern double Atomic_inv_b;

/* Sets Qnm function (part of ultrasfot pseudpotential*/

static inline double qval_inline (int lmax, int ih, int jh, int ic, double d0, double d1, double dm,  
        double * ptpr, int *nh_l2m,
        int *indv, double * ylm, double *ap, int *lpx, int *lpl, SPECIES * sp);

//void InitClebschGordan(int, int *[9][9], int*[9], int*[9][9]);
void GetQI (void)
{
    int idx, idx1, i, j,ion, size;
    int ix, iy, iz, species, num;
    int *lpx, *lpl;
    double *ap, *ylm;
    int nh;
    int ilow, jlow, klow, ihi, jhi, khi, map, icount;
    int *Aix, *Aiy, *Aiz;
    int icut, itmp, icenter, alloc;
    int *pvec, *ivec, *dvec;
    double x[3], cx[3], r;
    double xc, yc, zc;
    double hxxgrid, hyygrid, hzzgrid;
    ION *iptr;
    double *qnmlig, *QI_tpr;
    SPECIES *sp;

    hxxgrid = get_hxxgrid();
    hyygrid = get_hyygrid();
    hzzgrid = get_hzzgrid();

    // If norm conserving pp just return
    if(ct.norm_conserving_pp) return;

    int num_lm = (ct.max_l + 1) * (ct.max_l+1);
    int num_LM2 = (2*ct.max_l + 1) * (2*ct.max_l+1);

    lpx = new int[num_lm * num_lm];
    lpl = new int[num_lm * num_lm  * num_LM2];
    ap = new double[num_LM2 * num_lm * num_lm];
    ylm = new double[num_LM2];
    //ylm = new double[25];

    InitClebschGordan(ct.max_l, ap, lpx, lpl);


    alloc = ct.max_Qpoints;
    pvec = new int[2 * alloc];
    dvec = pvec + alloc;

    Aix = new int[get_FNX_GRID()];
    Aiy = new int[get_FNY_GRID()];
    Aiz = new int[get_FNZ_GRID()];

    for (ion = 0; ion < ct.num_ions; ion++)
    {

        /*Release memory first */
        if (pct.Qdvec[ion])
            delete [](pct.Qdvec[ion]);
        if (pct.Qindex[ion])
            delete [](pct.Qindex[ion]);
        if (pct.augfunc[ion])
            delete [](pct.augfunc[ion]);

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
                get_FNX_GRID(), get_FNY_GRID(), get_FNZ_GRID(),
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


            if (!icount) continue;
            pct.Qidxptrlen[ion] = icount;


            pct.Qdvec[ion] = new int[idx];

            ivec = pct.Qdvec[ion];
            for (idx1 = 0; idx1 < idx; idx1++)
                ivec[idx1] = (int) dvec[idx1];


            pct.Qindex[ion] = new int[icount + 128];

            ivec = pct.Qindex[ion];
            for (idx1 = 0; idx1 < icount; idx1++)
                ivec[idx1] = (int) pvec[idx1];


            size = nh * (nh + 1) / 2;
            pct.augfunc[ion] = new double[ size * icount + 128];


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
                            for(int l = 0; l <= 2*ct.max_l; l++)
                            for(int m = 0; m < 2*l + 1; m++)
                                ylm[l*l+m] = Ylm(l, m, cx);
                    
                            //ylmr2 (cx, ylm);

                            if((r < LOGGRID_START)) {
                                r = LOGGRID_START;
                            }

                            double d0, d1, dm, ic;
                            d0 = (log (r*Atomic_inv_a) * Atomic_inv_b);
                            ic = (int)d0;
                            ic = (ic > 0) ? ic : 1;

                            /* cubic interpolation using forward differences */
                            d0 -= (double) (ic);
                            d1 = (d0 - 1.0) * 0.5;
                            dm = (d0 - 2.0) / 3.0;


                            num = 0;
                            for (i = 0; i < nh; i++)
                            {

                                for (j = i; j < nh; j++)
                                {

                                    idx1 = num * pct.Qidxptrlen[ion] + icount;
                                    QI_tpr[idx1] = qval_inline (ct.max_l, i, j, ic, d0, d1, dm, qnmlig,sp->nh_l2m,
                                            sp->indv, ylm, ap, lpx, lpl, sp);

                                    num++;
                                }
                            }
                            icount++;
                        }

                        idx++;

                        zc += hzzgrid;

                    }           /* end for iz */
                    yc += hyygrid;

                }               /* end for iy */
                xc += hxxgrid;

            }                   /* end for ix */

        }                       /*end if map */

        // No longer needed so release it
        if (pct.Qdvec[ion]) delete [] pct.Qdvec[ion];
        pct.Qdvec[ion] = NULL;

    }                           /*end for ion */


    delete [](Aiz);
    delete [](Aiy);
    delete [](Aix);
    delete [](pvec);
    delete []lpx;
    delete []lpl;
    delete []ap;
    delete []ylm;


}

static inline double qval_inline (int lmax, int ih, int jh, int ic, double d0, double d1, double dm,  
        double * ptpr, int *nh_l2m,
        int *indv, double * ylm, double *ap, int *lpx, int *lpl, SPECIES * sp)
{
    int ivl, jvl;
    int nb, mb, nmb, lm, lp, l;
    double qrad, sum, *ptpr1;



    nb = indv[ih];
    mb = indv[jh];
    if (nb < mb)
        nmb = mb * (mb + 1) / 2 + nb;
    else
        nmb = nb * (nb + 1) / 2 + mb;

    int num_lm = (lmax +1 ) * (lmax +1);
    int num_LM2 = (2*lmax+1) *(2*lmax +1);

    ivl = nh_l2m[ih];
    jvl = nh_l2m[jh];

    sum = 0;

    for (lm = 0; lm < lpx[ivl*num_lm + jvl]; lm++)
    {
        lp = lpl[(ivl*num_lm + jvl) * num_LM2 + lm];
        if(lp == 0)
            l = 0;
        else if (lp < 4)
            l = 1;
        else if (lp < 9)
            l = 2;
        else 
            l = (int)sqrt(lp + 0.1);

        ptpr1 = ptpr + (nmb * sp->nlc + l) * MAX_LOGGRID;
        qrad = AtomicInterpolateInline_1 (ptpr1, ic, d0, d1, dm);
        sum += qrad * ap[lp*num_lm * num_lm + ivl * num_lm + jvl] * ylm[lp];
    }
    return (sum);
}
