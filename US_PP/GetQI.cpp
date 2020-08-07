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



#include "const.h"
#include "params.h"

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

/* Sets Qnm function (part of ultrasoft pseudpotential*/

static inline double qval_inline (int lmax, int ih, int jh, int ic, double d0, double d1, double dm,  
        double * ptpr, int *nh_l2m,
        int *indv, double * ylm, double *ap, int *lpx, int *lpl, SPECIES * sp);

static inline void map_qval_components(int ih, int jh, int *lpx, int *lpl, SPECIES * sp, std::multimap<size_t, qongrid> &aug_desc);

//void InitClebschGordan(int, int *[9][9], int*[9], int*[9][9]);
void GetQI (void)
{
    int idx, idx1, ion;
    int ix, iy, iz, species;
    int *lpx, *lpl;
    double *ylm;
    int nh;
    int ilow, jlow, klow, ihi, jhi, khi, map, icount;
    int *Aix, *Aiy, *Aiz;
    int icut, itmp, icenter, alloc;
    int *pvec, *dvec;
    double x[3], cx[3], r;
    double xc, yc, zc;
    double hxxgrid, hyygrid, hzzgrid;
    ION *iptr;
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
    ylm = new double[num_LM2];

    ct.cg_coeff.resize(boost::extents[num_lm][num_lm][num_LM2]);
    InitClebschGordan(ct.max_l, ct.cg_coeff.data(), lpx, lpl);



    alloc = ct.max_Qpoints;
    pvec = new int[2 * alloc];
    dvec = pvec + alloc;

    Aix = new int[get_FNX_GRID()];
    Aiy = new int[get_FNY_GRID()];
    Aiz = new int[get_FNZ_GRID()];

    ct.q_alloc[0] = 0;

    for (ion = 0; ion < ct.num_ions; ion++)
    {

        /*Release memory first */
        Atoms[ion].Qindex.clear();
        Atoms[ion].augfunc_desc.clear();
        Atoms[ion].stress_cx[0] = {};
        Atoms[ion].stress_cx[1] = {};
        Atoms[ion].stress_cx[2] = {};
        Atoms[ion].stress_cx[0].clear();
        Atoms[ion].stress_cx[1].clear();
        Atoms[ion].stress_cx[2].clear();

        Atoms[ion].grid_ylm.clear();
        Atoms[ion].grid_qr.clear();

        for (idx = 0; idx < alloc; idx++)
        {
            pvec[idx] = 0;
            dvec[idx] = 0;
        }

        icount = 0;


        iptr = &Atoms[ion];
        species = iptr->species;
        sp = &Species[species];

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

         //   Atoms[ion].Qindex.resize(icount);
            //for (idx1 = 0; idx1 < icount; idx1++) Atoms[ion].Qindex[idx1] = (int) pvec[idx1];
            for (idx1 = 0; idx1 < icount; idx1++) Atoms[ion].Qindex.push_back( pvec[idx1] );


            Atoms[ion].grid_ylm.resize((2*sp->max_l+1)*(2*sp->max_l+1));


            for(int l = 0; l <= 2*sp->max_l; l++)
            {
                for(int m = 0; m < 2*l + 1; m++)
                {
                    Atoms[ion].grid_ylm[lm_key(l, m)] = {};
                    Atoms[ion].grid_ylm[lm_key(l, m)].resize(icount);
                }
            }

            for (auto& qfunc: sp->qnmlig)
            {
                int nb = qnm_ival(qfunc.first);
                int mb = qnm_jval(qfunc.first);
                int l = qnm_lval(qfunc.first);
                Atoms[ion].grid_qr[qnm_key(nb, mb, l)] = {};
                Atoms[ion].grid_qr[qnm_key(nb, mb, l)].resize(icount);
            }

            for (int ih = 0; ih < nh; ih++)
            {
                for (int jh = ih; jh < nh; jh++)
                {
                    map_qval_components(ih, jh, lpx, lpl, sp, Atoms[ion].augfunc_desc);
                }
            }

            if(ct.stress)
            {
                Atoms[ion].stress_cx[0].resize(icount);
                Atoms[ion].stress_cx[1].resize(icount);
                Atoms[ion].stress_cx[2].resize(icount);
            }

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

                            if(ct.stress)
                            {
                                Atoms[ion].stress_cx[0][idx] = cx[0];
                                Atoms[ion].stress_cx[1][idx] = cx[1];
                                Atoms[ion].stress_cx[2][idx] = cx[2];
                            }

                            for(int l = 0; l <= 2*ct.max_l; l++)
                            {
                                for(int m = 0; m < 2*l + 1; m++)
                                {
                                    ylm[l*l+m] = Ylm(l, m, cx);
                                }
                            }

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


                            for(int l = 0; l <= 2*sp->max_l; l++)
                            {
                                for(int m = 0; m < 2*l + 1; m++)
                                {
                                    Atoms[ion].grid_ylm[lm_key(l,m)][icount] = Ylm(l, m, cx);
                                }
                            }

                            // Generate the q-functions on the 3D grid
                            for (auto& qfunc: sp->qnmlig)
                            {
                                int nb = qnm_ival(qfunc.first);
                                int mb = qnm_jval(qfunc.first);
                                int l = qnm_lval(qfunc.first);
                                double qrad = AtomicInterpolateInline_1 (qfunc.second.data(), ic, d0, d1, dm);
                                Atoms[ion].grid_qr[qnm_key(nb, mb, l)][icount] = qrad;
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

        ct.q_alloc[0] += (size_t)Atoms[ion].grid_ylm.size() * icount * sizeof(double) +
                         (size_t)Atoms[ion].grid_qr.size() * icount * sizeof(float);

    }                           /*end for ion */

    // Sum q-function memory usage over all nodes
    size_t talloc;
    MPI_Allreduce(&ct.q_alloc[0], &talloc, 1, MPI_LONG, MPI_MIN, pct.grid_comm);
    MPI_Allreduce(&talloc, &ct.q_alloc[1], 1, MPI_LONG, MPI_MIN, pct.kpsub_comm);

    MPI_Allreduce(&ct.q_alloc[0], &talloc, 1, MPI_LONG, MPI_MAX, pct.grid_comm);
    MPI_Allreduce(&talloc, &ct.q_alloc[2], 1, MPI_LONG, MPI_MAX, pct.kpsub_comm);

    MPI_Allreduce(MPI_IN_PLACE, &ct.q_alloc[0], 1, MPI_LONG, MPI_SUM, pct.grid_comm);
    MPI_Allreduce(MPI_IN_PLACE, &ct.q_alloc[0], 1, MPI_LONG, MPI_SUM, pct.kpsub_comm);


    delete [](Aiz);
    delete [](Aiy);
    delete [](Aix);
    delete [](pvec);
    delete []lpx;
    delete []lpl;
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

static inline void map_qval_components(int ih, int jh, int *lpx, int *lpl, SPECIES * sp, std::multimap<size_t, qongrid> &aug_desc)
{
    int lp, l;

    int nb = sp->indv[ih];
    int mb = sp->indv[jh];

    int num_lm = (ct.max_l +1 ) * (ct.max_l +1);
    int num_LM2 = (2*ct.max_l+1) *(2*ct.max_l +1);

    int ivl = sp->nh_l2m[ih];
    int jvl = sp->nh_l2m[jh];

    for (int lm = 0; lm < lpx[ivl*num_lm + jvl]; lm++)
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

        int cg_idx = lp*num_lm * num_lm + ivl * num_lm + jvl;
        int ylm_idx = lp;

//printf("Debug  %d  %d  %d\n",ih,jh,cg_idx);
        qongrid q;
        q.ival = ih;
        q.jval = jh;
        q.nb = nb;
        q.mb = mb;
        q.lval = l;
        q.cg_idx = cg_idx;
        q.ylm_idx = ylm_idx;
        aug_desc.insert(std::pair<size_t, qongrid>(qij_key(ih, jh), q));
    }
}
