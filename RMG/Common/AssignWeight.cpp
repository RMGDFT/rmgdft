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


#include "make_conf.h"
#include "const.h"
#include "grid.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include <complex>
#include "Kpoint.h"
#include "FiniteDiff.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"

template void AssignWeight<double> (Kpoint<double> *, SPECIES *, int, fftw_complex *, double *, double *, double *);
template void AssignWeight<std::complex<double> >(Kpoint<std::complex<double>> *, SPECIES *, int, fftw_complex *, double *, std::complex<double> *, std::complex<double> *);


template <typename KpointType>
void AssignWeight (Kpoint<KpointType> *kptr, SPECIES * sp, int ion, fftw_complex * beptr, double * rtptr, KpointType *Bweight, KpointType *Nlweight)
{

    if(!pct.nl_flag[ion]) return ;
    Lattice *L = kptr->L;
    TradeImages *T = kptr->T;
    ION *iptr = &ct.ions[ion];

    BaseGrid *G = kptr->G;
    int dimx = G->get_PX0_GRID(1);
    int dimy = G->get_PY0_GRID(1);
    int dimz = G->get_PZ0_GRID(1);
    int nxgrid = G->get_NX_GRID(1);
    int nygrid = G->get_NY_GRID(1);
    int nzgrid = G->get_NZ_GRID(1);


    int ilow, jlow, klow, ihi, jhi, khi;
    find_node_offsets(pct.gridpe, nxgrid, nygrid, nzgrid,
                      &ilow, &jlow, &klow);

    ihi = ilow + dimx - 1;
    jhi = jlow + dimy - 1;
    khi = klow + dimz - 1;


    int nlxdim = sp->nldim;
    int nlydim = sp->nldim;
    int nlzdim = sp->nldim;
    if(!ct.localize_projectors) {
        nlxdim = get_NX_GRID();
        nlydim = get_NY_GRID();
        nlzdim = get_NZ_GRID();
    }

    // These are probably not right for anything but simple cubic grids
    double wx = iptr->nlcrds[0] / ct.hmaxgrid;
    double wy = iptr->nlcrds[1] / ct.hmaxgrid;
    double wz = iptr->nlcrds[2] / ct.hmaxgrid;


    int pbasis = kptr->pbasis;
    KpointType ZERO_t(0.0);

    std::complex<double> *Nlweight_C = (std::complex<double> *)Nlweight;
    std::complex<double> *Bweight_C = (std::complex<double> *)Bweight;

    double *Nlweight_R = (double *)Nlweight;
    double *Bweight_R = (double *)Bweight;


    for(int idx = 0; idx < pbasis; idx++) rtptr[idx] = 0.0;
    for(int idx = 0; idx < pbasis; idx++) Bweight[idx] = ZERO_t;
    for(int idx = 0; idx < pbasis; idx++) Nlweight[idx] = ZERO_t;
    if(pct.idxptrlen[ion] == 0) return;

    std::complex<double> *nbeptr = (std::complex<double> *)beptr;


    KpointType *tem_array = new KpointType[nlxdim * nlydim * nlzdim]();
    KpointType *Btem_array = new KpointType[nlxdim * nlydim * nlzdim]();
    std::complex<double> *tem_array_C = (std::complex<double> *)tem_array;
    std::complex<double> *Btem_array_C = (std::complex<double> *)Btem_array;
    double *Btem_array_R = (double *)Btem_array;


    for(int ix = 0; ix < nlxdim * nlydim * nlzdim; ix++) {
        tem_array[ix] = std::real(nbeptr[ix]);
        if(!ct.is_gamma) {
            tem_array_C[ix] = nbeptr[ix];
        }
    }


    // Apply B operator then map weights back
    AppCirDriverBeta (L, T, tem_array, Btem_array, nlxdim, nlydim, nlzdim, ct.kohn_sham_fd_order);

    int ixstart = iptr->nl_global_grid_xstart;
    int iystart = iptr->nl_global_grid_ystart;
    int izstart = iptr->nl_global_grid_zstart;
    int igx, igy, igz;

    int icenter = sp->nldim / 2;
    int icut = (icenter + 1) * (icenter + 1);


    for (int ix = 0; ix < nlxdim; ix++)
    {

        // global index of nl grid , mode to unit cell  
        igx = (ixstart + ix +nxgrid) % nxgrid ;

        if(igx <ilow || igx > ihi) continue;

        double w1=1.0;
        if(ix==0) w1 = 0.5*(0.5 - wx);
        if(ix==nlxdim-1) w1 = 0.5*(0.5 + wx);
        for (int iy = 0; iy < nlydim; iy++)
        {
            // global index of nl grid , mode to unit cell  
            igy = (iystart + iy +nygrid) % nygrid ;
            if(igy <jlow || igy > jhi) continue;
            double w2=1.0;
            if(iy==0) w2 = 0.5*(0.5 - wy);
            if(iy==nlydim-1) w2 = 0.5*(0.5 + wy);

            for (int iz = 0; iz < nlzdim; iz++)
            {



                int itmp = (ix - icenter) * (ix - icenter) +
                    (iy - icenter) * (iy - icenter) + (iz - icenter) * (iz - icenter);

                if (ct.localize_projectors && itmp > icut ) continue;

                // global index of nl grid , mode to unit cell  
                igz = (izstart + iz +nzgrid) % nzgrid ;
                if(igz <klow || igz > khi) continue;


                double w3 = 1.0;
                if(iz==0) w3 = 0.5*(0.5 - wz);
                if(iz==nlzdim-1) w3 = 0.5*(0.5 + wz);
                // With localized projectors we need to know if this point maps
                // into this processors space. With delocalized we know it does
                // so the index array was never even generated.
                if(!ct.localize_projectors) 
                {
                    w1 = w2 = w3 = 1.0;
                }

                int idx1 = ix * nlydim * nlzdim + iy * nlzdim + iz;
                int idx2 = (igx-ilow) * dimy * dimz +(igy - jlow) * dimz + igz-klow;
                rtptr[idx2] = std::real(nbeptr[idx1]);
                if(ct.is_gamma) {
                    Nlweight_R[idx2] = std::real(nbeptr[idx1]);
                    Bweight_R[idx2] = Btem_array_R[idx1];
                    Nlweight_R[idx2] *= w1*w2*w3;
                    Bweight_R[idx2] *= w1*w2*w3;

                }
                else {
                    Nlweight_C[idx2] = nbeptr[idx1];
                    Bweight_C[idx2] = Btem_array_C[idx1];
                    Nlweight_C[idx2] *= w1*w2*w3;
                    Bweight_C[idx2] *= w1*w2*w3;
                }

            }
        }
    }



    delete [] Btem_array;
    delete [] tem_array;

}

