/*
 *
 * Copyright (c) 2014, Emil Briggs
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.

 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
*/



#include "packfuncs.h"
#include <complex>

template void rmg::pack_stop<double>(double*, double*, int, int, int);
template void rmg::pack_stop<float>(float*, float*, int, int, int);
template void rmg::pack_stop<std::complex<float> >(std::complex<float>*, std::complex<float>*, int, int, int);
template void rmg::pack_stop<std::complex<double> >(std::complex<double>*, std::complex<double>*, int, int, int);

template void rmg::pack_ptos<double>(double*, double*, int, int, int);
template void rmg::pack_ptos<float>(float*, float*, int, int, int);
template void rmg::pack_ptos<std::complex<float> >(std::complex<float> *, std::complex<float>*, int, int, int);
template void rmg::pack_ptos<std::complex<double> >(std::complex<double> *, std::complex<double>*, int, int, int);

template void rmg::pack_ptos_convert<double>(double*, double*, int, int, int);
template void rmg::pack_ptos_convert<float>(float*, float*, int, int, int);
template void rmg::pack_ptos_convert<std::complex<float> >(std::complex<float> *, std::complex<float>*, int, int, int);
template void rmg::pack_ptos_convert<std::complex<double> >(std::complex<double> *, std::complex<double>*, int, int, int);

template void rmg::pack_stop_axpy<double>(double*, double*, double, int, int, int);
template void rmg::pack_stop_axpy<float>(float*, float*, double, int, int, int);
template void rmg::pack_stop_axpy<std::complex<float> >(std::complex<float>*, std::complex<float>*, double, int, int, int);
template void rmg::pack_stop_axpy<std::complex<double> >(std::complex<double>*, std::complex<double>*, double, int, int, int);

template <typename RmgType>
void rmg::pack_stop (RmgType * sg, RmgType * pg, int dimx, int dimy, int dimz)
{

    int ix, iy, ixh, iyh;
    int incx, incy, incxs, incys;
    incy = dimz;
    incx = dimy * dimz;
    incys = dimz + 2;
    incxs = (dimy + 2) * (dimz + 2);


    /* Transfer pg into smoothing grid */
    for (ix = 0; ix < dimx; ix++)
    {

        ixh = ix + 1;
        for (iy = 0; iy < dimy; iy++)
        {

            iyh = iy + 1;
            for(int idx = 0;idx < dimz;idx++) pg[ix * incx + iy * incy + idx] = sg[ixh * incxs + iyh * incys + 1 + idx];

        }                       /* end for */

    }                           /* end for */

}                               /* end pack_stop */

void rmg::pack_stop_convert (float * sg, double * pg, int dimx, int dimy, int dimz)
{

    int ix, iy, ixh, iyh;
    int incx, incy, incxs, incys;
    incy = dimz;
    incx = dimy * dimz;
    incys = dimz + 2;
    incxs = (dimy + 2) * (dimz + 2);


    /* Transfer pg into smoothing grid */
    for (ix = 0; ix < dimx; ix++)
    {

        ixh = ix + 1;
        for (iy = 0; iy < dimy; iy++)
        {

            iyh = iy + 1;
            for(int idx = 0;idx < dimz;idx++) pg[ix * incx + iy * incy + idx] = (double)sg[ixh * incxs + iyh * incys + 1 + idx];

        }                       /* end for */

    }                           /* end for */

}                               /* end pack_stop */

void rmg::pack_stop_convert (std::complex<float> * sg, std::complex<double> * pg, int dimx, int dimy, int dimz)
{

    int ix, iy, ixh, iyh;
    int incx, incy, incxs, incys;
    incy = dimz;
    incx = dimy * dimz;
    incys = dimz + 2;
    incxs = (dimy + 2) * (dimz + 2);


    /* Transfer pg into smoothing grid */
    for (ix = 0; ix < dimx; ix++)
    {

        ixh = ix + 1;
        for (iy = 0; iy < dimy; iy++)
        {

            iyh = iy + 1;
            for(int idx = 0;idx < dimz;idx++) pg[ix * incx + iy * incy + idx] = (std::complex<double>)sg[ixh * incxs + iyh * incys + 1 + idx];

        }                       /* end for */

    }                           /* end for */

}                               /* end pack_stop */


template <typename RmgType>
void rmg::pack_ptos(RmgType * sg, RmgType * pg, int dimx, int dimy, int dimz)
{

    int ix, iy, ixh, iyh;
    int incx, incy, incxs, incys;
    incy = dimz;
    incx = dimy * dimz;
    incys = dimz + 2;
    incxs = (dimy + 2) * (dimz + 2);

    for (ix = 0; ix < (dimx + 2) * (dimy + 2) * (dimz + 2); ix++)
        sg[ix] = 0.0;

    /* Transfer pg into smoothing grid */
    for (ix = 0; ix < dimx; ix++)
    {

        ixh = ix + 1;
        for (iy = 0; iy < dimy; iy++)
        {

            iyh = iy + 1;
            for(int idx = 0;idx < dimz;idx++) sg[ixh * incxs + iyh * incys + 1 + idx] = pg[ix * incx + iy * incy + idx];

        }                       /* end for */

    }                           /* end for */

}                               /* end pack_ptos_f */

template <typename RmgType>
void rmg::pack_ptos_convert(RmgType * sg, RmgType * pg, int dimx, int dimy, int dimz)
{

    int ix, iy, ixh, iyh;
    int incx, incy, incxs, incys;
    incy = dimz;
    incx = dimy * dimz;
    incys = dimz + 2;
    incxs = (dimy + 2) * (dimz + 2);

    for (ix = 0; ix < (dimx + 2) * (dimy + 2) * (dimz + 2); ix++)
        sg[ix] = 0.0;

    /* Transfer pg into smoothing grid */
    for (ix = 0; ix < dimx; ix++)
    {

        ixh = ix + 1;
        for (iy = 0; iy < dimy; iy++)
        {

            iyh = iy + 1;
            for(int idx = 0;idx < dimz;idx++) sg[ixh * incxs + iyh * incys + 1 + idx] = pg[ix * incx + iy * incy + idx];

        }                       /* end for */

    }                           /* end for */

}                               /* end pack_ptos_f */

void rmg::pack_ptos_convert(float * sg, double * pg, int dimx, int dimy, int dimz)
{

    int ix, iy, ixh, iyh;
    int incx, incy, incxs, incys;
    incy = dimz;
    incx = dimy * dimz;
    incys = dimz + 2;
    incxs = (dimy + 2) * (dimz + 2);

    for (ix = 0; ix < (dimx + 2) * (dimy + 2) * (dimz + 2); ix++)
        sg[ix] = 0.0;

    /* Transfer pg into smoothing grid */
    for (ix = 0; ix < dimx; ix++)
    {

        ixh = ix + 1;
        for (iy = 0; iy < dimy; iy++)
        {

            iyh = iy + 1;
            for(int idx = 0;idx < dimz;idx++) sg[ixh * incxs + iyh * incys + 1 + idx] = (float)pg[ix * incx + iy * incy + idx];

        }                       /* end for */

    }                           /* end for */

}                               /* end pack_ptos_f */

void rmg::pack_ptos_convert(std::complex<float> * sg, std::complex<double> * pg, int dimx, int dimy, int dimz)
{

    int ix, iy, ixh, iyh;
    int incx, incy, incxs, incys;
    incy = dimz;
    incx = dimy * dimz;
    incys = dimz + 2;
    incxs = (dimy + 2) * (dimz + 2);

    for (ix = 0; ix < (dimx + 2) * (dimy + 2) * (dimz + 2); ix++)
        sg[ix] = 0.0;

    /* Transfer pg into smoothing grid */
    for (ix = 0; ix < dimx; ix++)
    {

        ixh = ix + 1;
        for (iy = 0; iy < dimy; iy++)
        {

            iyh = iy + 1;
            for(int idx = 0;idx < dimz;idx++) sg[ixh * incxs + iyh * incys + 1 + idx] = (std::complex<float>)pg[ix * incx + iy * incy + idx];

        }                       /* end for */

    }                           /* end for */

}                               /* end pack_ptos_f */



template <typename RmgType>
void rmg::pack_stop_axpy (RmgType * sg, RmgType * pg, double alpha, int dimx, int dimy, int dimz)
{

    int ix, iy, iz, ixh, iyh;
    int incx, incy, incxs, incys;

    incy = dimz;
    incx = dimy * dimz;
    incys = dimz + 2;
    incxs = (dimy + 2) * (dimz + 2);
    RmgType alpha1(alpha);

    /* Transfer pg into smoothing grid */
    for (ix = 0; ix < dimx; ix++)
    {

        ixh = ix + 1;
        for (iy = 0; iy < dimy; iy++)
        {

            iyh = iy + 1;
            for (iz = 0; iz < dimz; iz++)
            {

                pg[ix * incx + iy * incy + iz] = pg[ix * incx + iy * incy + iz] + alpha1 * sg[ixh * incxs + iyh * incys + iz + 1];

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */


} // end rmg::pack_stop_axpy


/*
 * FUNCTION
 *   void pack_vhstod(double *s, double *d, int dimx, int dimy, int dimz)
 *   Used to pack grids when computing the hartree potential 
 *   For periodic system, d = s
 *   For surface system,  d[x,y,Nz/2+1:Nz/2+Nz] = s[x,y,1:Nz]
 *   FOr cluster system:  d[Nx/2+1:Nx/2+Nx, ...] = s[1:Nx, 1:Ny, 1:Nz] 
 * INPUTS
 *   s[dmix*dimy*dimz]
 *   dimx, dimy, dimz: dimension of the wave function grid
 * OUTPUT
 *   d[ct.vh_pxgrid * ct.vh_pygrid * ct.vh_pzgrid]
 *       see init.c for ct.vh_pxgrid, ...
 *       array in the hatree potential grid
 * PARENTS
 *   get_vh.c
 * CHILDREN
 *   nothing
 * SOURCE
 */

#include "rmg_grid.h"
#include "boundary_conditions.h"
#include "packfuncs.h"


/* This function is used to pack grids when computing the hartree potential */
void rmg::pack_stod (rmg::grid *G, double * s, double * d, int dimx, int dimy, int dimz, int boundaryflag)
{
    int ix, iy, iz;
    int pex, pey, pez;
    int sxlo, sylo, szlo;
    int sxhi, syhi, szhi;
    int dxlo, dylo, dzlo;
    int dxhi, dyhi, dzhi;
    int loidx, hiidx;
    int hix, hiy, hiz;
    int hincx, hincy;

    if (boundaryflag == PERIODIC)
    {
        for (ix = 0; ix < (dimx * dimy * dimz); ix++)
            d[ix] = s[ix];
        return;
    }

    G->pe2xyz (G->get_rank(), &pex, &pey, &pez);
    sxlo = pex * dimx;
    sylo = pey * dimy;
    szlo = pez * dimz;
    sxhi = sxlo + dimx - 1;
    syhi = sylo + dimy - 1;
    szhi = szlo + dimz - 1;

    dxlo = pex * 2 * dimx - G->get_NX_GRID(1) / 2;
    dylo = pey * 2 * dimy - G->get_NY_GRID(1) / 2;
    dzlo = pez * 2 * dimz - G->get_NZ_GRID(1) / 2;
    dxhi = dxlo + 2 * dimx - 1;
    dyhi = dylo + 2 * dimy - 1;
    dzhi = dzlo + 2 * dimz - 1;
    hincy = 2 * dimz;
    hincx = 4 * dimy * dimz;

    if (boundaryflag == SURFACE)
    {
        dxlo = sxlo;
        dylo = sylo;
        dxhi = sxhi;
        dyhi = syhi;
        hincx = 2 * dimy * dimz;
    } 


    loidx = 0;
    for (ix = sxlo; ix <= sxhi; ix++)
    {
        for (iy = sylo; iy <= syhi; iy++)
        {
            for (iz = szlo; iz <= szhi; iz++)
            {
                /* Check if we map into the double grid */
                if ((ix >= dxlo) && (ix <= dxhi) &&
                    (iy >= dylo) && (iy <= dyhi) && (iz >= dzlo) && (iz <= dzhi))
                {

                    hix = ix - dxlo;
                    hiy = iy - dylo;
                    hiz = iz - dzlo;
                    hiidx = hix * hincx + hiy * hincy + hiz;

                    /* Yes there is a mapping so transfer the data */
                    d[hiidx] = s[loidx];
                }
                loidx++;
            }
        }
    }
}



/*
 * NAME
 * FUNCTION
 *   void pack_vhdtos(double *s, double *d, int dimx, int dimy, int dimz)
 *   Used to pack grids when computing the hartree potential 
 *   For periodic system, s = d
 *   For surface system, s[x,y,1:Nz] = d[x,y,Nz/2+1:Nz/2+Nz)
 *   FOr cluster system: s[1:Nx, 1:Ny, 1:Nz] = d[Nx/2+1:Nx/2+Nx, ...]
 * INPUTS
 *   d[ct.vh_pxgrid * ct.vh_pygrid * ct.vh_pzgrid]
 *     see init.c for ct.vh_pxgrid, ...
 *   dimx, dimy, dimz: dimension of the wave function grid
 * OUTPUT
 *   s[dimx*dimy*dimz] array in the wave function grid
 * PARENTS
 *   get_vh.c
 * CHILDREN
 *   nothing
 * SOURCE
 */
void rmg::pack_dtos (rmg::grid *G, double * s, double * d, int dimx, int dimy, int dimz, int boundaryflag)
{
    int ix, iy, iz;
    int pex, pey, pez;
    int sxlo, sylo, szlo;
    int sxhi, syhi, szhi;
    int dxlo, dylo, dzlo;
    int dxhi, dyhi, dzhi;
    int loidx, hiidx;
    int hix, hiy, hiz;
    int hincx, hincy;

    if (boundaryflag == PERIODIC)
    {
        for (ix = 0; ix < (dimx * dimy * dimz); ix++)
            s[ix] = d[ix];
        return;
    }

    G->pe2xyz (G->get_rank(), &pex, &pey, &pez);
    sxlo = pex * dimx;
    sylo = pey * dimy;
    szlo = pez * dimz;
    sxhi = sxlo + dimx - 1;
    syhi = sylo + dimy - 1;
    szhi = szlo + dimz - 1;

    dxlo = pex * 2 * dimx - G->get_NX_GRID(1) / 2;
    dylo = pey * 2 * dimy - G->get_NY_GRID(1) / 2;
    dzlo = pez * 2 * dimz - G->get_NZ_GRID(1) / 2;
    dxhi = dxlo + 2 * dimx - 1;
    dyhi = dylo + 2 * dimy - 1;
    dzhi = dzlo + 2 * dimz - 1;
    hincy = 2 * dimz;
    hincx = 4 * dimy * dimz;

    if (boundaryflag == SURFACE)
    {
        dxlo = sxlo;
        dylo = sylo;
        dxhi = sxhi;
        dyhi = syhi;
        hincx = 2 * dimy * dimz;
    }

    loidx = 0;
    for (ix = sxlo; ix <= sxhi; ix++)
    {
        for (iy = sylo; iy <= syhi; iy++)
        {
            for (iz = szlo; iz <= szhi; iz++)
            {
                /* Check if we map into the double grid */
                if ((ix >= dxlo) && (ix <= dxhi) &&
                    (iy >= dylo) && (iy <= dyhi) && (iz >= dzlo) && (iz <= dzhi))
                {
                    hix = ix - dxlo;
                    hiy = iy - dylo;
                    hiz = iz - dzlo;
                    hiidx = hix * hincx + hiy * hincy + hiz;

                    /* Yes there is a mapping so transfer the data */
                    s[loidx] = d[hiidx];
                }
                loidx++;
            }
        }
    }
} 

