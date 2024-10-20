/*
 *
 * Copyright 2019 The RMG Project Developers. See the COPYRIGHT file 
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

#include "Prolong.h"
#include "blas.h"
#include <complex>
#include <cmath>
//#include "transition.h"

// This is a high order interpolation/prolongation operators used when constructing the
// charge density on the high density grid.
template void Prolong::prolong (float *full, float *half, int dimx, int dimy, int dimz, int half_dimx,
                       int half_dimy, int half_dimz);
template void Prolong::prolong (double *full, double *half, int dimx, int dimy, int dimz, int half_dimx,
                       int half_dimy, int half_dimz);
template void Prolong::prolong (std::complex<double> *full, std::complex<double> *half, int dimx, int dimy, int dimz, int half_dimx,
                       int half_dimy, int half_dimz);

Prolong::Prolong(int ratio_in, int order_in, double cmix_in, TradeImages &TR_in, Lattice &L_in, BaseGrid &BG_in) : ratio(ratio_in), order(order_in), cmix(cmix_in), TR(TR_in), L(L_in), BG(BG_in)
{
    /*Order has to be even number */
    if (order % 2)
        rmg_error_handler (__FILE__, __LINE__, "This function works only for even orders.");

    ibrav = L.get_ibrav_type();

    for (int ix = 0; ix < MAX_PROLONG_ORDER; ix++)
    {
        for (int iy = 0; iy < MAX_PROLONG_ORDER; iy++)
        {
            a[ix][iy] = 0.0;
            af[ix][iy] = 0.0;
        }
    }


    for (int i = 0; i < ratio; i++)
    {
        double fraction = (double) i / (double)ratio;
        double c[MAX_PROLONG_ORDER] = {0.0};
        cgen_prolong (c, fraction);

        for (int iy = 0; iy < MAX_PROLONG_ORDER; iy++)
        {
            a[i][iy] = c[iy];
 //           rmg_printf("\n acac %d %f", iy, c[iy]);
            if(fabs(c[iy]) > 1.0e-20) af[i][iy] = (float)c[iy];
        }
    }

    c000.resize(1);
    c000[0].coeff = 1.0;
    c000[0].ix = 0;
    c000[0].iy = 0;
    c000[0].iz = 0;
    std::vector<double> fraction{0.0, 0.0, 0.0};

    fraction = {0.5,0.0,0.0};
    cgen_dist_inverse(c100, fraction);

    fraction = {0.0,0.5,0.0};
    cgen_dist_inverse(c010, fraction);

    fraction = {0.0,0.0,0.5};
    cgen_dist_inverse(c001, fraction);

    fraction = {0.5,0.5,0.0};
    cgen_dist_inverse(c110, fraction);

    fraction = {0.5,0.0,0.5};
    cgen_dist_inverse(c101, fraction);

    fraction = {0.0,0.5,0.5};
    cgen_dist_inverse(c011, fraction);

    fraction = {0.5,0.5,0.5};
    cgen_dist_inverse(c111, fraction);
}

Prolong::~Prolong(void)
{

}

void Prolong::cgen_dist_inverse(std::vector<coef_idx> &coeff_indx, std::vector<double> &fraction)
{
    int num_shells = order;
    double power_weight = 3.0;
    double xcry[3];
    double hx =  BG.get_hxgrid(1);
    double hy =  BG.get_hygrid(1);
    double hz =  BG.get_hzgrid(1);
    coef_idx one_item;
    coeff_indx.clear();
    std::vector<double> dist_list;
    for(int ix = -num_shells; ix <= num_shells; ix++)
    {
        for(int iy = -num_shells; iy <= num_shells; iy++)
        {
            for(int iz = -num_shells; iz <= num_shells; iz++)
            {
                xcry[0] = (ix - fraction[0]) * hx;
                xcry[1] = (iy - fraction[1]) * hy;
                xcry[2] = (iz - fraction[2]) * hz;
                double dist = L.metric(xcry);
                dist_list.push_back(dist);
            }
        }
    }
    std::sort(dist_list.begin(), dist_list.end() );
    double max_dist = 0.0;
    int count_shells = 0;
    for(auto it = dist_list.begin(); it != dist_list.end(); ++it)
    {
 //       rmg_printf("\n dist %f", *it);
        if(*it > max_dist + 1.0e-5)
        {
            count_shells++;
            max_dist = *it;
  //          rmg_printf("\n maxdist %d %f",count_shells, max_dist);
        }
        if(count_shells == num_shells)
        {
            break;
        }
    }

    double tot_weight = 0.0;
    for(int ix = -num_shells; ix <= num_shells; ix++)
    {
        for(int iy = -num_shells; iy <= num_shells; iy++)
        {
            for(int iz = -num_shells; iz <= num_shells; iz++)
            {
                xcry[0] = (ix - fraction[0]) * hx;
                xcry[1] = (iy - fraction[1]) * hy;
                xcry[2] = (iz - fraction[2]) * hz;
                double dist = L.metric(xcry);
                if( dist - max_dist <1.0e-5)
                {
 //                   rmg_printf("\n dist %f index %d %d %d", dist, ix, iy, iz);
                    double weight = std::pow(dist, -power_weight);
                    tot_weight += weight;
                    one_item.coeff = weight;
                    one_item.ix = ix;
                    one_item.iy = iy;
                    one_item.iz = iz;
                    coeff_indx.push_back(one_item);
                }
            }
        }
    }


    //rmg_printf("\n size of coef %d", coeff_indx.size());
    for(auto it = coeff_indx.begin(); it != coeff_indx.end(); ++it)
    {
        it->coeff /= tot_weight;
  //      rmg_printf("\n aaa %f %d %d %d", it->coeff, it->ix, it->iy, it->iz);
    }
}

void Prolong::cgen_prolong (double *coef, double fraction)
{

    if(order == 0 ) return;
    double *A = new double[order * order]();
    double *b = new double[order]();
    double *d = new double[order]();
    int *ipvt = new int[order]();
    int info;
    int ione = 1;

    // cmix is the mixing coefficient from order and order-2

    /*filling A , b and d (d is distance from different coarse grids to the interpolated pt, 
      coarse grid spacing is normalized to be ONE) */
    b[0] = 1.0;

    for (int iy = 0; iy < order; iy++)
    {
        d[iy] = iy + 1.0 - fraction - (double) order / 2;
        for (int ix = 0; ix < order; ix++)
        {
            A[iy * order + ix] = pow (d[iy], ix);
        }                       /* end for */
    }                           /* end for */

    /*  solving Ac=b for c using  b = A^(-1) * b  */
    dgesv (&order, &ione, A, &order, ipvt, b, &order, &info);

    if(order == 2) cmix = 0.0;
    for (int ix = 0; ix < order; ix++) coef[ix] = (1+cmix) * b[ix];
    if(order == 2) return;

    int order2 = order-2;

    for(int i = 0; i < order2; i++) b[i] = 0.0;
    b[0] = 1.0;

    for (int iy = 0; iy < order2; iy++)
    {
        d[iy] = iy + 1.0 - fraction - (double) order2 / 2;
        for (int ix = 0; ix < order2; ix++)
        {
            A[iy * order2 + ix] = pow (d[iy], ix);
        }                       /* end for */
    }                           /* end for */

    /*  solving Ac=b for c using  b = A^(-1) * b  */
    dgesv (&order2, &ione, A, &order2, ipvt, b, &order2, &info);

    for (int ix = 0; ix < order2; ix++) coef[ix+1] =  coef[ix+1] - cmix * b[ix];


    delete [] A;
    delete [] ipvt;
    delete [] d;
    delete [] b;

}




template <typename T> void Prolong::prolong (T *full, T *half, int dimx, int dimy, int dimz, int half_dimx,
        int half_dimy, int half_dimz)
{

    if(0)
    {
        prolong_any (full, half, dimx, dimy, dimz, half_dimx, half_dimy, half_dimz);
        return;
    }
    if(ibrav == HEXAGONAL && ratio == 2)
    {
        prolong_hex2 (full, half, dimx, dimy, dimz, half_dimx, half_dimy, half_dimz);
        return;
    }

    if(ibrav == HEXAGONAL2 && ratio == 2)
    {
        prolong_hex2a (full, half, dimx, dimy, dimz, half_dimx, half_dimy, half_dimz);
        return;
    }

    if( ibrav == CUBIC_BC && ratio == 2)
    {
        prolong_bcc (full, half, dimx, dimy, dimz, half_dimx, half_dimy, half_dimz);
        //prolong_bcc_other (full, half, dimx, dimy, dimz, half_dimx, half_dimy, half_dimz);
        return;
    }
    if( ibrav == CUBIC_FC && ratio == 2)
    {
        prolong_fcc (full, half, dimx, dimy, dimz, half_dimx, half_dimy, half_dimz);
        return;
    }

    T *sg_half = new T[(half_dimx + order) * (half_dimy + order) * (half_dimz + order)];

    TR.trade_imagesx (half, sg_half, half_dimx, half_dimy, half_dimz, order/2, FULL_TRADE);

    int incy = dimz / ratio + order;
    int incx = (dimz / ratio + order) * (dimy / ratio + order);

    int incy2 = dimz;
    int incx2 = dimz * dimy;

    int incx3 = (dimz / ratio + order) * dimy;

    // Optimized most common case
    if(order == 10 && ratio == 2)
    {
        T* fulla0 = new T[(dimy / ratio + order) * (dimz / ratio + order)];
        T* fulla1 = new T[(dimy / ratio + order) * (dimz / ratio + order)];
        T* fullb0 = new T[dimy * (dimz / ratio + order)];
        T* fullb1 = new T[dimy * (dimz / ratio + order)];

        for (int ix = 0; ix < dimx / 2; ix++)
        {
            for (int iy = 0; iy < dimy / 2 + 10; iy++)
            {
                for (int iz = 0; iz < dimz / 2 + 10; iz++)
                {
                    T sum = 0.0;
                    T *halfptr = &sg_half[(ix + 1) * incx + iy * incy + iz];
                    for(int k = 0;k < 10;k++) sum+= a[0][k] * halfptr[k*incx];
                    fulla0[iy * incy + iz] = sum;
                    sum = 0.0;
                    for(int k = 0;k < 10;k++) sum+= a[1][k] * halfptr[k*incx];
                    fulla1[iy * incy + iz] = sum;
                }
            }

            for (int iy = 0; iy < dimy / 2; iy++)
            {
                for (int iz = 0; iz < dimz / 2 + 10; iz++)
                {
                    T sum = 0.0;
                    T *full_tmp = &fulla0[(iy + 1) * incy + iz];

                    for(int k = 0;k < 10;k++) sum+= a[0][k] * full_tmp[k*incy];
                    fullb0[(2 * iy + 0) * incy + iz] = sum;

                    sum = 0.0;
                    for(int k = 0;k < 10;k++) sum+= a[1][k] * full_tmp[k*incy];
                    fullb0[(2 * iy + 1) * incy + iz] = sum;

                    sum = 0.0;
                    full_tmp = &fulla1[(iy + 1) * incy + iz];
                    for(int k = 0;k < 10;k++) sum+= a[0][k] * full_tmp[k*incy];
                    fullb1[(2 * iy + 0) * incy + iz] = sum;

                    sum = 0.0;
                    for(int k = 0;k < 10;k++) sum+= a[1][k] * full_tmp[k*incy];
                    fullb1[(2 * iy + 1) * incy + iz] = sum;

                }
            }

            for (int iy = 0; iy < dimy; iy++)
            {
                for (int iz = 0; iz < dimz / 2; iz++)
                {
                    T sum = 0.0;
                    T *full_tmp = &fullb0[iy * incy + iz + 1];
                    for(int k = 0;k < 10;k++) sum+= a[0][k] * full_tmp[k];
                    full[2*ix * incx2 + iy * incy2 + 2 * iz + 0] = sum;

                    sum = 0.0;
                    for(int k = 0;k < 10;k++) sum+= a[1][k] * full_tmp[k];
                    full[2*ix * incx2 + iy * incy2 + 2 * iz + 1] = sum;

                    sum = 0.0;
                    full_tmp = &fullb1[iy * incy + iz + 1];
                    for(int k = 0;k < 10;k++) sum+= a[0][k] * full_tmp[k];
                    full[(2*ix + 1) * incx2 + iy * incy2 + 2 * iz + 0] = sum;

                    sum = 0.0;
                    for(int k = 0;k < 10;k++) sum+= a[1][k] * full_tmp[k];
                    full[(2*ix + 1) * incx2 + iy * incy2 + 2 * iz + 1] = sum;

                }
            }
        }


        delete [] fullb1;
        delete [] fullb0;
        delete [] fulla1;
        delete [] fulla0;
        delete [] sg_half;
        return;
    }

    T* fulla = new T[dimx * (dimy / ratio + order) * (dimz / ratio + order)];
    T* fullb = new T[ dimx * dimy * (dimz / ratio + order)];

    for (int i = 0; i < ratio; i++)
    {
        for (int ix = 0; ix < dimx / ratio; ix++)
        {
            for (int iy = 0; iy < dimy / ratio + order; iy++)
            {
                for (int iz = 0; iz < dimz / ratio + order; iz++)
                {
                    T sum = 0.0;
                    T *halfptr = &sg_half[(ix + 1) * incx + iy * incy + iz];
                    for(int k = 0;k < order;k++) sum+= a[i][k] * halfptr[k*incx];
                    fulla[((ratio * ix) + i) * incx + iy * incy + iz] = sum;
                }
            }
        }
    }


    for (int i = 0; i < ratio; i++)
    {
        for (int ix = 0; ix < dimx; ix++)
        {
            for (int iy = 0; iy < dimy / ratio; iy++)
            {
                for (int iz = 0; iz < dimz / ratio + order; iz++)
                {
                    T sum = 0.0;
                    T *full_tmp = &fulla[ix * incx + (iy + 1) * incy + iz];
                    for(int k = 0;k < order;k++) sum+= a[i][k] * full_tmp[k*incy];
                    fullb[ix * incx3 + (ratio * iy + i) * incy + iz] = sum;
                }
            }
        }
    }


    for (int i = 0; i < ratio; i++)
    {
        for (int ix = 0; ix < dimx; ix++)
        {
            for (int iy = 0; iy < dimy; iy++)
            {
                for (int iz = 0; iz < dimz / ratio; iz++)
                {
                    T sum = 0.0;
                    T *full_tmp = &fullb[ix * incx3 + iy * incy + iz + 1];
                    for(int k = 0;k < order;k++) sum+= a[i][k] * full_tmp[k];
                    full[ix * incx2 + iy * incy2 + ratio * iz + i] = sum;
                }
            }
        }
    }


    delete [] fulla;
    delete [] fullb;
    delete [] sg_half;

}




template void Prolong::prolong_hex2 (float *full, float *half, int dimx, int dimy, int dimz, int half_dimx,
        int half_dimy, int half_dimz);
template void Prolong::prolong_hex2 (double *full, double *half, int dimx, int dimy, int dimz, int half_dimx,
        int half_dimy, int half_dimz);
template void Prolong::prolong_hex2 (std::complex<double> *full, std::complex<double> *half, int dimx, int dimy, int dimz, int half_dimx,
        int half_dimy, int half_dimz);
template <typename T> void Prolong::prolong_hex2 (T *full, T *half, int dimx, int dimy, int dimz, int half_dimx,
        int half_dimy, int half_dimz)
{

    if(ratio != 2)
        rmg_error_handler (__FILE__, __LINE__, "This function works only for 2-times fine grid.");

    T *sg_half = new T[(half_dimx + order) * (half_dimy + order) * (half_dimz + order)];

    TR.trade_imagesx (half, sg_half, half_dimx, half_dimy, half_dimz, order/2, FULL_TRADE);

    int incy = dimz / ratio + order;
    int incx = (dimz / ratio + order) * (dimy / ratio + order);


    int incx2 = (dimz / ratio + order) * dimy;
    int incy2 = dimz / ratio + order;

    // Optimized most common case

    T* fulla = new T[dimx * dimy * (dimz / ratio + order)];

    for (int ix = 0; ix < dimx / ratio; ix++)
    {
        for (int iy = 0; iy < dimy / ratio; iy++)
        {
            for (int iz = 0; iz < dimz / ratio + order; iz++)
            {
                T sum = 0.0;
                fulla[(ratio * ix) * incx2 + ratio * iy * incy2 + iz] = sg_half[(ix+order/2) * incx + (iy+order/2) * incy + iz];

                sum = 0.0;
                for(int k = 0;k < order;k++) sum+= a[1][k] * sg_half[(ix+k+1)*incx + (iy+order/2) * incy + iz];
                fulla[((ratio * ix) + 1) * incx2 + ratio * iy * incy2 + iz] = sum;

                sum = 0.0;
                for(int k = 0;k < order;k++) sum+= a[1][k] * sg_half[(ix+order/2)*incx + (iy+k+1) * incy + iz];
                fulla[(ratio * ix) * incx2 + (ratio * iy +1) * incy2 + iz] = sum;

                sum = 0.0;
                for(int k = 0;k < order;k++) sum+= a[1][k] * sg_half[(ix+k+1)*incx + (iy+k+1) * incy + iz];
                fulla[((ratio * ix)+1) * incx2 + (ratio * iy +1) * incy2 + iz] = sum;
            }
        }
    }


    for (int i = 0; i < ratio; i++)
    {
        for (int ix = 0; ix < dimx; ix++)
        {
            for (int iy = 0; iy < dimy; iy++)
            {
                for (int iz = 0; iz < dimz / ratio; iz++)
                {
                    T sum = 0.0;
                    T *full_tmp = &fulla[ix * incx2 + iy * incy + iz + 1];
                    for(int k = 0;k < order;k++) sum+= a[i][k] * full_tmp[k];
                    full[ix * dimy * dimz + iy * dimz + ratio * iz + i] = sum;
                }
            }
        }
    }


    delete [] fulla;
    delete [] sg_half;

}

template void Prolong::prolong_hex2a (float *full, float *half, int dimx, int dimy, int dimz, int half_dimx,
        int half_dimy, int half_dimz);
template void Prolong::prolong_hex2a (double *full, double *half, int dimx, int dimy, int dimz, int half_dimx,
        int half_dimy, int half_dimz);
template void Prolong::prolong_hex2a (std::complex<double> *full, std::complex<double> *half, int dimx, int dimy, int dimz, int half_dimx,
        int half_dimy, int half_dimz);
template <typename T> void Prolong::prolong_hex2a (T *full, T *half, int dimx, int dimy, int dimz, int half_dimx,
        int half_dimy, int half_dimz)
{

    if(ratio != 2)
        rmg_error_handler (__FILE__, __LINE__, "This function works only for 2-times fine grid.");

    T *sg_half = new T[(half_dimx + order) * (half_dimy + order) * (half_dimz + order)];

    TR.trade_imagesx (half, sg_half, half_dimx, half_dimy, half_dimz, order/2, FULL_TRADE);

    int incy = dimz / ratio + order;
    int incx = (dimz / ratio + order) * (dimy / ratio + order);


    int incx2 = (dimz / ratio + order) * dimy;
    int incy2 = dimz / ratio + order;

    // Optimized most common case

    T* fulla = new T[dimx * dimy * (dimz / ratio + order)];

    for (int ix = 0; ix < dimx / ratio; ix++)
    {
        for (int iy = 0; iy < dimy / ratio; iy++)
        {
            for (int iz = 0; iz < dimz / ratio + order; iz++)
            {
                T sum = 0.0;
                fulla[(ratio * ix) * incx2 + ratio * iy * incy2 + iz] = sg_half[(ix+order/2) * incx + (iy+order/2) * incy + iz];

                sum = 0.0;
                for(int k = 0;k < order;k++) sum+= a[1][k] * sg_half[(ix+k+1)*incx + (iy+order/2) * incy + iz];
                fulla[((ratio * ix) + 1) * incx2 + ratio * iy * incy2 + iz] = sum;

                sum = 0.0;
                for(int k = 0;k < order;k++) sum+= a[1][k] * sg_half[(ix+order/2)*incx + (iy+k+1) * incy + iz];
                fulla[(ratio * ix) * incx2 + (ratio * iy +1) * incy2 + iz] = sum;

                sum = 0.0;
                for(int k = 0;k < order;k++) sum+= a[1][k] * sg_half[(ix+k+1)*incx + (iy + order -k) * incy + iz];
                fulla[((ratio * ix)+1) * incx2 + (ratio * iy +1) * incy2 + iz] = sum;
            }
        }
    }


    for (int i = 0; i < ratio; i++)
    {
        for (int ix = 0; ix < dimx; ix++)
        {
            for (int iy = 0; iy < dimy; iy++)
            {
                for (int iz = 0; iz < dimz / ratio; iz++)
                {
                    T sum = 0.0;
                    T *full_tmp = &fulla[ix * incx2 + iy * incy + iz + 1];
                    for(int k = 0;k < order;k++) sum+= a[i][k] * full_tmp[k];
                    full[ix * dimy * dimz + iy * dimz + ratio * iz + i] = sum;
                }
            }
        }
    }


    delete [] fulla;
    delete [] sg_half;

}

template void Prolong::prolong_bcc (float *full, float *half, int dimx, int dimy, int dimz, int half_dimx,
        int half_dimy, int half_dimz);
template void Prolong::prolong_bcc (double *full, double *half, int dimx, int dimy, int dimz, int half_dimx,
        int half_dimy, int half_dimz);
template void Prolong::prolong_bcc (std::complex<double> *full, std::complex<double> *half, int dimx, int dimy, int dimz, int half_dimx,
        int half_dimy, int half_dimz);
template <typename T> void Prolong::prolong_bcc (T *full, T *half, int dimx, int dimy, int dimz, int half_dimx,
        int half_dimy, int half_dimz)
{

    if(ratio != 2)
        rmg_error_handler (__FILE__, __LINE__, "This function works only for 2-times fine grid.");

    // need to have one more points each side, order 10, each side need 6 points
    T *sg_half = new T[(half_dimx + order) * (half_dimy + order) * (half_dimz + order)];

    TR.trade_imagesx (half, sg_half, half_dimx, half_dimy, half_dimz, order/2, FULL_TRADE);

    int incy = dimz / ratio + order ;
    int incx = (dimz / ratio + order ) * (dimy / ratio + order);

    int order1 = order/2 ;

    int incx2 = dimz * dimy;
    int incy2 = dimz;

    // Optimized most common case

    for (int ix = 0; ix < dimx / ratio; ix++)
    {
        for (int iy = 0; iy < dimy / ratio; iy++)
        {
            for (int iz = 0; iz < dimz / ratio; iz++)
            {
                // original grid
                T sum = 0.0;
                int idx_full = ratio *(ix * incx2 +iy *incy2 + iz);
                int idx_half =(ix+order1) * incx + (iy+order1) *incy + iz + order1;
                full[idx_full] = sg_half[idx_half];
                // +1 grid in lattice a direction
                sum = 0.0;
                for(int k = -order/2+1;k < order/2+1;k++) sum+= a[1][k+order/2-1] * sg_half[idx_half + k*incx];
                full[idx_full + incx2] = sum;

                // +1 grid in lattice b direction
                sum = 0.0;
                for(int k = -order/2+1;k < order/2+1;k++) sum+= a[1][k+order/2-1] * sg_half[idx_half + k*incy];
                full[idx_full + incy2] = sum;

                // +1 grid in lattice c direction
                sum = 0.0;
                for(int k = -order/2+1;k < order/2+1;k++) sum+= a[1][k+order/2-1] * sg_half[idx_half + k];
                full[idx_full + 1] = sum;

                // +1 grid in lattice a and b directions
                sum = 0.0;
                for(int k = -order/2+1;k < order/2+1;k++) sum+= a[1][k+order/2-1] * sg_half[idx_half + k * (incx + incy)];
                full[idx_full + incx2 + incy2] = sum ;

                // +1 grid in lattice a and c directions
                sum = 0.0;
                for(int k = -order/2+1;k < order/2+1;k++) sum+= a[1][k+order/2-1] * sg_half[idx_half + k * (incx + 1) ];
                full[idx_full + incx2 + 1] = sum;

                // +1 grid in lattice b and c directions
                sum = 0.0;
                for(int k = -order/2+1;k < order/2+1;k++) sum+= a[1][k+order/2-1] * sg_half[idx_half + k * (incy + 1) ];
                full[idx_full + incy2 + 1] = sum;

                // +1 grid in lattice a,  b and c directions
                sum = 0.0;
                for(int k = -order/2+1;k < order/2+1;k++) sum+= a[1][k+order/2-1] * sg_half[idx_half + k*(incx + incy + 1)];
                full[idx_full + incx2 + incy2 + 1] = sum;
            }
        }
    }


    delete [] sg_half;

}

template void Prolong::prolong_bcc_other (float *full, float *half, int dimx, int dimy, int dimz, int half_dimx,
        int half_dimy, int half_dimz);
template void Prolong::prolong_bcc_other (double *full, double *half, int dimx, int dimy, int dimz, int half_dimx,
        int half_dimy, int half_dimz);
template void Prolong::prolong_bcc_other (std::complex<double> *full, std::complex<double> *half, int dimx, int dimy, int dimz, int half_dimx,
        int half_dimy, int half_dimz);
template <typename T> void Prolong::prolong_bcc_other (T *full, T *half, int dimx, int dimy, int dimz, int half_dimx,
        int half_dimy, int half_dimz)
{

    if(ratio != 2)
        rmg_error_handler (__FILE__, __LINE__, "This function works only for 2-times fine grid.");

    // need to have one more points each side, order 10, each side need 6 points
    T *sg_half = new T[(dimx + 2*order) * (dimy + 2*order) * (dimz + 2*order)];

    TR.trade_imagesx (half, sg_half, half_dimx, half_dimy, half_dimz, order/2, FULL_TRADE);

    int incy = dimz / ratio + order ;
    int incx = (dimz / ratio + order ) * (dimy / ratio + order);

    int order1 = order/2 ;

    int incx2 = dimz * dimy;
    int incy2 = dimz;

    // Optimized most common case

    for(int i = 0; i < dimx*dimy*dimz; i++) full[i] = 0.0;
    for (int ix = 0; ix < dimx / ratio; ix++)
    {
        for (int iy = 0; iy < dimy / ratio; iy++)
        {
            for (int iz = 0; iz < dimz / ratio; iz++)
            {
                // original grid
                T sum = 0.0;
                int idx_full = ratio *(ix * incx2 +iy *incy2 + iz);
                int idx_half =(ix+order1) * incx + (iy+order1) *incy + iz + order1;
                full[idx_full] = sg_half[idx_half];
                // +1 grid in lattice a direction
                sum = 0.0;
                for(int k = -order/2+1;k < order/2+1;k++) sum+= a[1][k+order/2-1] * sg_half[idx_half + k*incx];
                full[idx_full + incx2] = sum;

                // +1 grid in lattice b direction
                sum = 0.0;
                for(int k = -order/2+1;k < order/2+1;k++) sum+= a[1][k+order/2-1] * sg_half[idx_half + k*incy];
                full[idx_full + incy2] = sum;

                // +1 grid in lattice c direction
                sum = 0.0;
                for(int k = -order/2+1;k < order/2+1;k++) sum+= a[1][k+order/2-1] * sg_half[idx_half + k];
                full[idx_full + 1] = sum;

                // +1 grid in lattice a,  b and c directions
                sum = 0.0;
                for(int k = -order/2+1;k < order/2+1;k++) sum+= a[1][k+order/2-1] * sg_half[idx_half + k*(incx + incy + 1)];
                full[idx_full + incx2 + incy2 + 1] = sum;
            }
        }
    }

    TR.trade_imagesx (full, sg_half, dimx, dimy, dimz, order, FULL_TRADE);

    incy = dimz + 2*order ;
    incx = (dimz + 2*order ) * (dimy  + 2*order);
    int incz = 1;
    for (int ix = 0; ix < dimx / ratio; ix++)
    {
        for (int iy = 0; iy < dimy / ratio; iy++)
        {
            for (int iz = 0; iz < dimz / ratio; iz++)
            {
                int idx_full = ratio *(ix * incx2 +iy *incy2 + iz);
                int idx_full_s110 = ratio *(ix * incx +iy *incy + iz) + order *(incx + incy + incz) + incx + incy;
                int idx_full_s101 = ratio *(ix * incx +iy *incy + iz) + order *(incx + incy + incz) + incx + incz;
                int idx_full_s011 = ratio *(ix * incx +iy *incy + iz) + order *(incx + incy + incz) + incy + incz;

                // +1 grid in lattice a and b directions
                T sum = 0.0;
                for(int k = -order/2+1;k < order/2+1;k++) sum+= a[1][k+order/2-1] * sg_half[idx_full_s110 + (2*k-1)*incx];
                for(int k = -order/2+1;k < order/2+1;k++) sum+= a[1][k+order/2-1] * sg_half[idx_full_s110 + (2*k-1)*incy];
                for(int k = -order/2+1;k < order/2+1;k++) sum+= a[1][k+order/2-1] * sg_half[idx_full_s110 + (2*k-1)*incz];
                for(int k = -order/2+1;k < order/2+1;k++) sum+= a[1][k+order/2-1] * sg_half[idx_full_s110 + (2*k-1)*(incx + incy +incz)];
                full[idx_full + incx2 + incy2] = sum/4.0;

                // +1 grid in lattice a and c directions
                sum = 0.0;
                for(int k = -order/2+1;k < order/2+1;k++) sum+= a[1][k+order/2-1] * sg_half[idx_full_s101 + (2*k-1)*incx];
                for(int k = -order/2+1;k < order/2+1;k++) sum+= a[1][k+order/2-1] * sg_half[idx_full_s101 + (2*k-1)*incy];
                for(int k = -order/2+1;k < order/2+1;k++) sum+= a[1][k+order/2-1] * sg_half[idx_full_s101 + (2*k-1)*incz];
                for(int k = -order/2+1;k < order/2+1;k++) sum+= a[1][k+order/2-1] * sg_half[idx_full_s101 + (2*k-1)*(incx + incy +incz)];
                full[idx_full + incx2 + 1] = sum/4.0;

                // +1 grid in lattice b and c directions
                sum = 0.0;
                for(int k = -order/2+1;k < order/2+1;k++) sum+= a[1][k+order/2-1] * sg_half[idx_full_s011 + (2*k-1)*incx];
                for(int k = -order/2+1;k < order/2+1;k++) sum+= a[1][k+order/2-1] * sg_half[idx_full_s011 + (2*k-1)*incy];
                for(int k = -order/2+1;k < order/2+1;k++) sum+= a[1][k+order/2-1] * sg_half[idx_full_s011 + (2*k-1)*incz];
                for(int k = -order/2+1;k < order/2+1;k++) sum+= a[1][k+order/2-1] * sg_half[idx_full_s011 + (2*k-1)*(incx + incy +incz)];
                full[idx_full + incy2 + 1] = sum/4.0;
            }
        }
    }

    delete [] sg_half;

}


template void Prolong::prolong_any (float *full, float *half, int dimx, int dimy, int dimz, int half_dimx,
        int half_dimy, int half_dimz);
template void Prolong::prolong_any (double *full, double *half, int dimx, int dimy, int dimz, int half_dimx,
        int half_dimy, int half_dimz);
template void Prolong::prolong_any (std::complex<double> *full, std::complex<double> *half, int dimx, int dimy, int dimz, int half_dimx,
        int half_dimy, int half_dimz);
template <typename T> void Prolong::prolong_any (T *full, T *half, int dimx, int dimy, int dimz, int half_dimx,
        int half_dimy, int half_dimz)
{

    if(ratio != 2)
        rmg_error_handler (__FILE__, __LINE__, "This function works only for 2-times fine grid.");

    // need to have one more points each side, order 10, each side need 6 points
    T *sg_half = new T[(half_dimx + order) * (half_dimy + order) * (half_dimz + order)];

    TR.trade_imagesx (half, sg_half, half_dimx, half_dimy, half_dimz, order/2, FULL_TRADE);

    int incy = dimz / ratio + order ;
    int incx = (dimz / ratio + order ) * (dimy / ratio + order);

    int order1 = order/2 ;

    int incx2 = dimz * dimy;
    int incy2 = dimz;

    // Optimized most common case

    for (int ix = 0; ix < dimx / ratio; ix++)
    {
        for (int iy = 0; iy < dimy / ratio; iy++)
        {
            for (int iz = 0; iz < dimz / ratio; iz++)
            {
                // original grid
                T sum = 0.0;
                int idx_full = ratio *(ix * incx2 +iy *incy2 + iz);
                int idx_half =(ix+order1) * incx + (iy+order1) *incy + iz + order1;
                full[idx_full] = sg_half[idx_half];
                // +1 grid in lattice a direction
                sum = 0.0;
                for(auto it = c100.begin(); it != c100.end(); ++it)
                {
                    sum+=  it->coeff * sg_half[idx_half + it->ix*incx + it->iy * incy + it->iz];
                }
                full[idx_full + incx2] = sum;

                // +1 grid in lattice b direction
                sum = 0.0;
                for(auto it = c010.begin(); it != c010.end(); ++it)
                {
                    sum+=  it->coeff * sg_half[idx_half + it->ix*incx + it->iy * incy + it->iz];
                }
                full[idx_full + incy2] = sum;

                // +1 grid in lattice c direction
                sum = 0.0;
                for(auto it = c001.begin(); it != c001.end(); ++it)
                {
                    sum+=  it->coeff * sg_half[idx_half + it->ix*incx + it->iy * incy + it->iz];
                }
                full[idx_full + 1] = sum;

                // +1 grid in lattice a and b directions
                sum = 0.0;
                for(auto it = c110.begin(); it != c110.end(); ++it)
                {
                    sum+=  it->coeff * sg_half[idx_half + it->ix*incx + it->iy * incy + it->iz];
                }
                full[idx_full + incx2 + incy2] = sum ;

                // +1 grid in lattice a and c directions
                sum = 0.0;
                for(auto it = c101.begin(); it != c101.end(); ++it)
                {
                    sum+=  it->coeff * sg_half[idx_half + it->ix*incx + it->iy * incy + it->iz];
                }
                full[idx_full + incx2 + 1] = sum;

                // +1 grid in lattice b and c directions
                sum = 0.0;
                for(auto it = c011.begin(); it != c011.end(); ++it)
                {
                    sum+=  it->coeff * sg_half[idx_half + it->ix*incx + it->iy * incy + it->iz];
                }
                full[idx_full + incy2 + 1] = sum;

                // +1 grid in lattice a,  b and c directions
                sum = 0.0;
                for(auto it = c111.begin(); it != c111.end(); ++it)
                {
                    sum+=  it->coeff * sg_half[idx_half + it->ix*incx + it->iy * incy + it->iz];
                }
                full[idx_full + incx2 + incy2 + 1] = sum;
            }
        }
    }


    delete [] sg_half;

}

template void Prolong::prolong_fcc (float *full, float *half, int dimx, int dimy, int dimz, int half_dimx,
        int half_dimy, int half_dimz);
template void Prolong::prolong_fcc (double *full, double *half, int dimx, int dimy, int dimz, int half_dimx,
        int half_dimy, int half_dimz);
template void Prolong::prolong_fcc (std::complex<double> *full, std::complex<double> *half, int dimx, int dimy, int dimz, int half_dimx,
        int half_dimy, int half_dimz);
template <typename T> void Prolong::prolong_fcc (T *full, T *half, int dimx, int dimy, int dimz, int half_dimx,
        int half_dimy, int half_dimz)
{

    if(ratio != 2)
        rmg_error_handler (__FILE__, __LINE__, "This function works only for 2-times fine grid.");

    // need to have one more points each side, order 10, each side need 6 points
    T *sg_half = new T[(half_dimx + order) * (half_dimy + order) * (half_dimz + order)];

    TR.trade_imagesx (half, sg_half, half_dimx, half_dimy, half_dimz, order/2, FULL_TRADE);

    int incy = dimz / ratio + order ;
    int incx = (dimz / ratio + order ) * (dimy / ratio + order);

    int order1 = order/2 ;

    int incx2 = dimz * dimy;
    int incy2 = dimz;

    // Optimized most common case

    for (int ix = 0; ix < dimx / ratio; ix++)
    {
        for (int iy = 0; iy < dimy / ratio; iy++)
        {
            for (int iz = 0; iz < dimz / ratio; iz++)
            {
                // original grid
                T sum = 0.0;
                int idx_full = ratio *(ix * incx2 +iy *incy2 + iz);
                int idx_half =(ix+order1) * incx + (iy+order1) *incy + iz + order1;
                full[idx_full] = sg_half[idx_half];
                // +1 grid in lattice a direction
                sum = 0.0;
                for(int k = -order/2+1;k < order/2+1;k++) sum+= a[1][k+order/2-1] * sg_half[idx_half + k*incx];
                full[idx_full + incx2] = sum;

                // +1 grid in lattice b direction
                sum = 0.0;
                for(int k = -order/2+1;k < order/2+1;k++) sum+= a[1][k+order/2-1] * sg_half[idx_half + k*incy];
                full[idx_full + incy2] = sum;

                // +1 grid in lattice c direction
                sum = 0.0;
                for(int k = -order/2+1;k < order/2+1;k++) sum+= a[1][k+order/2-1] * sg_half[idx_half + k];
                full[idx_full + 1] = sum;

                // +1 grid in lattice a and b directions
                sum = 0.0;
                for(int k = -order/2+1;k < order/2+1;k++) sum+= a[1][k+order/2-1] * sg_half[idx_half + k * incx + (-k+1) * incy];
                full[idx_full + incx2 + incy2] = sum ;

                // +1 grid in lattice a and c directions
                sum = 0.0;
                for(int k = -order/2+1;k < order/2+1;k++) sum+= a[1][k+order/2-1] * sg_half[idx_half + k * incx + (-k+1) ];
                full[idx_full + incx2 + 1] = sum;

                // +1 grid in lattice b and c directions
                sum = 0.0;
                for(int k = -order/2+1;k < order/2+1;k++) sum+= a[1][k+order/2-1] * sg_half[idx_half + k * incy + (-k+1) ];
                full[idx_full + incy2 + 1] = sum;

                // +1 grid in lattice a,  b and c directions
                sum = 0.0;
                for(int k = -order/2+1;k < order/2+1;k++) sum+= a[1][k+order/2-1] * sg_half[idx_half + k*(incx + incy) + (-k+1) ];
                for(int k = -order/2+1;k < order/2+1;k++) sum+= a[1][k+order/2-1] * sg_half[idx_half + k*(incx + 1) + (-k+1) * incy ];
                for(int k = -order/2+1;k < order/2+1;k++) sum+= a[1][k+order/2-1] * sg_half[idx_half + k*(incy + 1) + (-k+1) * incx ];
                full[idx_full + incx2 + incy2 + 1] = sum/3.0;
            }
        }
    }


    delete [] sg_half;

}
