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
#include "RmgTimer.h"
//#include "transition.h"

// This is a high order interpolation/prolongation operators used when constructing the
// charge density on the high density grid.
template void Prolong::prolong (float *full, float *half, int dimx, int dimy, int dimz, int half_dimx,
                       int half_dimy, int half_dimz);
template void Prolong::prolong (double *full, double *half, int dimx, int dimy, int dimz, int half_dimx,
                       int half_dimy, int half_dimz);
template void Prolong::prolong (std::complex<double> *full, std::complex<double> *half, int dimx, int dimy, int dimz, int half_dimx,
                       int half_dimy, int half_dimz);

template void Prolong::prolong<float, 6> (float *full, float *half, int half_dimx, int half_dimy, int half_dimz);
template void Prolong::prolong<float, 8> (float *full, float *half, int half_dimx, int half_dimy, int half_dimz);
template void Prolong::prolong<float, 10> (float *full, float *half, int half_dimx, int half_dimy, int half_dimz);
template void Prolong::prolong<float, 12> (float *full, float *half, int half_dimx, int half_dimy, int half_dimz);

template void Prolong::prolong<double, 6> (double *full, double *half, int half_dimx, int half_dimy, int half_dimz);
template void Prolong::prolong<double, 8> (double *full, double *half, int half_dimx, int half_dimy, int half_dimz);
template void Prolong::prolong<double, 10> (double *full, double *half, int half_dimx, int half_dimy, int half_dimz);
template void Prolong::prolong<double, 12> (double *full, double *half, int half_dimx, int half_dimy, int half_dimz);

template void Prolong::prolong_hex2<float, 6> (float *full, float *half, int half_dimx, int half_dimy, int half_dimz);
template void Prolong::prolong_hex2<float, 8> (float *full, float *half, int half_dimx, int half_dimy, int half_dimz);
template void Prolong::prolong_hex2<float, 10> (float *full, float *half, int half_dimx, int half_dimy, int half_dimz);
template void Prolong::prolong_hex2<float, 12> (float *full, float *half, int half_dimx, int half_dimy, int half_dimz);

template void Prolong::prolong_hex2<double, 6> (double *full, double *half, int half_dimx, int half_dimy, int half_dimz);
template void Prolong::prolong_hex2<double, 8> (double *full, double *half, int half_dimx, int half_dimy, int half_dimz);
template void Prolong::prolong_hex2<double, 10> (double *full, double *half, int half_dimx, int half_dimy, int half_dimz);
template void Prolong::prolong_hex2<double, 12> (double *full, double *half, int half_dimx, int half_dimy, int half_dimz);

template void Prolong::prolong<std::complex<double>, 6>
  (std::complex<double> *full, std::complex<double> *half, int half_dimx, int half_dimy, int half_dimz);
template void Prolong::prolong<std::complex<double>, 8>
  (std::complex<double> *full, std::complex<double> *half, int half_dimx, int half_dimy, int half_dimz);
template void Prolong::prolong<std::complex<double>, 10>
  (std::complex<double> *full, std::complex<double> *half, int half_dimx, int half_dimy, int half_dimz);
template void Prolong::prolong<std::complex<double>, 12>
  (std::complex<double> *full, std::complex<double> *half, int half_dimx, int half_dimy, int half_dimz);



Prolong::Prolong(int ratio_in, int order_in, double cmix_in, TradeImages &TR_in, Lattice &L_in, BaseGrid &BG_in) : ratio(ratio_in), order(order_in), cmix(cmix_in), TR(TR_in), L(L_in), BG(BG_in)
{
    /*Order has to be even number */
    if (order % 2)
        rmg_error_handler (__FILE__, __LINE__, "This function works only for even orders.");

    ibrav = L.get_ibrav_type();

    for (int ix = 0; ix < MAX_PROLONG_RATIO; ix++)
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
        //prolong_hex2 (full, half, dimx, dimy, dimz, half_dimx, half_dimy, half_dimz);
        if (order == 6)
        {
           prolong_hex2<T, 6>(full, half, half_dimx, half_dimy, half_dimz);
        }
        if (order == 8)
        {
           prolong_hex2<T, 8>(full, half, half_dimx, half_dimy, half_dimz);
        }
        if (order == 10)
        {
           prolong_hex2<T, 10>(full, half, half_dimx, half_dimy, half_dimz);
        }
        if (order == 12)
        {
           prolong_hex2<T, 12>(full, half, half_dimx, half_dimy, half_dimz);
        }
        return;
    }

    if(ibrav == HEXAGONAL2 && ratio == 2)
    {
//        prolong_hex2a (full, half, dimx, dimy, dimz, half_dimx, half_dimy, half_dimz);
        if (order == 6)
        {
           prolong_hex2a<T, 6>(full, half, half_dimx, half_dimy, half_dimz);
        }
        if (order == 8)
        {
           prolong_hex2a<T, 8>(full, half, half_dimx, half_dimy, half_dimz);
        }
        if (order == 10)
        {
           prolong_hex2a<T, 10>(full, half, half_dimx, half_dimy, half_dimz);
        }
        if (order == 12)
        {
           prolong_hex2a<T, 12>(full, half, half_dimx, half_dimy, half_dimz);
        }
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

    if(ratio == 2)
    {
        if (order == 6)
        {
           prolong<T, 6>(full, half, half_dimx, half_dimy, half_dimz);
        }
        if (order == 8)
        {
           prolong<T, 8>(full, half, half_dimx, half_dimy, half_dimz);
        }
        if (order == 10)
        {
           prolong<T, 10>(full, half, half_dimx, half_dimy, half_dimz);
        }
        if (order == 12)
        {
           prolong<T, 12>(full, half, half_dimx, half_dimy, half_dimz);
        }
        return;
    }

    T *sg_half = new T[(half_dimx + order) * (half_dimy + order) * (half_dimz + order)];

    TR.trade_imagesx (half, sg_half, half_dimx, half_dimy, half_dimz, order/2, FULL_TRADE);

    int incy = dimz / ratio + order;
    int incx = (dimz / ratio + order) * (dimy / ratio + order);

    int incy2 = dimz;
    int incx2 = dimz * dimy;

    int incx3 = (dimz / ratio + order) * dimy;


    // Fall through for ratio != 2
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




template <typename T, int ord> void Prolong::prolong_hex2 (T *full, T *half, int half_dimx, int half_dimy, int half_dimz)
{

    size_t n1 = half_dimx*half_dimy*half_dimz;
    size_t n2 = 8*half_dimx*half_dimy*half_dimz;

    if constexpr (std::is_same_v<T, double>)
    {
        std::vector<float> h1(n1);
        std::vector<float> f1(n2);
        for(size_t i=0;i < n1;i++) h1[i] = (float)half[i];
        if constexpr(ord == 6)
            prolong_hex_internal<float, 6, HEXAGONAL> (f1.data(), h1.data(), half_dimx, half_dimy, half_dimz);
        if constexpr(ord == 8)
            prolong_hex_internal<float, 8, HEXAGONAL> (f1.data(), h1.data(), half_dimx, half_dimy, half_dimz);
        if constexpr(ord == 10)
            prolong_hex_internal<float, 10, HEXAGONAL> (f1.data(), h1.data(), half_dimx, half_dimy, half_dimz);
        if constexpr(ord == 12)
            prolong_hex_internal<float, 12, HEXAGONAL> (f1.data(), h1.data(), half_dimx, half_dimy, half_dimz);
        for(size_t i=0;i < n2;i++) full[i] = (double)f1[i];
    }
    if constexpr (std::is_same_v<T, std::complex<double>>)
    {
        std::vector<std::complex<float>> h1(half_dimx*half_dimy*half_dimz);
        std::vector<std::complex<float>> f1(8*half_dimx*half_dimy*half_dimz);
        for(size_t i=0;i < n1;i++) h1[i] = std::complex<float>(std::real(half[i]), std::imag(half[i]));

        if constexpr(ord == 6)
            prolong_hex_internal<std::complex<float>, 6, HEXAGONAL> (f1.data(), h1.data(), half_dimx, half_dimy, half_dimz);
        if constexpr(ord == 8)
            prolong_hex_internal<std::complex<float>, 8, HEXAGONAL> (f1.data(), h1.data(), half_dimx, half_dimy, half_dimz);
        if constexpr(ord == 10)
            prolong_hex_internal<std::complex<float>, 10, HEXAGONAL> (f1.data(), h1.data(), half_dimx, half_dimy, half_dimz);
        if constexpr(ord == 12)
            prolong_hex_internal<std::complex<float>, 12, HEXAGONAL> (f1.data(), h1.data(), half_dimx, half_dimy, half_dimz);
        for(size_t i=0;i < n2;i++) full[i] = std::complex<double>(std::real(f1[i]), std::imag(f1[i]));
    }
}

template <typename T, int ord> void Prolong::prolong_hex2a (T *full, T *half, int half_dimx, int half_dimy, int half_dimz)
{

    size_t n1 = half_dimx*half_dimy*half_dimz;
    size_t n2 = 8*half_dimx*half_dimy*half_dimz;

    if constexpr (std::is_same_v<T, double>)
    {
        std::vector<float> h1(n1);
        std::vector<float> f1(n2);
        for(size_t i=0;i < n1;i++) h1[i] = (float)half[i];

        if constexpr(ord == 6)
            prolong_hex_internal<float, 6, HEXAGONAL2> (f1.data(), h1.data(), half_dimx, half_dimy, half_dimz);
        if constexpr(ord == 8)
            prolong_hex_internal<float, 8, HEXAGONAL2> (f1.data(), h1.data(), half_dimx, half_dimy, half_dimz);
        if constexpr(ord == 10)
            prolong_hex_internal<float, 10, HEXAGONAL2> (f1.data(), h1.data(), half_dimx, half_dimy, half_dimz);
        if constexpr(ord == 12)
            prolong_hex_internal<float, 12, HEXAGONAL2> (f1.data(), h1.data(), half_dimx, half_dimy, half_dimz);
        for(size_t i=0;i < n2;i++) full[i] = (double)f1[i];
    }
    if constexpr (std::is_same_v<T, std::complex<double>>)
    {
        std::vector<std::complex<float>> h1(half_dimx*half_dimy*half_dimz);
        std::vector<std::complex<float>> f1(8*half_dimx*half_dimy*half_dimz);
        for(size_t i=0;i < n1;i++) h1[i] = std::complex<float>(std::real(half[i]), std::imag(half[i]));

        if constexpr(ord == 6)
            prolong_hex_internal<std::complex<float>, 6, HEXAGONAL2> (f1.data(), h1.data(), half_dimx, half_dimy, half_dimz);
        if constexpr(ord == 8)
            prolong_hex_internal<std::complex<float>, 8, HEXAGONAL2> (f1.data(), h1.data(), half_dimx, half_dimy, half_dimz);
        if constexpr(ord == 10)
            prolong_hex_internal<std::complex<float>, 10, HEXAGONAL2> (f1.data(), h1.data(), half_dimx, half_dimy, half_dimz);
        if constexpr(ord == 12)
            prolong_hex_internal<std::complex<float>, 12, HEXAGONAL2> (f1.data(), h1.data(), half_dimx, half_dimy, half_dimz);
        for(size_t i=0;i < n2;i++) full[i] = std::complex<double>(std::real(f1[i]), std::imag(f1[i]));
    }
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

template <typename T, int ord>
void Prolong::prolong (T *full, T *half, int half_dimx, int half_dimy, int half_dimz)
{
    size_t n1 = half_dimx*half_dimy*half_dimz;
    size_t n2 = 8*half_dimx*half_dimy*half_dimz;

    if constexpr (std::is_same_v<T, double>)
    {
        std::vector<float> h1(n1);
        std::vector<float> f1(n2);
        for(size_t i=0;i < n1;i++) h1[i] = (float)half[i];

        if constexpr(ord == 6)
            prolong_internal<float, 6> (f1.data(), h1.data(), half_dimx, half_dimy, half_dimz);
        if constexpr(ord == 8)
            prolong_internal<float, 8> (f1.data(), h1.data(), half_dimx, half_dimy, half_dimz);
        if constexpr(ord == 10)
            prolong_internal<float, 10> (f1.data(), h1.data(), half_dimx, half_dimy, half_dimz);
        if constexpr(ord == 12)
            prolong_internal<float, 12> (f1.data(), h1.data(), half_dimx, half_dimy, half_dimz);
        for(size_t i=0;i < n2;i++) full[i] = (double)f1[i];
    }
    if constexpr (std::is_same_v<T, std::complex<double>>)
    {
        std::vector<std::complex<float>> h1(half_dimx*half_dimy*half_dimz);
        std::vector<std::complex<float>> f1(8*half_dimx*half_dimy*half_dimz);
        for(size_t i=0;i < n1;i++) h1[i] = std::complex<float>(std::real(half[i]), std::imag(half[i]));

        if constexpr(ord == 6)
            prolong_internal<std::complex<float>, 6> (f1.data(), h1.data(), half_dimx, half_dimy, half_dimz);
        if constexpr(ord == 8)
            prolong_internal<std::complex<float>, 8> (f1.data(), h1.data(), half_dimx, half_dimy, half_dimz);
        if constexpr(ord == 10)
            prolong_internal<std::complex<float>, 10> (f1.data(), h1.data(), half_dimx, half_dimy, half_dimz);
        if constexpr(ord == 12)
            prolong_internal<std::complex<float>, 12> (f1.data(), h1.data(), half_dimx, half_dimy, half_dimz);
        for(size_t i=0;i < n2;i++) full[i] = std::complex<double>(std::real(f1[i]), std::imag(f1[i]));
    }
}

template <typename T, int ord>
void Prolong::prolong_internal (T *full, T *half, int half_dimx, int half_dimy, int half_dimz)
{
    size_t sg_hbasis = (half_dimx + ord) * (half_dimy + ord) * (half_dimz + ord);

    RmgTimer *RT = new RmgTimer("Prolong tradeimages");
    std::vector<T> sg_half(sg_hbasis);
    TR.trade_imagesx (half, sg_half.data(), half_dimx, half_dimy, half_dimz, ord/2, FULL_TRADE);
    delete RT;

#if HIP_ENABLED
    if constexpr (std::is_same_v<T, float>)
    {
        if constexpr(ord == 6)
            prolong_ortho_gpu<float, 6> (full, sg_half.data(), half_dimx, half_dimy, half_dimz);
        if constexpr(ord == 8)
            prolong_ortho_gpu<float, 8> (full, sg_half.data(), half_dimx, half_dimy, half_dimz);
        if constexpr(ord == 10)
            prolong_ortho_gpu<float, 10> (full, sg_half.data(), half_dimx, half_dimy, half_dimz);
        if constexpr(ord == 12)
            prolong_ortho_gpu<float, 12> (full, sg_half.data(), half_dimx, half_dimy, half_dimz);
    }
    if constexpr (std::is_same_v<T, std::complex<float>>)
    {
        if constexpr(ord == 6)
            prolong_ortho_gpu<std::complex<float>, 6> (full, sg_half.data(), half_dimx, half_dimy, half_dimz);
        if constexpr(ord == 8)
            prolong_ortho_gpu<std::complex<float>, 8> (full, sg_half.data(), half_dimx, half_dimy, half_dimz);
        if constexpr(ord == 10)
            prolong_ortho_gpu<std::complex<float>, 10> (full, sg_half.data(), half_dimx, half_dimy, half_dimz);
        if constexpr(ord == 12)
            prolong_ortho_gpu<std::complex<float>, 12> (full, sg_half.data(), half_dimx, half_dimy, half_dimz);
    }
    return;
#endif

    int ic = ord/2 - 1;
    int dimx = 2 * half_dimx;
    int dimy = 2 * half_dimy;
    int dimz = 2 * half_dimz;

    int incy = half_dimz + ord;
    int incx = (half_dimy + ord) * incy;

    int incy2 = dimz;
    int incx2 = dimz * dimy;

    std::vector<T> fulla0((dimy / 2 + ord) * (dimz / 2 + ord));
    std::vector<T> fulla1((dimy / 2 + ord) * (dimz / 2 + ord));
    std::vector<T> fullb0(dimy * (dimz / 2 + ord));
    std::vector<T> fullb1(dimy * (dimz / 2 + ord));

    // lambda to clean up code
    auto stencil = [&](const T *ptr, const int stride) {
        T sum(0.0);
        if constexpr (std::is_same_v<T, double> || std::is_same_v<T, std::complex<double>>)
            for(int k = 0;k < ord;k++) sum = sum + a[1][k] * ptr[k*stride];
        if constexpr (std::is_same_v<T, float> || std::is_same_v<T, std::complex<float>>)
            for(int k = 0;k < ord;k++) sum = sum + af[1][k] * ptr[k*stride];
        return sum;
    };

    T *baseptr = sg_half.data();
    for (int ix = 0; ix < dimx / 2; ix++)
    {
        for (int iy = 0; iy < dimy / 2 + ord; iy++)
        {
            
            T *halfptr = &baseptr[(ix + 1) * incx + iy * incy];
            for (int iz = 0; iz < dimz / 2 + ord; iz++)
            {
                if constexpr (std::is_same_v<T, double> || std::is_same_v<T, std::complex<double>>)
                    fulla0[iy * incy + iz] = a[0][ic] * halfptr[ic*incx + iz];
                if constexpr (std::is_same_v<T, float> || std::is_same_v<T, std::complex<float>>)
                    fulla0[iy * incy + iz] = af[0][ic] * halfptr[ic*incx + iz];

                fulla1[iy * incy + iz] = stencil(halfptr + iz, incx);
            }
        }

        for (int iy = 0; iy < dimy / 2; iy++)
        {
            for (int iz = 0; iz < dimz / 2 + ord; iz++)
            {
                T *full_tmp = &fulla0[(iy + 1) * incy + iz];
                if constexpr (std::is_same_v<T, double> || std::is_same_v<T, std::complex<double>>)
                    fullb0[(2 * iy + 0) * incy + iz] = a[0][ic] * full_tmp[ic*incy];
                if constexpr (std::is_same_v<T, float> || std::is_same_v<T, std::complex<float>>)
                    fullb0[(2 * iy + 0) * incy + iz] = af[0][ic] * full_tmp[ic*incy];

                fullb0[(2 * iy + 1) * incy + iz] = stencil(full_tmp, incy);

                full_tmp = &fulla1[(iy + 1) * incy + iz];
                if constexpr (std::is_same_v<T, double> || std::is_same_v<T, std::complex<double>>)
                    fullb1[(2 * iy + 0) * incy + iz] = a[0][ic] * full_tmp[ic*incy];
                if constexpr (std::is_same_v<T, float> || std::is_same_v<T, std::complex<float>>)
                    fullb1[(2 * iy + 0) * incy + iz] = af[0][ic] * full_tmp[ic*incy];

                fullb1[(2 * iy + 1) * incy + iz] = stencil(full_tmp, incy);
            }
        }

        for (int iy = 0; iy < dimy; iy++)
        {
            for (int iz = 0; iz < dimz / 2; iz++)
            {
                T *full_tmp = &fullb0[iy * incy + iz + 1];
                if constexpr (std::is_same_v<T, double> || std::is_same_v<T, std::complex<double>>)
                    full[2*ix * incx2 + iy * incy2 + 2 * iz + 0] = a[0][ic] * full_tmp[ic];
                if constexpr (std::is_same_v<T, float> || std::is_same_v<T, std::complex<float>>)
                    full[2*ix * incx2 + iy * incy2 + 2 * iz + 0] = af[0][ic] * full_tmp[ic];

                full[2*ix * incx2 + iy * incy2 + 2 * iz + 1] = stencil(full_tmp, 1);

                full_tmp = &fullb1[iy * incy + iz + 1];
                if constexpr (std::is_same_v<T, double> || std::is_same_v<T, std::complex<double>>)
                    full[(2*ix + 1) * incx2 + iy * incy2 + 2 * iz + 0] = a[0][ic] * full_tmp[ic];
                if constexpr (std::is_same_v<T, float> || std::is_same_v<T, std::complex<float>>)
                    full[(2*ix + 1) * incx2 + iy * incy2 + 2 * iz + 0] = af[0][ic] * full_tmp[ic];

                full[(2*ix + 1) * incx2 + iy * incy2 + 2 * iz + 1] = stencil(full_tmp, 1);
            }
        }
    }
}


template <typename T, int ord, int htype> void Prolong::prolong_hex_internal (T *full, T *half, int half_dimx,
        int half_dimy, int half_dimz)
{
    size_t sg_hbasis = (half_dimx + ord) * (half_dimy + ord) * (half_dimz + ord);
    std::vector<T> sg_half(sg_hbasis);
    TR.trade_imagesx (half, sg_half.data(), half_dimx, half_dimy, half_dimz, ord/2, FULL_TRADE);

    int dimx = 2 * half_dimx;
    int dimy = 2 * half_dimy;
    int dimz = 2 * half_dimz;

    int incy = half_dimz + ord;
    int incx = (half_dimy + ord)*incy;

    int incy2 = half_dimz + ord;
    int incx2 = dimy * incy2;

    std::vector<T> fulla(2 * dimy * (dimz / 2 + ord));

    for (int ix = 0; ix < half_dimx; ix++)
    {
        for (int iy = 0; iy < half_dimy; iy++)
        {
            for (int iz = 0; iz < half_dimz + ord; iz++)
            {
                T sum = 0.0;
                fulla[(0) * incx2 + 2 * iy * incy2 + iz] = sg_half[(ix+ord/2) * incx + (iy+ord/2) * incy + iz];

                sum = 0.0;
                if constexpr (std::is_same_v<T, double> || std::is_same_v<T, std::complex<double>>)
                    for(int k = 0;k < ord;k++) sum+= a[1][k] * sg_half[(ix+k+1)*incx + (iy+ord/2) * incy + iz];
                if constexpr (std::is_same_v<T, float> || std::is_same_v<T, std::complex<float>>)
                    for(int k = 0;k < ord;k++) sum+= af[1][k] * sg_half[(ix+k+1)*incx + (iy+ord/2) * incy + iz];
                fulla[(1) * incx2 + 2 * iy * incy2 + iz] = sum;

                sum = 0.0;
                if constexpr (std::is_same_v<T, double> || std::is_same_v<T, std::complex<double>>)
                    for(int k = 0;k < ord;k++) sum+= a[1][k] * sg_half[(ix+ord/2)*incx + (iy+k+1) * incy + iz];
                if constexpr (std::is_same_v<T, float> || std::is_same_v<T, std::complex<float>>)
                    for(int k = 0;k < ord;k++) sum+= af[1][k] * sg_half[(ix+ord/2)*incx + (iy+k+1) * incy + iz];
                fulla[(0) * incx2 + (2 * iy +1) * incy2 + iz] = sum;

                sum = 0.0;
                // This is for hex2
                if constexpr (htype == HEXAGONAL && (std::is_same_v<T, double> || std::is_same_v<T, std::complex<double>>))
                    for(int k = 0;k < ord;k++) sum+= a[1][k] * sg_half[(ix+k+1)*incx + (iy+k+1) * incy + iz];
                if constexpr (htype == HEXAGONAL && (std::is_same_v<T, float> || std::is_same_v<T, std::complex<float>>))
                    for(int k = 0;k < ord;k++) sum+= af[1][k] * sg_half[(ix+k+1)*incx + (iy+k+1) * incy + iz];

                // This is for hex2a
                if constexpr (htype == HEXAGONAL2 && (std::is_same_v<T, double> || std::is_same_v<T, std::complex<double>>))
                    for(int k = 0;k < ord;k++) sum+= a[1][k] * sg_half[(ix+k+1)*incx + (iy + ord -k) * incy + iz];
                if constexpr (htype == HEXAGONAL2 && (std::is_same_v<T, float> || std::is_same_v<T, std::complex<float>>))
                    for(int k = 0;k < ord;k++) sum+= af[1][k] * sg_half[(ix+k+1)*incx + (iy + ord -k) * incy + iz];

                fulla[(1) * incx2 + (2 * iy +1) * incy2 + iz] = sum;
            }
        }

        for (int iy = 0; iy < dimy; iy++)
        {
            for (int iz = 0; iz < dimz / 2; iz++)
            {
                T sum = 0.0;
                T *full_tmp = &fulla[0*ix * incx2 + iy * incy + iz + 1];
                if constexpr (std::is_same_v<T, double> || std::is_same_v<T, std::complex<double>>)
                    for(int k = 0;k < ord;k++) sum+= a[0][k] * full_tmp[k];
                if constexpr (std::is_same_v<T, float> || std::is_same_v<T, std::complex<float>>)
                    for(int k = 0;k < ord;k++) sum+= af[0][k] * full_tmp[k];

                full[2*ix * dimy * dimz + iy * dimz + 2 * iz + 0] = sum;

                sum = 0.0;
                if constexpr (std::is_same_v<T, double> || std::is_same_v<T, std::complex<double>>)
                    for(int k = 0;k < ord;k++) sum+= a[1][k] * full_tmp[k];
                if constexpr (std::is_same_v<T, float> || std::is_same_v<T, std::complex<float>>)
                    for(int k = 0;k < ord;k++) sum+= af[1][k] * full_tmp[k];
                full[2*ix * dimy * dimz + iy * dimz + 2 * iz + 1] = sum;

                sum = 0.0;
                full_tmp = &fulla[(1) * incx2 + iy * incy + iz + 1];
                if constexpr (std::is_same_v<T, double> || std::is_same_v<T, std::complex<double>>)
                    for(int k = 0;k < ord;k++) sum+= a[0][k] * full_tmp[k];
                if constexpr (std::is_same_v<T, float> || std::is_same_v<T, std::complex<float>>)
                    for(int k = 0;k < ord;k++) sum+= af[0][k] * full_tmp[k];
                full[(2*ix+1) * dimy * dimz + iy * dimz + 2 * iz + 0] = sum;

                sum = 0.0;
                if constexpr (std::is_same_v<T, double> || std::is_same_v<T, std::complex<double>>)
                    for(int k = 0;k < ord;k++) sum+= a[1][k] * full_tmp[k];
                if constexpr (std::is_same_v<T, float> || std::is_same_v<T, std::complex<float>>)
                    for(int k = 0;k < ord;k++) sum+= af[1][k] * full_tmp[k];
                full[(2*ix+1) * dimy * dimz + iy * dimz + 2 * iz + 1] = sum;

            }
        }
    }
}
