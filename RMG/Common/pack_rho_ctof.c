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

#include "portability.h"
#include <stdio.h>
#include <stdlib.h>
#include "grid.h"
#include "common_prototypes.h"
#include "main.h"
#include <float.h>
#include "TradeImages.h"

void pack_rho_ctof (double * rho, double * rho_f)
{
    int i, j, k, ii, jj, kk, basis1, basis2, basis3;
    int pbasis, dimx, dimy, dimz, in, jn, kn;
    int ifxs, ifys;
    int icxs, icys, alloc;
    int PX0_GRID, PY0_GRID, PZ0_GRID, FG_RATIO;
    double tmp1, tmp2, tmp3, frac, cc[10][4];
    double *rho_c, sum_rho, sum_rhof, coef;

    PX0_GRID = get_PX0_GRID();
    PY0_GRID = get_PY0_GRID();
    PZ0_GRID = get_PZ0_GRID();
    FG_RATIO = get_FG_RATIO();

    alloc = (PX0_GRID + 4) * (PY0_GRID + 4) * (PZ0_GRID + 4);
    my_malloc(rho_c, alloc, double);

    ifxs = get_FPY0_GRID() * get_FPZ0_GRID();
    ifys = get_FPZ0_GRID();

    icxs = (PY0_GRID + 4) * (PZ0_GRID + 4);
    icys = PZ0_GRID + 4;

    int num = 0;


    dimx = PX0_GRID;
    dimy = PY0_GRID;
    dimz = PZ0_GRID;
    pbasis =get_P0_BASIS();

    sum_rho = 0.0;
    sum_rhof = 0.0;
    for (i = 0; i < pbasis; i++)
        sum_rho += rho[i];
    sum_rho = real_sum_all (sum_rho, pct.grid_comm);
    sum_rho *= get_vel();

    for (i = 0; i < FG_RATIO; i++)
    {
        frac = (double) i / (double) FG_RATIO;
        cc[i][0] = -frac * (1.0 - frac) * (2.0 - frac) / 6.0;
        cc[i][1] = (1.0 + frac) * (1.0 - frac) * (2.0 - frac) / 2.0;
        cc[i][2] = (1.0 + frac) * frac * (2.0 - frac) / 2.0;
        cc[i][3] = -(1.0 + frac) * frac * (1.0 - frac) / 6.0;
    }

    trade_imagesx (rho, rho_c, dimx, dimy, dimz, 2, FULL_TRADE);

    for (i = 2; i < PX0_GRID + 2; i++)
    {
        for (j = 2; j < PY0_GRID + 2; j++)
        {
            for (k = 2; k < PZ0_GRID + 2; k++)
            {
                rho_f[(FG_RATIO * (i - 2))*ifxs + (FG_RATIO * (j - 2))*ifys + FG_RATIO * (k - 2)] =
                    rho_c[i*icxs + j*icys + k];
                ++num;
            }
        }
    }

    for (i = 2; i < PX0_GRID + 2; i++)
    {
        for (j = 2; j < PY0_GRID + 2; j++)
        {
            for (k = 2; k < PZ0_GRID + 2; k++)
            {

                for (in = 1; in < FG_RATIO; in++)
                {
                    tmp1 = 0.0;
                    tmp2 = 0.0;
                    tmp3 = 0.0;
                    basis1 = -1;

                    for (ii = 0; ii < 4; ii++)
                    {
                        tmp1 += cc[in][ii] * rho_c[(i + basis1)*icxs + j*icys + k];
                        tmp2 += cc[in][ii] * rho_c[i*icxs + (j + basis1)*icys + k];
                        tmp3 += cc[in][ii] * rho_c[i*icxs + j*icys + k + basis1];
                        ++basis1;
                    }

                    rho_f[(FG_RATIO * (i - 2) + in)*ifxs + (FG_RATIO * (j - 2))*ifys + FG_RATIO * (k - 2)] = tmp1;
                    rho_f[(FG_RATIO * (i - 2))*ifxs + (FG_RATIO * (j - 2) + in)*ifys + FG_RATIO * (k - 2)] = tmp2;
                    rho_f[(FG_RATIO * (i - 2))*ifxs + (FG_RATIO * (j - 2))*ifys + FG_RATIO * (k - 2) + in] = tmp3;
                    num += 3;
                }

            }
        }
    }

    for (i = 2; i < PX0_GRID + 2; i++)
    {
        for (j = 2; j < PY0_GRID + 2; j++)
        {
            for (k = 2; k < PZ0_GRID + 2; k++)
            {

                for (in = 1; in < FG_RATIO; in++)
                {
                    for (jn = 1; jn < FG_RATIO; jn++)
                    {

                        tmp1 = 0.0;
                        tmp2 = 0.0;
                        tmp3 = 0.0;
                        basis1 = -1;
                        for (ii = 0; ii < 4; ii++)
                        {
                            basis2 = -1;
                            for (jj = 0; jj < 4; jj++)
                            {
                                tmp1 +=
                                    cc[in][ii] * cc[jn][jj] * rho_c[(i + basis1)*icxs +(j + basis2)*icys + k];
                                tmp2 +=
                                    cc[in][ii] * cc[jn][jj] * rho_c[(i + basis1)*icxs + j*icys + k + basis2];
                                tmp3 +=
                                    cc[in][ii] * cc[jn][jj] * rho_c[i*icxs + (j + basis1)*icys + k + basis2];
                                ++basis2;
                            }
                            ++basis1;
                        }
                        rho_f[(FG_RATIO * (i - 2) + in)*ifxs + (FG_RATIO * (j - 2) + jn)*ifys + FG_RATIO * (k - 2)] = tmp1;
                        rho_f[(FG_RATIO * (i - 2) + in)*ifxs + (FG_RATIO * (j - 2))*ifys + FG_RATIO * (k - 2) + jn] = tmp2;
                        rho_f[(FG_RATIO * (i - 2))*ifxs + (FG_RATIO * (j - 2) + in)*ifys + FG_RATIO * (k - 2) + jn] = tmp3;
                        num += 3;
                    }
                }
            }
        }
    }


    for (i = 2; i < PX0_GRID + 2; i++)
    {
        for (j = 2; j < PY0_GRID + 2; j++)
        {
            for (k = 2; k < PZ0_GRID + 2; k++)
            {

                for (in = 1; in < FG_RATIO; in++)
                {
                    for (jn = 1; jn < FG_RATIO; jn++)
                    {
                        for (kn = 1; kn < FG_RATIO; kn++)
                        {

                            tmp1 = 0.0;
                            basis1 = -1;
                            for (ii = 0; ii < 4; ii++)
                            {
                                basis2 = -1;
                                for (jj = 0; jj < 4; jj++)
                                {
                                    basis3 = -1;
                                    for (kk = 0; kk < 4; kk++)
                                    {
                                        tmp1 +=
                                            cc[in][ii] * cc[jn][jj] * cc[kn][kk] * rho_c[(i + basis1)*icxs + (j + basis2)*icys + k + basis3];
                                        ++basis3;
                                    }
                                    ++basis2;
                                }
                                ++basis1;
                            }
                            rho_f[(FG_RATIO * (i - 2) + in)*ifxs + (FG_RATIO * (j - 2) + jn)*ifys + FG_RATIO * (k - 2) + kn] = tmp1;
                            num += 1;

                        }
                    }
                }
            }
        }
    }

    if (num != get_FP0_BASIS())
        error_handler ("there is something wrong here");

    for (i = 0; i < num; i++)
        sum_rhof += rho_f[i];
    sum_rhof = real_sum_all (sum_rhof, pct.grid_comm);
    sum_rhof *= get_vel_f();
    coef = sum_rho / sum_rhof;
    for (i = 0; i < num; i++)
        rho_f[i] *= coef;

    my_free(rho_c);

}
