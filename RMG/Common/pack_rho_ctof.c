/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include <float.h>

void pack_rho_ctof (REAL * rho, REAL * rho_f)
{
    int i, j, k, ii, jj, kk, basis1, basis2, basis3;
    int pbasis, dimx, dimy, dimz, in, jn, kn;
    int ifxs, ifys;
    int icxs, icys, alloc;
    REAL tmp1, tmp2, tmp3, frac, cc[10][4];
    REAL *rho_c, sum_rho, sum_rhof, coef;

    alloc = (pct.PX0_GRID + 4) * (pct.PY0_GRID + 4) * (pct.PZ0_GRID + 4);
    my_malloc(rho_c, alloc, REAL);

    ifxs = pct.FPY0_GRID * pct.FPZ0_GRID;
    ifys = pct.FPZ0_GRID;

    icxs = (pct.PY0_GRID + 4) * (pct.PZ0_GRID + 4);
    icys = pct.PZ0_GRID + 4;

    int num = 0;

    dimx = pct.PX0_GRID;
    dimy = pct.PY0_GRID;
    dimz = pct.PZ0_GRID;
    pbasis =pct.P0_BASIS;

    sum_rho = 0.0;
    sum_rhof = 0.0;
    for (i = 0; i < pbasis; i++)
        sum_rho += rho[i];
    sum_rho = real_sum_all (sum_rho, pct.grid_comm);
    sum_rho *= ct.vel;

    for (i = 0; i < FG_NX; i++)
    {
        frac = (REAL) i / (REAL) FG_NX;
        cc[i][0] = -frac * (1.0 - frac) * (2.0 - frac) / 6.0;
        cc[i][1] = (1.0 + frac) * (1.0 - frac) * (2.0 - frac) / 2.0;
        cc[i][2] = (1.0 + frac) * frac * (2.0 - frac) / 2.0;
        cc[i][3] = -(1.0 + frac) * frac * (1.0 - frac) / 6.0;
    }

    trade_imagesx (rho, rho_c, dimx, dimy, dimz, 2, FULL_FD);

    for (i = 2; i < pct.PX0_GRID + 2; i++)
    {
        for (j = 2; j < pct.PY0_GRID + 2; j++)
        {
            for (k = 2; k < pct.PZ0_GRID + 2; k++)
            {
                rho_f[(FG_NX * (i - 2))*ifxs + (FG_NX * (j - 2))*ifys + FG_NX * (k - 2)] =
                    rho_c[i*icxs + j*icys + k];
                ++num;
            }
        }
    }

    for (i = 2; i < pct.PX0_GRID + 2; i++)
    {
        for (j = 2; j < pct.PY0_GRID + 2; j++)
        {
            for (k = 2; k < pct.PZ0_GRID + 2; k++)
            {

                for (in = 1; in < FG_NX; in++)
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

                    rho_f[(FG_NX * (i - 2) + in)*ifxs + (FG_NX * (j - 2))*ifys + FG_NX * (k - 2)] = tmp1;
                    rho_f[(FG_NX * (i - 2))*ifxs + (FG_NX * (j - 2) + in)*ifys + FG_NX * (k - 2)] = tmp2;
                    rho_f[(FG_NX * (i - 2))*ifxs + (FG_NX * (j - 2))*ifys + FG_NX * (k - 2) + in] = tmp3;
                    num += 3;
                }

            }
        }
    }

    for (i = 2; i < pct.PX0_GRID + 2; i++)
    {
        for (j = 2; j < pct.PY0_GRID + 2; j++)
        {
            for (k = 2; k < pct.PZ0_GRID + 2; k++)
            {

                for (in = 1; in < FG_NX; in++)
                {
                    for (jn = 1; jn < FG_NX; jn++)
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
                        rho_f[(FG_NX * (i - 2) + in)*ifxs + (FG_NX * (j - 2) + jn)*ifys + FG_NX * (k - 2)] = tmp1;
                        rho_f[(FG_NX * (i - 2) + in)*ifxs + (FG_NX * (j - 2))*ifys + FG_NX * (k - 2) + jn] = tmp2;
                        rho_f[(FG_NX * (i - 2))*ifxs + (FG_NX * (j - 2) + in)*ifys + FG_NX * (k - 2) + jn] = tmp3;
                        num += 3;
                    }
                }
            }
        }
    }


    for (i = 2; i < pct.PX0_GRID + 2; i++)
    {
        for (j = 2; j < pct.PY0_GRID + 2; j++)
        {
            for (k = 2; k < pct.PZ0_GRID + 2; k++)
            {

                for (in = 1; in < FG_NX; in++)
                {
                    for (jn = 1; jn < FG_NX; jn++)
                    {
                        for (kn = 1; kn < FG_NX; kn++)
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
                            rho_f[(FG_NX * (i - 2) + in)*ifxs + (FG_NX * (j - 2) + jn)*ifys + FG_NX * (k - 2) + kn] = tmp1;
                            num += 1;

                        }
                    }
                }
            }
        }
    }

    if (num != pct.FP0_BASIS)
        error_handler ("there is something wrong here");

    for (i = 0; i < num; i++)
        sum_rhof += rho_f[i];
    sum_rhof = real_sum_all (sum_rhof, pct.grid_comm);
    sum_rhof *= ct.vel_f;
    coef = sum_rho / sum_rhof;
    for (i = 0; i < num; i++)
        rho_f[i] *= coef;

    my_free(rho_c);

}
