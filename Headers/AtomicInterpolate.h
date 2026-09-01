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

#ifndef RMG_AtomicInterpolate_H
#define RMG_AtomicInterpolate_H 1

#include "params.h"
#include "math.h"

extern double Atomic_inv_a;
extern double Atomic_inv_b;
extern double GlogGrid_inv_a;
extern double GlogGrid_inv_b;

// Interpolates function f that is defined on the shared logarithmic grid
static inline double AtomicInterpolateInline(double *f, double r)
{
    double d0, d1, dm;
    double f0, g0, g1, g2, h1, h2, i2;
    int ic;

    // truncate to nearest integer index into r_filtered array
    if((r < LOGGRID_START)) {
        r = LOGGRID_START;
        if(fabs(f[0]) < 1.0e-5) return 0.0;
    }

    d0 = (log (r*Atomic_inv_a) * Atomic_inv_b);
    ic = (int)d0;
    ic = (ic > 0) ? ic : 1;

    /* cubic interpolation using forward differences */
    d0 -= (double) (ic);
    d1 = (d0 - 1.0) * 0.5;
    dm = (d0 - 2.0) / 3.0;

    f0 = f[ic];
    g0 = f[ic] - f[ic - 1];
    g1 = f[ic + 1] - f[ic];
    g2 = f[ic + 2] - f[ic + 1];
    h1 = g1 - g0;
    h2 = g2 - g1;
    i2 = h2 - h1;

    return f0 + d0 * (g1 + d1 * (h2 + dm * i2));

}

static inline double AtomicInterpolateInline_Ggrid(double *f, double g)
{
    double d0, d1, dm;
    double f0, g0, g1, g2, h1, h2, i2;
    int ic;
//   first Gvec = 0, then start log grid.

    if((g < LOGGRID_START)) {
        return f[0];
    }

    d0 = (log (g*GlogGrid_inv_a) * GlogGrid_inv_b);
    ic = (int)d0;
    ic = (ic > 0) ? ic : 1;

    /* cubic interpolation using forward differences */
    d0 -= (double) (ic);
    d1 = (d0 - 1.0) * 0.5;
    dm = (d0 - 2.0) / 3.0;

    f0 = f[ic];
    g0 = f[ic] - f[ic - 1];
    g1 = f[ic + 1] - f[ic];
    g2 = f[ic + 2] - f[ic + 1];
    h1 = g1 - g0;
    h2 = g2 - g1;
    i2 = h2 - h1;

    return f0 + d0 * (g1 + d1 * (h2 + dm * i2));

}

static inline double AtomicInterpolateInline_1(double *f, int ic, double d0, double d1, double dm)
{
    double f0, g0, g1, g2, h1, h2, i2;

    f0 = f[ic];
    g0 = f[ic] - f[ic - 1];
    g1 = f[ic + 1] - f[ic];
    g2 = f[ic + 2] - f[ic + 1];
    h1 = g1 - g0;
    h2 = g2 - g1;
    i2 = h2 - h1;

    return f0 + d0 * (g1 + d1 * (h2 + dm * i2));

}

#endif
