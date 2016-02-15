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

#ifndef RMG_vdW_H
#define RMG_vdW_H 1


#include "BaseGrid.h"
#include "Lattice.h"
#include "TradeImages.h"
#include "FiniteDiff.h"

/* Maximum number of points to use in the log interpolation grid */
#define         MAX_LOGGRID      (1400)

/* Starting radius for the log interpolation grid */
#define         LOGGRID_START    (0.00001)

/* Log interpolation grid mesh parameter */
#define         LOGGRID_MESHPARM  (1.01)

#define SMALL 1.0e-8

#define RADIAL_GVECS 1500


#ifdef __cplusplus

class Atomic {

private:

    double Gcutoff (double g1, double gcut, double width);

    static int Log_grid_initialized;

    // BaseGrid class
    BaseGrid *Grid;

    // Lattice object
    Lattice *L;

public:
    Atomic (void);
    ~Atomic (void);
    void RftToLogGrid (double cparm, double * f, double * r, double * ffil, double *rab, int rg_points, int lval, double width);
    void FilterPotential (double *potential, double *r, int rg_points, double rmax,
                              double offset, double parm, double* potential_lgrid, double *rab,
                              int l_value, double gwidth, double rcut, double rwidth, double *drpotential_lgrid);

    static double r_filtered[MAX_LOGGRID];
    static double log_r_filtered[MAX_LOGGRID];

};

#endif
#endif
