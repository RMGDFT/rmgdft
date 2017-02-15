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

#ifndef RMG_Atomic_H
#define RMG_Atomic_H 1


#include <complex>
#include "BaseGrid.h"
#include "Lattice.h"
#include "TradeImages.h"
#include "FiniteDiff.h"

#define RADIAL_GVECS 1000


#ifdef __cplusplus

void InitLocalObject (double *sumobject, double * &lobject, int object_type, bool compute_lobject);
void InitDelocalizedObject (double *sumobject, double * &lobject, int object_type, bool compute_lobject);
void InitLocalObject(double *sumobject, int object_type);
void LcaoGetAtomicRho(double *arho);

class Atomic {

private:

    double Gcutoff (double g1, double gcut, double width);

    static int Log_grid_initialized;
    static double a;
    static double b;

    // BaseGrid class
    BaseGrid *Grid;

    // Lattice object
    Lattice *L;

public:
    Atomic (void);
    ~Atomic (void);
    void RftToLogGrid (double cparm, double * f, double * r, double * ffil, double *rab, int rg_points, int lval, double width);
    void InitBessel(double *r, int rg_points, int lmax, double *bessel_rg);
    void RLogGridToGLogGrid (double * f, double * r, double *rab, double * ffil, int rg_points, int lval, double *bessel_rg);
    double BesselToLogGrid (double cparm, double * f, double * r, double * ffil, double *rab, int rg_points, int lval, double rcut, double hmin);
    void FilterPotential (double *potential, double *r, int rg_points, double rmax,
                              double parm, double* potential_lgrid, double *rab,
                              int l_value, double gwidth, double rcut, double rwidth, double hmin);
    double *GetRgrid(void);
    double Interpolate(double *f, double r);
    double GetRange(double *f, double *r, double *rab, int rg_points);


    static double r_filtered[MAX_LOGGRID];
    static double log_r_filtered[MAX_LOGGRID];
    static double gvec[RADIAL_GVECS];
    void PackFine2Rhogrid(std::complex<double> *gweight_fine, int ngrid_fine, std::complex<double> *gweigh, int ngrid);

};

#endif
#endif
