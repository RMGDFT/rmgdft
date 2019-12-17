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


#ifndef ON_Exx_on_H
#define ON_Exx_on_H 1


#include <string>
#include <vector>
#include <set>
#include <complex>
#include <mutex>
#include "BaseGrid.h"
#include "Lattice.h"
#include "fftw3.h"
#include "Pw.h"
#include "Functional.h"
#include "Exxbase.h"
#include "LocalObject.h"



template <typename T> class Exx_on {

private:
    // BaseGrid class (distributed) and half grid
    BaseGrid &G;
    BaseGrid &G_h;

    // Lattice object
    Lattice &L;

    // Don't need to keep recomputing these
    double tpiba;
    double tpiba2;
    double alpha;

    // File path for wavefunction file. Spin and kpoint identifiers should be added by parent.
    const std::string &wavefile;
    LocalObject<T> &Phi;

    double *occ;
    // Exx mode
    int mode;

    // Grid points on this processing node
    int pbasis;
    int pbasis_h;

    // Number of  orbitals
    int nstates;
    Exxbase<T> *Exxb;
    //  Dephi_k = sum_l  K_kl * Phi_l 
    T *DePhi;

public:
    Exx_on (
            BaseGrid &G, 
            BaseGrid &G_h, 
            Lattice &L, 
            const std::string &wavefile,
            LocalObject<T> &Phi, double *occ, int mode_in);

    ~Exx_on(void);
    void Omega(T *rho_matrix, bool use_fft_float);
    void Omega_rmg(T *Cij_local, T *Cij_glob, bool use_fft_float);
    void Xij(T *Sij_inverse, LocalObject<T> &Phi);
    double Exxenergy(T *rho_matrix);
    void HijExx(T *Hij_glob, LocalObject<T> &Phi);
    void OmegaSinv(T *Sij_inverse, LocalObject<T> &Phi);
    void OmegaRes(T *res, LocalObject<T> &Phi);
    T *Xij_mat;
    T *Xij_mat_Sinv;
    T *Omega_j;
    T *PreOrbital;

};

#endif

extern Exx_on<double> *Exx_onscf;

