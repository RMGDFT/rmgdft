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


#ifndef RMG_Pdos_H
#define RMG_Pdos_H 1


#include <string>
#include <vector>
#include <set>
#include <complex>
#include <mutex>
#include "BaseGrid.h"
#include "Lattice.h"
#include "Symmetry.h"
#include "Exxbase.h"
#include "Kpoint.h"


template <typename T> class Pdos {

private:
    // BaseGrid class (distributed) and half grid

    Exxbase<T> *Exx;
    BaseGrid &G;

    // Lattice object
    Lattice &L;

    // File path for wavefunction file. Spin and kpoint identifiers should be added by parent.
    const std::string wavefile;
    // Number of states
    int nstates;

    // Base of domain distributed wavefunction array
    T *psi;


    // Grid points on this processing node
    int ngrid, ngrid_noncoll;

    // Mmapped serial wavefunction array
    T *psi_s;

    // File descriptor for serial wavefile
    int serial_fd;

    int nx_grid, ny_grid, nz_grid;
    int px0_grid, py0_grid, pz0_grid;
    int px_offset, py_offset, pz_offset;
    int pbasis, pbasis_noncoll;

public:
    Pdos (
            BaseGrid &G, 
            Lattice &L, 
            const std::string wavefile,
            int nstates,
            T *psi_in, Kpoint<T> **Kptr);

    ~Pdos(void);

    void ReadRotatePsi(int ik, int isy, int isya, std::string wavefile, T *psi_k);
    void Pdos_calc(Kpoint<T> **Kptr, std::vector<double> eigs);

};
#endif
