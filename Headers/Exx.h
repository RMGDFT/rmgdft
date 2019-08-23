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


#ifndef RMG_Exx_H
#define RMG_Exx_H 1


#include <string>
#include <complex>
#include "BaseGrid.h"
#include "Lattice.h"
#include "fftw3.h"

template <typename T> class Exx {

private:
    // BaseGrid class
    BaseGrid &G;

    // Lattice object
    Lattice &L;

    // File path for wavefunction file. Spin and kpoint identifiers should be added by parent.
    std::string &wavefile;

    // Grid points on this processing node
    int pbasis;

    // Total number of grid points
    size_t N;

    // Number of orbitals
    int nstates;

    // Base of domain distributed wavefunction array
    T *psi;

    // Mmapped serial wavefunction array
    T *psi_s;

    // File descriptor for serial wavefile
    int serial_fd;

    // Each MPI process keeps a portion of the orbitals resident in memory and
    // these two class members control that.
    int st_start;
    int st_count;

public:
    Exx (BaseGrid &G, 
         Lattice &L, 
         std::string &wavefile,
         int nstates,
         T *psi_in);

    ~Exx(void);

};

#endif
