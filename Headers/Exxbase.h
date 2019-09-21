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


#ifndef RMG_Exxbase_H
#define RMG_Exxbase_H 1


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

// Screening types
#define NO_SCREENING 0
#define GAU_SCREENING 1
#define ERF_SCREENING 2
#define ERFC_SCREENING 3


template <typename T> class Exxbase {

private:
    // BaseGrid class (distributed)
    BaseGrid &G;

    // Lattice object
    Lattice &L;

    // Don't need to keep recomputing these
    double tpiba;
    double tpiba2;
    double alpha;

    // File path for wavefunction file. Spin and kpoint identifiers should be added by parent.
    const std::string &wavefile;

    // Exx mode
    int mode;

    // Grid points on this processing node
    int pbasis;

    // Total number of grid points
    size_t N;

    // Number of occupied orbitals
    int nstates;

    // Occupations for the orbitals
    double *occ;

    // Base of domain distributed wavefunction array
    T *psi;

    // Mmapped serial wavefunction array
    T *psi_s;

    // File descriptor for serial wavefile
    int serial_fd;

    // Each MPI process keeps a portion of the orbitals resident in memory and
    // these two class members control that.
    int pair_start;
    int pair_count;

    // Local MPI communicator
    MPI_Comm lcomm;

    // BaseGrid instance for local grids
    BaseGrid *LG;

    // Plane wave object for local grids
    Pw *pwave;

    // <psi_i, psi_j> pairs that this MPI task is responsible for
    std::vector< std::pair <int,int> > pairs;

    std::mutex pair_mutex;
    int scr_type =1;

    double *gfac;

    void fftpair(T *psi_i, T*psi_j, std::complex<double> *p);
    void setup_gfac(void);

    std::vector< std::pair <int,int> > wf_pairs;
    int block_size = 64;
    double *kl_pair;
    double *ij_pair;
    double *Exxints;
    double *Summedints;
    std::complex<double> *wf_fft;


public:
    Exxbase (
            BaseGrid &G, 
            Lattice &L, 
            const std::string &wavefile,
            int nstates,
            double *occ,
            T *psi_in);

    ~Exxbase(void);

    void Vexx(T *vexx);
    void Vexx_integrals(std::string &ifile);
    void Vexx_integrals_block(FILE *fp, int ij_start, int ij_end, int kl_start, int kl_end);
    void WriteWfsToSingleFile(void);

};

#endif


