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
#include "Symmetry.h"

// Screening types
#define NO_SCREENING 0
#define GAU_SCREENING 1
#define ERF_SCREENING 2
#define ERFC_SCREENING 3
#define YUKAWA_SCREENING 4


template <typename T> class Exxbase {

private:
    // BaseGrid class (distributed) and half grid
    int num_q;
    double *qvec;
    double *kqvec;
    int *q_to_kindex;
    int *q_to_k_symindex;
    int *kq_index;

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
    // Number of occupied orbitals
    int nstates;
    // Occupations for the orbitals
    double *init_occ;

    // Base of domain distributed wavefunction array
    T *psi;

    // Exx mode
    int mode;

    // Grid points on this processing node
    int pbasis;
    int pbasis_h;


    std::vector<double> occ;

    // Mmapped serial wavefunction array
    T *psi_s;
    T *vexx_global;

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

    // <psi_i, psi_j> pairs that this MPI task is responsible for
    std::vector< std::pair <int,int> > pairs;

    std::mutex pair_mutex;
    double erfc_scrlen=0.0;
    double gau_scrlen=0.0;
    int scr_type = ERFC_SCREENING;


    void fftpair(T *psi_i, T*psi_j, std::complex<double> *p, double *);
    void fftpair(T *psi_i, T*psi_j, std::complex<double> *p, std::complex<float> *workbuf, double *);
    void setup_gfac(double *kq);

    std::vector< std::pair <int,int> > wf_pairs;
    int block_size = 64;
    double *kl_pair;
    double *ij_pair;
    double *Exxints;
    double *Summedints;
    std::complex<double> *wf_fft;

    double eps_qdiv = 1.0e-8;

public:
    Exxbase (
            BaseGrid &G, 
            BaseGrid &G_h, 
            Lattice &L, 
            const std::string &wavefile,
            int nstates,
            double *init_occ,
            T *psi_in, int mode_in);

    ~Exxbase(void);

    void Vexx(T *vexx, bool use_float_fft);
    double Exxenergy(T *vexx);
    void Vexx_integrals(std::string &ifile);
    void Vexx_integrals_block(FILE *fp, int ij_start, int ij_end, int kl_start, int kl_end);
    void WriteWfsToSingleFile(void);
    void ReadWfsFromSingleFile(void);
    void kpoints_setup();

    // Plane wave object for local grids
    Pw *pwave;

    // Plane wave object for half density grids
    Pw *pwave_h;

    double *gfac;

};

#endif


