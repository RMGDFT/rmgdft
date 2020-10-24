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
#include "HdfHelpers.h"

// Screening types
#define NO_SCREENING 0
#define GAU_SCREENING 1
#define ERF_SCREENING 2
#define ERFC_SCREENING 3
#define YUKAWA_SCREENING 4


template <typename T> class Exxbase {

private:

    bool gamma_extrapolation;

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
    int nstates_occ;
    int Nchol;
    // Occupations for the orbitals
    double *init_occ;

    // Base of domain distributed wavefunction array
    T *psi;

    // Exx mode
    int mode;

    // Grid points on this processing node
    int pbasis;
    int pbasis_h;

    // Data for other nodes needed to remap
    std::vector<int> dimsx;
    std::vector<int> dimsy;
    std::vector<int> dimsz;
    std::vector<int> xoffsets;
    std::vector<int> yoffsets;
    std::vector<int> zoffsets;
    std::vector<int> recvcounts;
    std::vector<int> irecvoffsets;
    std::vector<size_t> recvoffsets;


    std::vector<double> occ;

    // Mmapped serial wavefunction array
    T *psi_s;

    // File descriptor for serial wavefile
    int serial_fd;
    int exxint_fd;
    std::vector<T> ExxInt;
    std::vector<T> ExxCholVec;

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
    double yukawa = 0.0;
    double exxdiv = 0.0;


    void fftpair(T *psi_i, T*psi_j, std::complex<double> *p, double *);
    void fftpair(T *psi_i, T*psi_j, std::complex<double> *p, std::complex<float> *workbuf, double *);
    void setup_gfac(double *kq);
    void setup_exxdiv();

    std::vector< std::pair <int,int> > wf_pairs;
    int block_size = 64;
    double *kl_pair;
    double *ij_pair;
    double *Exxints;
    double *Summedints;
    std::complex<double> *wf_fft;

    double eps_qdiv = 1.0e-8;

    int VxxIntChol(std::vector<T> &Exxint, std::vector<T> &ExxCholVec, int cmax, int nstates_occ);
    void waves_pair_and_fft(int k1, int k2, std::complex<double> *Xaolj_one, std::complex<double> *Xaoik_one);
    int Vexx_int_oneQ(int iq, int_2d_array QKtoK2, std::complex<double> *Cholvec,
        std::complex<double> *phase_Gr, std::complex<double> *Xaoik, std::complex<double> *Xaolj,
        double *residual, int ij_tot, int Ncho_max, int pbasis, MPI_Comm comm);
    void write_basics(hid_t h_grp, int_2d_array QKtoK2, std::vector<int> kminus);



public:
    // BaseGrid class (distributed) and half grid

    Exxbase (
            BaseGrid &G, 
            BaseGrid &G_h, 
            Lattice &L, 
            const std::string &wavefile,
            int nstates,
            double *init_occ,
            T *psi_in, int mode_in);

    ~Exxbase(void);

    std::vector<T> Hcore;
    void Vexx(T *vexx, bool use_float_fft);
    double Exxenergy(T *vexx);
    void SetHcore(T *Hij, int lda);
    void Vexx_integrals(std::string &ifile);
    void Vexx_integrals_block(FILE *fp, int ij_start, int ij_end, int kl_start, int kl_end);
    void WriteWfsToSingleFile(void);
    void ReadWfsFromSingleFile(void);
    void kpoints_setup();
    void Remap(T *inbuf, T *rbuf);

    // Plane wave object for local grids
    Pw *pwave;

    // Plane wave object for half density grids
    Pw *pwave_h;

    double *gfac;
    double *gfac_packed;
    double *gfac_dev;
    double *gfac_dev_packed;

    std::vector<double> vexx_RMS;
};

#endif


