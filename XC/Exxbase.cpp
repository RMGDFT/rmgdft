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

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <iterator>
#include <omp.h>


#include "const.h"
#include "Exxbase.h"
#include "RmgTimer.h"
#include "RmgException.h"
#include "RmgGemm.h"
#include "transition.h"
#include "rmgtypedefs.h"
#include "pe_control.h"
#include "GpuAlloc.h"

// This class implements exact exchange for delocalized orbitals.
// The wavefunctions are stored in a single file and are not domain
// decomposed. The file name is given by wavefile_in which is
// mmapped to an array. We only access the array in read-only mode.


template Exxbase<double>::Exxbase(BaseGrid &, BaseGrid &, Lattice &, const std::string &, int, double *, double *, int);
template Exxbase<std::complex<double>>::Exxbase(BaseGrid &, BaseGrid &, Lattice &, const std::string &, int, double *, std::complex<double> *, int );

template Exxbase<double>::~Exxbase(void);
template Exxbase<std::complex<double>>::~Exxbase(void);

template <class T> Exxbase<T>::Exxbase (
        BaseGrid &G_in,
        BaseGrid &G_h_in,
        Lattice &L_in,
        const std::string &wavefile_in,
        int nstates_in,
        double *occ_in,
        T *psi_in, int mode_in) : G(G_in), G_h(G_h_in), L(L_in), wavefile(wavefile_in), nstates(nstates_in), init_occ(occ_in), psi(psi_in), mode(mode_in)
{
    RmgTimer RT0("5-Functional: Exx init");
    for(int st=0;st < nstates;st++) occ.push_back(init_occ[st]);
    if(ct.qpoint_mesh[0] <= 0 || ct.qpoint_mesh[1] <= 0 || ct.qpoint_mesh[2] <=0)
    {
        ct.qpoint_mesh[0] = ct.kpoint_mesh[0];
        ct.qpoint_mesh[1] = ct.kpoint_mesh[1];
        ct.qpoint_mesh[2] = ct.kpoint_mesh[2];
    }
      
    tpiba = 2.0 * PI / L.celldm[0];
    tpiba2 = tpiba * tpiba;
    alpha = L.get_omega() / ((double)(G.get_NX_GRID(1) * G.get_NY_GRID(1) * G.get_NZ_GRID(1)));
    pbasis = G.get_P0_BASIS(1);

    kpoints_setup();
    if(mode == EXX_DIST_FFT) 
    {
        if(!ct.is_gamma) 
            throw RmgFatalException() << "EXX_DIST_FFT mode is only for Gamma point  \n";

        pbasis_h = G_h.get_P0_BASIS(1);
        pwave = coarse_pwaves;
        pwave_h = half_pwaves;
        psi_s = psi;
        LG = &G_in;
    }
    else
    {
        LG = new BaseGrid(G.get_NX_GRID(1), G.get_NY_GRID(1), G.get_NZ_GRID(1), 1, 1, 1, 0, 1);
        int rank = G.get_rank();
        MPI_Comm_split(G.comm, rank+1, rank, &lcomm);
        LG->set_rank(0, lcomm);
        pwave = new Pw(*LG, L, 1, false);
        size_t alloc = (size_t)pwave->pbasis * (size_t)nstates;
        vexx_global = new T[alloc]();
    }

    gfac = (double *)GpuMallocManaged((size_t)pwave->pbasis*sizeof(double));
    double kq[3] = {0.0, 0.0, 0.0};
    erfc_scrlen = Functional::get_screening_parameter_rmg();
    gau_scrlen = Functional::get_gau_parameter_rmg();
    gamma_extrapolation = true;
    setup_exxdiv();
    setup_gfac(kq);

}

template void Exxbase<double>::fftpair(double *psi_i, double *psi_j, 
        std::complex<double> *p, double *coul_fac);
template void Exxbase<std::complex<double>>::fftpair(std::complex<double> *psi_i, 
        std::complex<double> *psi_j, std::complex<double> *p,double *coul_fac);
template <class T> void Exxbase<T>::fftpair(T *psi_i, T *psi_j, std::complex<double> *p, double *coul_fac)
{
    for(size_t idx=0;idx < pwave->pbasis;idx++) p[idx] = psi_i[idx] * std::conj(psi_j[idx]);
    pwave->FftForward(p, p);
    for(size_t ig=0;ig < pwave->pbasis;ig++) p[ig] *= coul_fac[ig];
    pwave->FftInverse(p, p);
}

template void Exxbase<double>::fftpair(double *psi_i, double *psi_j, std::complex<double> *p, 
        std::complex<float> *workbuf, double *coul_fac);
template void Exxbase<std::complex<double>>::fftpair(std::complex<double> *psi_i, 
        std::complex<double> *psi_j, std::complex<double> *p,
        std::complex<float> *workbuf, double *coul_fac);
template <class T> void Exxbase<T>::fftpair(T *psi_i, T *psi_j, std::complex<double> *p, 
        std::complex<float> *workbuf, double *coul_fac)
{
    for(size_t idx=0;idx < pwave->pbasis;idx++) workbuf[idx] = std::complex<float>(psi_i[idx] * std::conj(psi_j[idx]));
    pwave->FftForward(workbuf, workbuf);
    for(size_t ig=0;ig < pwave->pbasis;ig++) workbuf[ig] *= coul_fac[ig];
    pwave->FftInverse(workbuf, workbuf);
    for(size_t idx=0;idx < pwave->pbasis;idx++) p[idx] = (std::complex<double>)workbuf[idx];
}


// This implements different ways of handling the divergence at G=0
template void Exxbase<double>::setup_gfac(double *kq);
template void Exxbase<std::complex<double>>::setup_gfac(double *kq);
template <class T> void Exxbase<T>::setup_gfac(double *kq)
{
    std::fill(gfac, gfac+pwave->pbasis, 0.0);


    double eps = 1.0e-5;
    // for HSE
    if(std::abs(erfc_scrlen) > eps) scr_type = ERFC_SCREENING;
    // for Gaupbe
    if(std::abs(gau_scrlen) > eps)
    {
         scr_type = GAU_SCREENING;
         gamma_extrapolation = false;
    }

    double a0 = 1.0;
    if(scr_type == GAU_SCREENING)
    {
        if (gau_scrlen < eps) gau_scrlen = 0.15;
        a0 = pow(PI / gau_scrlen, 1.5);
    }

    for(size_t ig=0;ig < pwave->pbasis;ig++)
    {
        double qq, v0, v1, v2;
        v0 = kq[0] + pwave->g[ig].a[0] * tpiba; 
        v1 = kq[1] + pwave->g[ig].a[1] * tpiba; 
        v2 = kq[2] + pwave->g[ig].a[2] * tpiba; 
        qq = v0* v0 + v1 * v1 + v2 * v2;
        if(!pwave->gmask[ig]) continue;
        double fac = 1.0;
        if (gamma_extrapolation)
        {
            bool on_double_grid = true;
            double qa0 = (v0 * Rmg_L.a0[0] + v1 * Rmg_L.a0[1] +v2 * Rmg_L.a0[2])/twoPI;
            double qa1 = (v0 * Rmg_L.a1[0] + v1 * Rmg_L.a1[1] +v2 * Rmg_L.a1[2])/twoPI;
            double qa2 = (v0 * Rmg_L.a2[0] + v1 * Rmg_L.a2[1] +v2 * Rmg_L.a2[2])/twoPI;
            qa0 = qa0/2.0 * ct.qpoint_mesh[0];
            qa1 = qa1/2.0 * ct.qpoint_mesh[1];
            qa2 = qa2/2.0 * ct.qpoint_mesh[2];
            if( std::abs(qa0 - std::round(qa0)) > eps )
                on_double_grid = false;
            if( std::abs(qa1 - std::round(qa1)) > eps )
                on_double_grid = false;
            if( std::abs(qa2 - std::round(qa2)) > eps )
                on_double_grid = false;

            fac = 8.0/7.0;
            if(on_double_grid) fac = 0.0;
        }
        if(scr_type == GAU_SCREENING)
        {
            gfac[ig] = a0 * exp(-qq / 4.0 / gau_scrlen);
        }
        else if (qq < eps)
        {
            gfac[ig] = -exxdiv;
        }
        else if(scr_type == ERFC_SCREENING)
        {
            
            gfac[ig] = fourPI * (1.0 - exp(-qq / 4.0 / (erfc_scrlen*erfc_scrlen))) /qq * fac;
        }
        else
        {
            gfac[ig] = fourPI/(qq + yukawa) * fac; 
        }

    }
}

// These compute the action of the exact exchange operator on all wavefunctions
// and writes the result into vfile.
template <> void Exxbase<double>::Vexx(double *vexx, bool use_float_fft)
{
    RmgTimer RT0("5-Functional: Exx potential");
    double scale = - 1.0 / (double)pwave->global_basis;

    if(!ct.is_gamma) 
        throw RmgFatalException() << "EXX_DIST_FFT mode is only for Gamma point  \n";

    // Clear vexx
    size_t stop = (size_t)nstates * (size_t)pbasis * (size_t)ct.num_kpts_pe;
    for(size_t idx=0;idx < stop;idx++) vexx[idx] = 0.0;

    int nstates_occ = 0;
    for(int st=0;st < nstates;st++) if(occ[st] > 1.0e-6) nstates_occ++;
    MPI_Allreduce(MPI_IN_PLACE, &nstates_occ, 1, MPI_INT, MPI_MAX, G.comm);

    std::vector<std::complex<double> *> pvec;
    std::vector<std::complex<float> *> wvec;
    pvec.resize(ct.OMP_THREADS_PER_NODE);
    wvec.resize(ct.OMP_THREADS_PER_NODE);
    for(int tid=0;tid < ct.OMP_THREADS_PER_NODE;tid++)
    {
        pvec[tid] = (std::complex<double> *)fftw_malloc(sizeof(std::complex<double>) * pwave->pbasis);
        wvec[tid] = (std::complex<float> *)fftw_malloc(sizeof(std::complex<float>) * pwave->pbasis);
    }

    if(mode == EXX_DIST_FFT)
    {
        // Loop over fft pairs and compute Kij(r) 
        for(int i=0;i < nstates_occ;i++)
        {
            double *psi_i = (double *)&psi_s[i*pbasis];
#pragma omp parallel for schedule(dynamic)
            for(int j=i;j < nstates;j++)
            {   
                double *psi_j = (double *)&psi_s[j*pbasis];
                int omp_tid = omp_get_thread_num();
                std::complex<double> *p = pvec[omp_tid];
                std::complex<float> *w = wvec[omp_tid];
                RmgTimer RT1("5-Functional: Exx potential fft");

                if(use_float_fft)
                    fftpair(psi_i, psi_j, p, w, gfac);
                else
                    fftpair(psi_i, psi_j, p, gfac);

#pragma omp critical(part1)
                {
                    if(j < nstates_occ)
                        for(int idx = 0;idx < pbasis;idx++)vexx[i*pbasis +idx] += scale * std::real(p[idx]) * psi_s[j*pbasis + idx];
                }
#pragma omp critical(part2)
                {
                    if(i!=j)
                        for(int idx = 0;idx < pbasis;idx++)vexx[j*pbasis +idx] += scale * std::real(p[idx]) * psi_s[i*pbasis + idx];
                }

            }
        }
    }
    else
    {
        // Write serial wavefunction files. May need to do some numa optimization here at some point
        size_t pstop = (size_t)nstates * (size_t)pwave->pbasis;
        for(size_t idx=0;idx < pstop;idx++) vexx_global[idx] = 0.0;
        RmgTimer *RT1 = new RmgTimer("5-Functional: Exx writewfs");
        WriteWfsToSingleFile();
        delete RT1;

        // Mmap wavefunction array
        std::string filename = wavefile + "_spin"+std::to_string(pct.spinpe) + "_kpt0";

        RT1 = new RmgTimer("5-Functional: mmap");
        serial_fd = open(filename.c_str(), O_RDONLY, (mode_t)0600);
        if(serial_fd < 0)
            throw RmgFatalException() << "Error! Could not open " << filename << " . Terminating.\n";
        delete RT1;

        size_t length = (size_t)nstates * (size_t)pwave->pbasis * sizeof(double);
        psi_s = (double *)mmap(NULL, length, PROT_READ, MAP_PRIVATE, serial_fd, 0);

        MPI_Barrier(G.comm);


        // Loop over blocks and process fft pairs I am responsible for
        int my_rank = G.get_rank();
        int npes = G.get_NPES();
        int rows = 0, trows=0, reqcount=0, max_rows=40;
        MPI_Request *reqs = new MPI_Request[nstates/max_rows+1];
        MPI_Status *mrstatus = new MPI_Status[nstates/max_rows+1];

        int flag=0;
        double *arptr = vexx_global;

        for(int i=0;i < nstates_occ;i++)
        {
            double *psi_i = (double *)&psi_s[(size_t)i*(size_t)pwave->pbasis];
            RmgTimer *RT1 = new RmgTimer("5-Functional: Exx potential fft");
#pragma omp parallel for schedule(dynamic)
            for(int j=i;j < nstates;j++)
            {
                int omp_tid = omp_get_thread_num();
                std::complex<double> *p = pvec[omp_tid];
                std::complex<float> *w = wvec[omp_tid];
                int fft_index = i*nstates + j;
                if(my_rank == (fft_index % npes))
                {
                    double *psi_j = (double *)&psi_s[(size_t)j*(size_t)pwave->pbasis];
                    if(use_float_fft)
                        fftpair(psi_i, psi_j, p, w, gfac);
                    else
                        fftpair(psi_i, psi_j, p, gfac);
                    // We can speed this up by adding more critical sections if it proves to be a bottleneck
#pragma omp critical(part3)
                    {
                        if(j < nstates_occ)
                            for(size_t idx = 0;idx < (size_t)pwave->pbasis;idx++) 
                                vexx_global[(size_t)i*(size_t)pwave->pbasis + idx] += scale * std::real(p[idx]) * psi_j[idx];
                    }
#pragma omp critical(part4)
                    {
                        if(i!=j)
                            for(size_t idx = 0;idx < (size_t)pwave->pbasis;idx++) 
                                vexx_global[(size_t)j*(size_t)pwave->pbasis +idx] += scale * std::real(p[idx]) * psi_i[idx];
                    }
                }
            }
            delete RT1;

            rows++;
            // This enables some concurrency with calculations and communication
            if(rows == max_rows)
            {
                MPI_Barrier(G.comm);
                MPI_Iallreduce(MPI_IN_PLACE, arptr, rows*pwave->pbasis, MPI_DOUBLE, MPI_SUM, G.comm, &reqs[reqcount]);
                reqcount++;
                arptr += rows*pwave->pbasis;
                trows += rows;
                rows = 0;
            }
            MPI_Testall(reqcount, reqs, &flag, mrstatus);
        }
        MPI_Waitall(reqcount, reqs, mrstatus);
        MPI_Barrier(G.comm);
        int count = nstates - trows;

        // This timer only picks up the last MPI_Allreduce since the ones above are asynchronous
        RmgTimer *RT3 = new RmgTimer("5-Functional: Exx allreduce");
        MPI_Allreduce(MPI_IN_PLACE, arptr, count*pwave->pbasis, MPI_DOUBLE, MPI_SUM, G.comm);
        delete RT3;

        delete [] mrstatus;
        delete [] reqs;

        // Map my portion of vexx into 
        int xoffset, yoffset, zoffset;
        int dimx = G.get_PX0_GRID(1);
        int dimy = G.get_PY0_GRID(1);
        int dimz = G.get_PZ0_GRID(1);
        int gdimy = G.get_NY_GRID(1);
        int gdimz = G.get_NZ_GRID(1);
        G.find_node_offsets(G.get_rank(), G.get_NX_GRID(1), G.get_NY_GRID(1), G.get_NZ_GRID(1), &xoffset, &yoffset, &zoffset);

        RmgTimer *RT2 = new RmgTimer("5-Functional: Exx remap");
        for(size_t i=0;i < (size_t)nstates;i++)
        {
            for(size_t ix=0;ix < (size_t)dimx;ix++)
            {
                for(size_t iy=0;iy < (size_t)dimy;iy++)
                {
                    double *tvexx = &vexx[i*(size_t)pbasis];
                    double *tvexx_global = &vexx_global[i*(size_t)pwave->pbasis];
                    for(size_t iz=0;iz < (size_t)dimz;iz++)
                    {
                        tvexx[ix*(size_t)dimy*(size_t)dimz + iy*(size_t)dimz + iz] = 
                        tvexx_global[(ix+(size_t)xoffset)*(size_t)gdimy*(size_t)gdimz + (iy+(size_t)yoffset)*(size_t)gdimz + iz + (size_t)zoffset];
                    }
                }
            }
        }
        delete RT2;

        MPI_Barrier(G.comm);

        // munmap wavefunction array
        munmap(psi_s, length);
        close(serial_fd);

    }

    for(int tid=0;tid < ct.OMP_THREADS_PER_NODE;tid++)
    {
        fftw_free(wvec[tid]);
        fftw_free(pvec[tid]);
    }

}


template <> double Exxbase<double>::Exxenergy(double *vexx)
{
    double energy = 0.0;
    for(int st=0;st < nstates;st++)
    {
        double scale = 0.5 * occ[st];
        for(int i=0;i < pbasis;i++) energy += scale*vexx[st*pbasis + i]*psi[st*pbasis + i];
    }
    energy = ct.exx_fraction*energy * this->L.get_omega() / (double)this->G.get_GLOBAL_BASIS(1);
    MPI_Allreduce(MPI_IN_PLACE, &energy, 1, MPI_DOUBLE, MPI_SUM, this->G.comm);
    MPI_Allreduce(MPI_IN_PLACE, &energy, 1, MPI_DOUBLE, MPI_SUM, pct.spin_comm);

    return energy;
}

template <> double Exxbase<std::complex<double>>::Exxenergy(std::complex<double> *vexx)
{
    double energy = 0.0;
    for(int ik = 0; ik < ct.num_kpts_pe; ik++)
    {
        int ik_glob = ik + pct.kstart;
        for(int st=0;st < nstates;st++)
        {
            double scale = 0.5 * occ[st] * ct.kp[ik_glob].kweight;
            for(int i=0;i < pbasis;i++) 
                energy += scale*std::real(vexx[ik * ct.run_states * pbasis + st*pbasis + i]*std::conj(psi[ik * ct.max_states * pbasis + st*pbasis + i]));
        }
    }
    energy = ct.exx_fraction*energy * this->L.get_omega() / (double)this->G.get_GLOBAL_BASIS(1);
    MPI_Allreduce(MPI_IN_PLACE, &energy, 1, MPI_DOUBLE, MPI_SUM, this->G.comm);
    MPI_Allreduce(MPI_IN_PLACE, &energy, 1, MPI_DOUBLE, MPI_SUM, pct.spin_comm);
    MPI_Allreduce(MPI_IN_PLACE, &energy, 1, MPI_DOUBLE, MPI_SUM, pct.kpsub_comm);

    return energy;
}

template <> void Exxbase<std::complex<double>>::Vexx_integrals_block(FILE *fp, int ij_start, int ij_end, int kl_start, int kl_end)
{
}
template <> void Exxbase<double>::Vexx_integrals_block(FILE *fp,  int ij_start, int ij_end, int kl_start, int kl_end)
{

    double beta = 0.0;

    char *trans_a = "t";
    char *trans_n = "n";

    RmgTimer *RT0 = new RmgTimer("5-Functional: Exx: kl");

    // calculate kl pairs
    for(int kl = kl_start; kl < kl_end; kl++)
    {
        int ie = kl- kl_start;
        int k = wf_pairs[kl].first;
        int l = wf_pairs[kl].second;
        double *psi_k = (double *)&psi_s[k*pbasis];
        double *psi_l = (double *)&psi_s[l*pbasis];
        for(int idx=0;idx < pbasis;idx++) kl_pair[ie*pbasis + idx] = psi_k[idx]*psi_l[idx];
    }

    delete RT0;
    // Now compute integrals for (i, j, k, l)

    int ij_length = ij_end - ij_start;
    int kl_length = kl_end - kl_start;
    // Now matrix multiply to produce a block of (1,jblocks, 1, nstates_occ) results
    RT0 = new RmgTimer("5-Functional: Exx: gemm");
    RmgGemm(trans_a, trans_n, ij_length, kl_length, pbasis, alpha, ij_pair, pbasis, kl_pair, pbasis, beta, Exxints, ij_length);
    delete RT0;
    int pairsize = ij_length * kl_length;
    RT0 = new RmgTimer("5-Functional: Exx: reduce");
    MPI_Reduce(Exxints, Summedints, pairsize, MPI_DOUBLE, MPI_SUM, 0, LG->comm);
    MPI_Barrier(LG->comm);
    delete RT0;

    RT0 = new RmgTimer("5-Functional: Exx: print");
    for(int ij = ij_start; ij < ij_end; ij++)
    {
        int ic = ij - ij_start;
        int i = wf_pairs[ij].first;
        int j = wf_pairs[ij].second;
        for(int kl = kl_start; kl < kl_end; kl++)
        {
            int ie = kl- kl_start;
            int k = wf_pairs[kl].first;
            int l = wf_pairs[kl].second;

            if(LG->get_rank()==0 && kl >= ij)fprintf(fp, "%14.8e %d %d %d %d\n", Summedints[ie * ij_length + ic], i, j, k, l);
        }
    }

    delete RT0;
}

// This computes exact exchange integrals
// and writes the result into vfile.
template <> void Exxbase<double>::Vexx_integrals(std::string &vfile)
{

    double scale = 1.0 / (double)pwave->global_basis;
    int nstates_occ = 0;
    for(int st=0;st < nstates;st++) if(occ[st] > 1.0e-6) nstates_occ++;

    RmgTimer RT0("5-Functional: Exx integrals");

    for(int i=0;i < nstates_occ;i++)
    {
        for(int j=i;j < nstates_occ;j++)
        {
            wf_pairs.push_back(std::make_pair(i, j));
        }
    }

    //block_size = 16;
    int num_blocks = (wf_pairs.size() + block_size -1)/block_size;


    // The wave function array is always at least double the number of states so use
    // the upper part for storage of our kl pairs

    kl_pair = new double[block_size * pbasis];
    ij_pair = new double[block_size * pbasis];

    Exxints = new double[block_size * block_size];
    Summedints = new double[block_size * block_size];

    wf_fft = new std::complex<double> [pbasis];


    char *buf=NULL;
    FILE *fp=NULL;
    if(LG->get_rank()==0)
    {
        int rank = G.get_rank();
        size_t size = 1<<18;
        buf = new char[size];
        std::string filename = vfile + "_spin"+std::to_string(pct.spinpe)+".fcidump_"+std::to_string(rank);
        fp= fopen(filename.c_str(), "w");
        if(int nn = setvbuf (fp, buf, _IOFBF, size) != 0) printf("\n failed buff %d\n", nn);

        if(rank == 0)
        {
            fprintf(fp, "&FCI\n");
            fprintf(fp, "NORB=%d,\n", nstates_occ);
            fprintf(fp, "NELEC=%d,\n", (int)(ct.nel+1.0e-6));
            int ms2 = pct.spinpe;
            if(ct.spin_flag && pct.spinpe == 0) ms2 = -1;
            fprintf(fp, "MS2=%d,\n", ms2);
            fprintf(fp, "\n&END\n");
        }
    }

    std::vector<int> blocks_in_thispe;

    if(mode == EXX_DIST_FFT)
    {
        for(int ib =0; ib < num_blocks; ib++)
            blocks_in_thispe.push_back(ib);
    }
    else
    {
        std::string filename = wavefile + "_spin"+std::to_string(pct.spinpe);
        serial_fd = open(filename.c_str(), O_RDONLY, (mode_t)0600);
        if(serial_fd < 0)
            throw RmgFatalException() << "Error! Could not open " << filename << " . Terminating.\n";

        size_t length = nstates * pbasis *sizeof(double);
        psi_s = (double *)mmap(NULL, length, PROT_READ, MAP_PRIVATE, serial_fd, 0);

        MPI_Barrier(G.comm);

        int npes = G.get_NPES();
        int rank = G.get_rank();
        //   first npes blocks assigned to processor 0,1,2,... npes-1
        //  second npes blocks assigned to processor npes-1, npes-2, ... 2, 1, 0
        // thus wuill balance the blocks on each processor
        for(int ib =0; ib < num_blocks; ib++)
        {
            int ibb = ib%(2*npes);
            if(ibb == rank) 
            {
                blocks_in_thispe.push_back(ib);
            }
            else if( (ibb >= npes) && (npes - 1 - ibb%npes == rank) )
            {

                blocks_in_thispe.push_back(ib);
            }

        }

    }

    for(std::vector<int>::iterator it = blocks_in_thispe.begin(); it != blocks_in_thispe.end(); it++)
    {
        int ib = *it;
        int ij_start = ib * block_size;
        int ij_end = (ib+1) * block_size;
        if ((size_t)ij_end > wf_pairs.size()) ij_end = wf_pairs.size();
        RmgTimer *RT1 = new RmgTimer("5-Functional: Exx: ij+fft");
        for(int ij = ij_start; ij < ij_end; ij++)
        {

            int ic = ij - ij_start;
            int i = wf_pairs[ij].first;
            int j = wf_pairs[ij].second;
            double *psi_i = (double *)&psi_s[i*pbasis];
            double *psi_j = (double *)&psi_s[j*pbasis];
            fftpair(psi_i, psi_j, wf_fft, gfac);

            // store (i,j) fft pair in ij_pair
            // forward and then backward fft should not have the scale, Wenchang
            for(int idx=0;idx < pbasis;idx++) ij_pair[ic * pbasis + idx] = scale * std::real(wf_fft[idx]);

            double tem1 =0.0, tem2 = 0.0, tem3 = 0.0;
            for(int idx=0;idx < pbasis;idx++) 
            {
                tem1 += ij_pair[ic * pbasis + idx];
                tem2 += psi_i[ic * pbasis + idx];
                tem3 += psi_j[ic * pbasis + idx];
            }
        }
        delete RT1;

        for(int ic = ib; ic < num_blocks; ic++)
        {

            int kl_start = ic * block_size;
            int kl_end = (ic+1) * block_size;
            if ((size_t)kl_end > wf_pairs.size()) kl_end = wf_pairs.size();
            Vexx_integrals_block(fp, ij_start, ij_end, kl_start, kl_end);
        }
    }

    if(LG->get_rank()==0)
    {
        fclose(fp);
        delete []buf;
    }


    delete [] ij_pair;
    delete [] kl_pair;
    delete [] Exxints;
    delete [] Summedints;
    delete [] wf_fft;

}

template <> void Exxbase<std::complex<double>>::Vexx_integrals(std::string &vfile)
{
    printf("Exx mode not programmed yet\n");
}

template <class T> Exxbase<T>::~Exxbase(void)
{

    if(mode == EXX_DIST_FFT) return;

    close(serial_fd);
    size_t length = nstates * pbasis * sizeof(T);
    munmap(psi_s, length);
    //std::string filename= wavefile + "_spin"+ std::to_string(pct.spinpe);
    //unlink(filename.c_str());

    GpuFreeManaged(gfac);
    delete LG;
    MPI_Comm_free(&lcomm); 
    delete pwave;
    delete [] vexx_global;
}

template void Exxbase<double>::ReadWfsFromSingleFile();
template void Exxbase<std::complex<double>>::ReadWfsFromSingleFile();
template <class T> void Exxbase<T>::ReadWfsFromSingleFile()
{
    // Write the domain distributed wavefunction array and map it to psi_s
    MPI_Datatype wftype = MPI_DOUBLE;
    if(typeid(T) == typeid(std::complex<double>)) wftype = MPI_DOUBLE_COMPLEX;

    int sizes_c[3];
    int subsizes_c[3];
    int starts_c[3];

    sizes_c[0] = G.get_NX_GRID(1);
    sizes_c[1] = G.get_NY_GRID(1);
    sizes_c[2] = G.get_NZ_GRID(1);

    subsizes_c[0] = G.get_PX0_GRID(1);
    subsizes_c[1] = G.get_PY0_GRID(1);
    subsizes_c[2] = G.get_PZ0_GRID(1);

    starts_c[0] = G.get_PX_OFFSET(1);
    starts_c[1] = G.get_PY_OFFSET(1);
    starts_c[2] = G.get_PZ_OFFSET(1);

    int order = MPI_ORDER_C;
    MPI_Info fileinfo;
    MPI_Datatype grid_c;
    MPI_Status status;

    MPI_Type_create_subarray(3, sizes_c, subsizes_c, starts_c, order, wftype, &grid_c);
    MPI_Type_commit(&grid_c);

    MPI_Info_create(&fileinfo);

    int amode = MPI_MODE_RDONLY;
    MPI_File mpi_fhand ;

    MPI_Barrier(G.comm);

    std::string filename = wavefile + "_spin"+std::to_string(pct.spinpe);
    int flag = MPI_File_open(G.comm, filename.c_str(), amode, fileinfo, &mpi_fhand);
    if(flag) 
        throw RmgFatalException() << "Error! Could not open " << filename << " . Terminating.\n";
    MPI_Offset disp = 0;

    T *wfptr;
    MPI_File_set_view(mpi_fhand, disp, wftype, grid_c, "native", MPI_INFO_NULL);
    int dis_dim = G.get_P0_BASIS(1);
    for(int st=0;st < nstates;st++)
    {
        wfptr = &psi[st * dis_dim];
        MPI_File_read_all(mpi_fhand, wfptr, dis_dim, MPI_DOUBLE, &status);
    }
    MPI_Barrier(G.comm);
    MPI_File_close(&mpi_fhand);
    MPI_Type_free(&grid_c);
    fflush(NULL);
    MPI_Barrier(G.comm);
}

template void Exxbase<double>::WriteWfsToSingleFile();
template void Exxbase<std::complex<double>>::WriteWfsToSingleFile();
template <class T> void Exxbase<T>::WriteWfsToSingleFile()
{
    // Write the domain distributed wavefunction array and map it to psi_s
    MPI_Datatype wftype = MPI_DOUBLE;
    if(typeid(T) == typeid(std::complex<double>)) wftype = MPI_DOUBLE_COMPLEX;

    int sizes_c[3];
    int subsizes_c[3];
    int starts_c[3];

    sizes_c[0] = G.get_NX_GRID(1);
    sizes_c[1] = G.get_NY_GRID(1);
    sizes_c[2] = G.get_NZ_GRID(1);

    subsizes_c[0] = G.get_PX0_GRID(1);
    subsizes_c[1] = G.get_PY0_GRID(1);
    subsizes_c[2] = G.get_PZ0_GRID(1);

    starts_c[0] = G.get_PX_OFFSET(1);
    starts_c[1] = G.get_PY_OFFSET(1);
    starts_c[2] = G.get_PZ_OFFSET(1);

    int order = MPI_ORDER_C;
    MPI_Info fileinfo;
    MPI_Datatype grid_c;
    MPI_Status status;

    MPI_Type_create_subarray(3, sizes_c, subsizes_c, starts_c, order, wftype, &grid_c);
    MPI_Type_commit(&grid_c);

    MPI_Info_create(&fileinfo);

    int amode = MPI_MODE_RDWR|MPI_MODE_CREATE;
    MPI_File mpi_fhand ;

    MPI_Barrier(G.comm);

    for(int ik = 0; ik < ct.num_kpts_pe; ik++)
    {
        int ik_glob = ik + pct.kstart;
        std::string filename = wavefile + "_spin"+std::to_string(pct.spinpe) + "_kpt" + std::to_string(ik_glob);
        MPI_File_open(G.comm, filename.c_str(), amode, fileinfo, &mpi_fhand);
        MPI_Offset disp = 0;

        T *wfptr;
        MPI_File_set_view(mpi_fhand, disp, wftype, grid_c, "native", MPI_INFO_NULL);
        int dis_dim = G.get_P0_BASIS(1);
        for(int st=0;st < nstates;st++)
        {
            wfptr = &psi[ik * ct.max_states * dis_dim + st * dis_dim];
            MPI_File_write_all(mpi_fhand, wfptr, dis_dim, wftype, &status);
        }
        MPI_Barrier(G.comm);
        MPI_File_close(&mpi_fhand);
    }
    MPI_Type_free(&grid_c);
    fflush(NULL);
    MPI_Barrier(G.comm);
}


template void Exxbase<double>::kpoints_setup();
template void Exxbase<std::complex<double>>::kpoints_setup();
template <class T> void Exxbase<T>::kpoints_setup()
{
    num_q = ct.qpoint_mesh[0] * ct.qpoint_mesh[1] * ct.qpoint_mesh[2];
    qvec = new double[3* num_q];
    kqvec = new double[3* num_q];
    q_to_kindex = new int[num_q];
    q_to_k_symindex = new int[num_q];
    kq_index = new int[ct.num_kpts * num_q];


    if(ct.is_gamma) 
    {
        qvec[0] = 0.0;
        qvec[1] = 0.0;
        qvec[2] = 0.0;
        q_to_kindex[0] = 0;
        q_to_k_symindex[0] = 0;

        kqvec[0] = 0.0;
        kqvec[1] = 0.0;
        kqvec[2] = 0.0;
        kq_index[0] = 0;

        return;
    }

    if(!Rmg_Symm)
        throw RmgFatalException() << " Symmetry not defined for non-gamma point in exx \n";

    // get all kpoints in the first BZ
    int num_q_temp = ct.num_kpts;
    for(int ik = 0; ik < ct.num_kpts; ik++)
    {
        qvec[ik * 3 + 0] = ct.kp[ik].kpt[0];
        qvec[ik * 3 + 1] = ct.kp[ik].kpt[1];
        qvec[ik * 3 + 2] = ct.kp[ik].kpt[2];
        q_to_kindex[ik] = ik;
        q_to_k_symindex[ik] = 0;
    }

    double sym_qvec[3], dk[3];
    for(int ik = 0; ik < ct.num_kpts; ik++)
        for(int isym = 0; isym < Rmg_Symm->nsym; isym++)
        {
            sym_qvec[0] = Rmg_Symm->sym_rotate[isym * 9 + 0 * 3 + 0 ] * ct.kp[ik].kpt[0] +
                Rmg_Symm->sym_rotate[isym * 9 + 0 * 3 + 1 ] * ct.kp[ik].kpt[1] +
                Rmg_Symm->sym_rotate[isym * 9 + 0 * 3 + 2 ] * ct.kp[ik].kpt[2];
            sym_qvec[1] = Rmg_Symm->sym_rotate[isym * 9 + 1 * 3 + 0 ] * ct.kp[ik].kpt[0] +
                Rmg_Symm->sym_rotate[isym * 9 + 1 * 3 + 1 ] * ct.kp[ik].kpt[1] +
                Rmg_Symm->sym_rotate[isym * 9 + 1 * 3 + 2 ] * ct.kp[ik].kpt[2];
            sym_qvec[2] = Rmg_Symm->sym_rotate[isym * 9 + 2 * 3 + 0 ] * ct.kp[ik].kpt[0] +
                Rmg_Symm->sym_rotate[isym * 9 + 2 * 3 + 1 ] * ct.kp[ik].kpt[1] +
                Rmg_Symm->sym_rotate[isym * 9 + 2 * 3 + 2 ] * ct.kp[ik].kpt[2];

            //  check if this new qvec in the qvec list
            bool is_assigned = false;
            for(int iq = 0; iq < num_q_temp; iq++)
            {
                dk[0] = sym_qvec[0] - qvec[iq * 3 + 0];
                dk[1] = sym_qvec[1] - qvec[iq * 3 + 1];
                dk[2] = sym_qvec[2] - qvec[iq * 3 + 2];
                dk[0] = dk[0] - std::round(dk[0]);
                dk[1] = dk[1] - std::round(dk[1]);
                dk[2] = dk[2] - std::round(dk[2]);
                if( (std::abs(dk[0]) < 1.0e-10) && (std::abs(dk[1]) < 1.0e-10) && (std::abs(dk[2]) < 1.0e-10) )
                {
                    is_assigned = true;
                    break;
                }
            }

            if(!is_assigned)
            {
                if(num_q_temp + 1 > num_q)
                {
                    throw RmgFatalException() << num_q_temp << " num_q larger than mesh in exx " << num_q <<"\n";
                }
                qvec[num_q_temp * 3 + 0] = sym_qvec[0] ;
                qvec[num_q_temp * 3 + 1] = sym_qvec[1] ;
                qvec[num_q_temp * 3 + 2] = sym_qvec[2] ;
                q_to_kindex[num_q_temp] = ik;
                q_to_k_symindex[num_q_temp] = isym;
                num_q_temp++;
            }

            //  check if this new -qvec in the qvec list
            is_assigned = false;
            for(int iq = 0; iq < num_q_temp; iq++)
            {
                dk[0] = -sym_qvec[0] - qvec[iq * 3 + 0];
                dk[1] = -sym_qvec[1] - qvec[iq * 3 + 1];
                dk[2] = -sym_qvec[2] - qvec[iq * 3 + 2];
                dk[0] = dk[0] - std::round(dk[0]);
                dk[1] = dk[1] - std::round(dk[1]);
                dk[2] = dk[2] - std::round(dk[2]);
                if( (std::abs(dk[0]) < 1.0e-10) && (std::abs(dk[1]) < 1.0e-10) && (std::abs(dk[2]) < 1.0e-10) )
                {
                    is_assigned = true;
                    break;
                }
            }

            if(!is_assigned)
            {
                if(num_q_temp + 1 > num_q)
                    throw RmgFatalException() << num_q_temp << " num_q larger than mesh in -exx " << num_q <<"\n";
                qvec[num_q_temp * 3 + 0] = -sym_qvec[0];
                qvec[num_q_temp * 3 + 1] = -sym_qvec[1];
                qvec[num_q_temp * 3 + 2] = -sym_qvec[2];
                q_to_kindex[num_q_temp] = ik;
                q_to_k_symindex[num_q_temp] = -isym;
                num_q_temp++;
            }


        }

    if(ct.verbose && pct.imgpe == 0)
    {
        for(int iq = 0; iq < num_q_temp; iq++)
            printf("\n qvec %f %f %f %d %d", qvec[iq*3], qvec[iq*3+1], qvec[iq*3+2], q_to_kindex[iq], q_to_k_symindex[iq]);
    }

    if(num_q_temp != num_q)
        throw RmgFatalException() << num_q_temp << " num_q wrong in exx " << num_q <<"\n";

    for(int iq = 0; iq < num_q; iq++)
    {
        double v1, v2, v3;

        v1 = qvec[iq * 3 + 0] *Rmg_L.b0[0]
            +qvec[iq * 3 + 1] *Rmg_L.b1[0] 
            +qvec[iq * 3 + 2] *Rmg_L.b2[0];
        v2 = qvec[iq * 3 + 0] *Rmg_L.b0[1]
            +qvec[iq * 3 + 1] *Rmg_L.b1[1] 
            +qvec[iq * 3 + 2] *Rmg_L.b2[1];
        v3 = qvec[iq * 3 + 0] *Rmg_L.b0[2]
            +qvec[iq * 3 + 1] *Rmg_L.b1[2] 
            +qvec[iq * 3 + 2] *Rmg_L.b2[2];

        qvec[iq * 3 + 0] = v1 * twoPI;
        qvec[iq * 3 + 1] = v2 * twoPI;
        qvec[iq * 3 + 2] = v3 * twoPI;
    }
}
template <> void Exxbase<std::complex<double>>::Vexx(std::complex<double> *vexx, bool use_float_fft)
{
    RmgTimer RT0("5-Functional: Exx potential");
    double scale = - 1.0 / (double)pwave->global_basis;

    if(mode == EXX_DIST_FFT)
        throw RmgFatalException() << "EXX_DIST_FFT mode is only for Gamma point  \n";

    int nstates_occ = 0;
    for(int st=0;st < nstates;st++) if(std::abs(occ[st]) > 1.0e-6) nstates_occ++;
    MPI_Allreduce(MPI_IN_PLACE, &nstates_occ, 1, MPI_INT, MPI_MAX, G.comm);

    std::vector<std::complex<double> *> pvec;
    std::vector<std::complex<float> *> wvec;
    pvec.resize(ct.OMP_THREADS_PER_NODE);
    wvec.resize(ct.OMP_THREADS_PER_NODE);
    for(int tid=0;tid < ct.OMP_THREADS_PER_NODE;tid++)
    {
        pvec[tid] = (std::complex<double> *)fftw_malloc(sizeof(std::complex<double>) * pwave->pbasis);
        wvec[tid] = (std::complex<float> *)fftw_malloc(sizeof(std::complex<float>) * pwave->pbasis);
    }

    // Write serial wavefunction files. May need to do some numa optimization here at some point
    RmgTimer *RT1 = new RmgTimer("5-Functional: Exx writewfs");
    WriteWfsToSingleFile();
    delete RT1;

    // Loop over blocks and process fft pairs I am responsible for
    int my_rank = G.get_rank();
    int npes = G.get_NPES();

    size_t length = (size_t)nstates * (size_t)pwave->pbasis * sizeof(std::complex<double>);
    std::complex<double> *psi_q = (std::complex<double> *)GpuMallocManaged(length);
    std::complex<double> *psi_q_map;
    double kq[3];
    for(int ik = 0; ik < ct.num_kpts_pe; ik++)
    {
        int ik_glob = ik + pct.kstart;
        // Mmap wavefunction array
        std::string filename = wavefile + "_spin"+std::to_string(pct.spinpe) + "_kpt" + std::to_string(ik_glob);
        RT1 = new RmgTimer("5-Functional: mmap");
        serial_fd = open(filename.c_str(), O_RDONLY, (mode_t)0600);
        if(serial_fd < 0)
            throw RmgFatalException() << "Error! Could not open " << filename << " . Terminating.\n";
        delete RT1;

        psi_s = (std::complex<double> *)mmap(NULL, length, PROT_READ, MAP_PRIVATE, serial_fd, 0);

        MPI_Barrier(G.comm);

        for(size_t idx=0;idx < (size_t)nstates*pwave->pbasis;idx++) vexx_global[idx] = 0.0;
        for(int iq = 0; iq < num_q; iq++)
        {

            int ikindex = q_to_kindex[iq];
            int isym = q_to_k_symindex[iq];
            int isyma = std::abs(isym);


            kq[0] = ct.kp[ik_glob].kvec[0] - qvec[iq * 3 +0];
            kq[1] = ct.kp[ik_glob].kvec[1] - qvec[iq * 3 +1];
            kq[2] = ct.kp[ik_glob].kvec[2] - qvec[iq * 3 +2];

            setup_gfac(kq);

            std::string filename_q = wavefile + "_spin"+std::to_string(pct.spinpe) + "_kpt" + std::to_string(ikindex);
            RT1 = new RmgTimer("5-Functional: mmap");
            int serial_fd_q = open(filename_q.c_str(), O_RDONLY, (mode_t)0600);
            if(serial_fd_q < 0)
                throw RmgFatalException() << "Error! Could not open " << filename << " . Terminating.\n";
            delete RT1;

            psi_q_map = (std::complex<double> *)mmap(NULL, length, PROT_READ, MAP_PRIVATE, serial_fd_q, 0);

            MPI_Barrier(G.comm);


            // rotate wavefunctions for q point from symmetry-related k point.
            // define exp(-i (k-q) r) 
            int nx_grid = G.get_NX_GRID(1);
            int ny_grid = G.get_NY_GRID(1);
            int nz_grid = G.get_NZ_GRID(1);

            int nbasis = nx_grid * ny_grid * nz_grid;
            int ixx, iyy, izz;
            for (int ix = 0; ix < nx_grid; ix++) {
                for (int iy = 0; iy < ny_grid; iy++) {
                    for (int iz = 0; iz < nz_grid; iz++) {

                        symm_ijk(&Rmg_Symm->s[isyma *9], &Rmg_Symm->ftau_wave[isyma*3], ix, iy, iz, ixx, iyy, izz, nx_grid, ny_grid, nz_grid);

                        if(isym >= 0)
                        {
                            for(int st = 0; st < nstates_occ; st++)
                            {
                                psi_q[st * nbasis + ixx * ny_grid * nz_grid + iyy * nz_grid + izz]
                                    = (psi_q_map[st * nbasis + ix * ny_grid * nz_grid + iy * nz_grid + iz]);
                            }
                        }
                        else
                        {
                            for(int st = 0; st < nstates_occ; st++)
                            {
                                psi_q[st * nbasis + ixx * ny_grid * nz_grid + iyy * nz_grid + izz]
                                    = std::conj(psi_q_map[st * nbasis + ix * ny_grid * nz_grid + iy * nz_grid + iz]);
                            }
                        }
                    }
                }
            }


            for(int i=0;i < nstates;i++)
            {
                std::complex<double> *psi_i = &psi_s[i*pwave->pbasis];
                RmgTimer *RT1 = new RmgTimer("5-Functional: Exx potential fft");
#pragma omp parallel for schedule(dynamic)
                for(int j=0;j < nstates_occ;j++)
                {
                    int omp_tid = omp_get_thread_num();
                    std::complex<double> *p = pvec[omp_tid];
                    std::complex<float> *w = wvec[omp_tid];
                    int fft_index = i*nstates + j;
                    if(my_rank == (fft_index % npes))
                    {
                        std::complex<double> *psi_j = &psi_q[j*pwave->pbasis];
                        if(use_float_fft)
                            fftpair(psi_i, psi_j, p, w, gfac);
                        else
                            fftpair(psi_i, psi_j, p, gfac);
                        // We can speed this up by adding more critical sections if it proves to be a bottleneck

#pragma omp critical(part3)
                        {
                            for(size_t idx = 0;idx < pwave->pbasis;idx++) 
                                vexx_global[i*pwave->pbasis +idx] += scale * p[idx] * psi_j[idx];
                        }
                    }
                }
                delete RT1;
            }

            munmap(psi_q_map, length);
            close(serial_fd_q);
        }

        MPI_Barrier(G.comm);
        RmgTimer *RT3 = new RmgTimer("5-Functional: Exx allreduce");
        MPI_Allreduce(MPI_IN_PLACE, vexx_global, 2*nstates*pwave->pbasis, MPI_DOUBLE, MPI_SUM, G.comm);
        delete RT3;

        // Map my portion of vexx into 
        int xoffset, yoffset, zoffset;
        int dimx = G.get_PX0_GRID(1);
        int dimy = G.get_PY0_GRID(1);
        int dimz = G.get_PZ0_GRID(1);
        int gdimy = G.get_NY_GRID(1);
        int gdimz = G.get_NZ_GRID(1);
        G.find_node_offsets(G.get_rank(), G.get_NX_GRID(1), G.get_NY_GRID(1), G.get_NZ_GRID(1), &xoffset, &yoffset, &zoffset);

        RmgTimer *RT2 = new RmgTimer("5-Functional: Exx remap");
        for(int i=0;i < nstates;i++)
        {
            std::complex<double> *tvexx = &vexx[ik * ct.run_states * pbasis + i*pbasis];
            std::complex<double> *tvexx_global = &vexx_global[i*pwave->pbasis];
            for(int ix=0;ix < dimx;ix++)
            {
                for(int iy=0;iy < dimy;iy++)
                {
                    for(int iz=0;iz < dimz;iz++)
                    {
                        tvexx[ix*dimy*dimz + iy*dimz + iz] = 
                            tvexx_global[(ix+xoffset)*gdimy*gdimz + (iy+yoffset)*gdimz + iz + zoffset]/(double)num_q;
                    }
                }
            }
        }
        delete RT2;

        MPI_Barrier(G.comm);

        // munmap wavefunction array
        munmap(psi_s, length);
        close(serial_fd);

    }

    GpuFreeManaged(psi_q);
    for(int tid=0;tid < ct.OMP_THREADS_PER_NODE;tid++)
    {
        fftw_free(wvec[tid]);
        fftw_free(pvec[tid]);
    }

}

template void Exxbase<double>::setup_exxdiv();
template void Exxbase<std::complex<double>>::setup_exxdiv();
template <class T> void Exxbase<T>::setup_exxdiv()
{
    double dqx = 1.0/(double)ct.qpoint_mesh[0]; 
    double dqy = 1.0/(double)ct.qpoint_mesh[1]; 
    double dqz = 1.0/(double)ct.qpoint_mesh[2]; 
    double xq[3];
    double eps = 1.0e-5;
    double grid_factor = 8.0/7.0;
    bool on_double_grid = false;
    double alpha = 10.0/(coarse_pwaves->gcut * ct.filter_factor * tpiba2);
//    alpha = 0.833333333333;
    exxdiv = 0.0;
    for(int iqx = 0; iqx < ct.qpoint_mesh[0];iqx++)
        for(int iqy = 0; iqy < ct.qpoint_mesh[1];iqy++)
            for(int iqz = 0; iqz < ct.qpoint_mesh[2];iqz++)
            {
                xq[0] = Rmg_L.b0[0] * iqx * dqx + 
                    Rmg_L.b1[0] * iqy * dqy + 
                    Rmg_L.b2[0] * iqz * dqz;
                xq[1] = Rmg_L.b0[1] * iqx * dqx + 
                    Rmg_L.b1[1] * iqy * dqy + 
                    Rmg_L.b2[1] * iqz * dqz;
                xq[2] = Rmg_L.b0[2] * iqx * dqx + 
                    Rmg_L.b1[2] * iqy * dqy + 
                    Rmg_L.b2[2] * iqz * dqz;

                for(size_t ig=0;ig < coarse_pwaves->pbasis;ig++)
                {
                    double qq, v0, v1, v2;
                    v0 = xq[0]*twoPI + coarse_pwaves->g[ig].a[0] * tpiba;
                    v1 = xq[1]*twoPI + coarse_pwaves->g[ig].a[1] * tpiba;
                    v2 = xq[2]*twoPI + coarse_pwaves->g[ig].a[2] * tpiba;
                    qq = v0*v0 + v1 * v1 + v2 * v2;
                    if (gamma_extrapolation)
                    {
                        on_double_grid = true;
                        double qa0 = (v0 * Rmg_L.a0[0] + v1 * Rmg_L.a0[1] +v2 * Rmg_L.a0[2])/twoPI;
                        double qa1 = (v0 * Rmg_L.a1[0] + v1 * Rmg_L.a1[1] +v2 * Rmg_L.a1[2])/twoPI;
                        double qa2 = (v0 * Rmg_L.a2[0] + v1 * Rmg_L.a2[1] +v2 * Rmg_L.a2[2])/twoPI;
                        qa0 = qa0/2.0 * ct.qpoint_mesh[0];
                        qa1 = qa1/2.0 * ct.qpoint_mesh[1];
                        qa2 = qa2/2.0 * ct.qpoint_mesh[2];
                        if( std::abs(qa0 - std::round(qa0)) > eps )
                            on_double_grid = false;
                        if( std::abs(qa1 - std::round(qa1)) > eps )
                            on_double_grid = false;
                        if( std::abs(qa2 - std::round(qa2)) > eps )
                            on_double_grid = false;
                    }
                
                    if(!on_double_grid && qq > eps)
                    {
                        if(scr_type == ERFC_SCREENING)
                        {
                            exxdiv += grid_factor * std::exp(-alpha * qq)/qq * 
                                (1.0 - std::exp(-qq /4.0/erfc_scrlen/erfc_scrlen));
                        }
                        else
                        {
                            exxdiv += grid_factor * std::exp(-alpha * qq)/(qq + yukawa);
                        }

                    }
                }


            }
    MPI_Allreduce(MPI_IN_PLACE, &exxdiv, 1, MPI_DOUBLE, MPI_SUM, this->G.comm);

    if(!gamma_extrapolation)
    {
        if(scr_type == ERFC_SCREENING)
            exxdiv += 1.0/4.0/erfc_scrlen/erfc_scrlen;
        else if( std::abs(yukawa) > eps )
            exxdiv += 1.0/yukawa;
        else
            exxdiv += -alpha;
    }
    exxdiv *= fourPI /num_q;

    

    int nqq = 100000;
    double dq = 5.0/std::sqrt(alpha) /nqq;
    double aa = 0.0;
    for(int iq = 0; iq < nqq; iq++)
    {
        double qq = dq *(iq+0.5) ;
        qq = qq * qq;
        if(scr_type == ERFC_SCREENING)
        {
            aa -= std::exp(-alpha * qq) * std::exp(-qq/4.0/erfc_scrlen/erfc_scrlen) * dq;
        }
        else
        {
            aa -= std::exp(-alpha * qq) * yukawa /(qq+yukawa) * dq;
        }


    }
    aa = aa * 2.0/PI;
    aa += 1.0/std::sqrt(alpha * PI);
    exxdiv -= aa * L.get_omega();

    exxdiv *= num_q;
    if(ct.verbose && pct.gridpe == 0) 
    {
        printf("\n exxdiv = %f %f %f", exxdiv, aa, alpha);
        printf("\n erfc_scrlen = %f", erfc_scrlen);
    }
}
