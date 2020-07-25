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

#include "WriteEshdf.h"


#include "const.h"
#include "Exxbase.h"
#include "RmgTimer.h"
#include "RmgException.h"
#include "RmgGemm.h"
#include "transition.h"
#include "rmgtypedefs.h"
#include "pe_control.h"
#include "GpuAlloc.h"
#include "blas.h"

// This class implements exact exchange for delocalized orbitals.
// The wavefunctions are stored in a single file and are not domain
// decomposed. The file name is given by wavefile_in which is
// mmapped to an array. We only access the array in read-only mode.

using namespace std;
//using namespace hdfHelper;

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
    }

    gfac = (double *)GpuMallocManaged((size_t)pwave->pbasis*sizeof(double));
    double kq[3] = {0.0, 0.0, 0.0};
    erfc_scrlen = Functional::get_screening_parameter_rmg();
    gau_scrlen = Functional::get_gau_parameter_rmg();
    setup_exxdiv();
    setup_gfac(kq);

    int npes = G.get_NPES();
    int my_rank = G.get_rank();
    dimsx.resize(npes, 0);
    dimsy.resize(npes, 0);
    dimsz.resize(npes, 0);

    // Need to know these on other nodes
    dimsx[my_rank] = G.get_PX0_GRID(1);
    dimsy[my_rank] = G.get_PY0_GRID(1);
    dimsz[my_rank] = G.get_PZ0_GRID(1);
    MPI_Allreduce(MPI_IN_PLACE, dimsx.data(), npes, MPI_INT, MPI_SUM, G.comm);
    MPI_Allreduce(MPI_IN_PLACE, dimsy.data(), npes, MPI_INT, MPI_SUM, G.comm);
    MPI_Allreduce(MPI_IN_PLACE, dimsz.data(), npes, MPI_INT, MPI_SUM, G.comm);

    int xoffset, yoffset, zoffset;
    xoffsets.resize(npes, 0);
    yoffsets.resize(npes, 0);
    zoffsets.resize(npes, 0);
    G.find_node_offsets(G.get_rank(), G.get_NX_GRID(1), G.get_NY_GRID(1), G.get_NZ_GRID(1), &xoffset, &yoffset, &zoffset);
    xoffsets[my_rank] = xoffset;
    yoffsets[my_rank] = yoffset;
    zoffsets[my_rank] = zoffset;
    MPI_Allreduce(MPI_IN_PLACE, xoffsets.data(), npes, MPI_INT, MPI_SUM, G.comm);
    MPI_Allreduce(MPI_IN_PLACE, yoffsets.data(), npes, MPI_INT, MPI_SUM, G.comm);
    MPI_Allreduce(MPI_IN_PLACE, zoffsets.data(), npes, MPI_INT, MPI_SUM, G.comm);

    recvcounts.resize(npes, 0);
    recvcounts[my_rank] = pbasis;
    MPI_Allreduce(MPI_IN_PLACE, recvcounts.data(), npes, MPI_INT, MPI_SUM, G.comm);

    recvoffsets.resize(npes, 0);
    irecvoffsets.resize(npes, 0);

    for(int idx=1;idx < npes;idx++) recvoffsets[idx] = recvoffsets[idx-1] + (size_t)recvcounts[idx-1];
    for(int idx=1;idx < npes;idx++) irecvoffsets[idx] = irecvoffsets[idx-1] + recvcounts[idx-1];

    vexx_RMS.resize(ct.max_exx_steps, 0.0);
    nstates_occ = 0;
    for(int st=0;st < nstates;st++) if(occ[st] > 1.0e-6) nstates_occ++;
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
         ct.gamma_extrapolation = false;
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

      //  if(!pwave->gmask[ig]) continue;
        double fac = 1.0;
        if (ct.gamma_extrapolation)
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
        size_t pstop = (size_t)nstates * (size_t)pbasis;
        double *prev_vexx = &psi[pstop];
        for(size_t idx=0;idx < pstop;idx++) prev_vexx[idx] = vexx[idx];
        for(size_t idx=0;idx < pstop;idx++) vexx[idx] = 0.0;

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

        double tvexx_RMS = 0.0;
        pstop = (size_t)nstates_occ * (size_t)pbasis;
        for(size_t idx=0;idx < pstop;idx++) tvexx_RMS += (prev_vexx[idx] - vexx[idx])*(prev_vexx[idx] - vexx[idx]);
        MPI_Allreduce(MPI_IN_PLACE, &tvexx_RMS, 1, MPI_DOUBLE, MPI_SUM, this->G.comm);
        double vscale = (double)nstates_occ*(double)G.get_GLOBAL_BASIS(1);
        vexx_RMS[ct.exx_steps] = sqrt(tvexx_RMS / vscale);
        ct.vexx_rms = vexx_RMS[ct.exx_steps];

    }
    else
    {

        // Map my portion of vexx into 
        int my_rank = G.get_rank();
        int npes = G.get_NPES();

        double *arbuf, *atbuf;
        MPI_Alloc_mem(pwave->pbasis*sizeof(double), MPI_INFO_NULL, &arbuf);
        MPI_Alloc_mem(pbasis*sizeof(double), MPI_INFO_NULL, &atbuf);

        double *vexx_global = new double[pwave->pbasis]();

        // Write serial wavefunction files. May need to do some numa optimization here at some point
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
        MPI_Request req=MPI_REQUEST_NULL;
        MPI_Status mrstatus;
        double tvexx_RMS = 0.0;
        int flag=0;

        for(int i=0;i < nstates;i++)
        {
            double *psi_i = (double *)&psi_s[(size_t)i*(size_t)pwave->pbasis];
            RmgTimer *RT1 = new RmgTimer("5-Functional: Exx potential fft");
#pragma omp parallel for schedule(dynamic)
            for(int j=0;j<nstates_occ;j++)
            {
                int omp_tid = omp_get_thread_num();
                std::complex<double> *p = pvec[omp_tid];
                std::complex<float> *w = wvec[omp_tid];
                if(my_rank == (j % npes))
                {
                    double *psi_j = (double *)&psi_s[(size_t)j*(size_t)pwave->pbasis];
                    if(use_float_fft)
                        fftpair(psi_i, psi_j, p, w, gfac);
                    else
                        fftpair(psi_i, psi_j, p, gfac);
                    // We can speed this up by adding more critical sections if it proves to be a bottleneck
#pragma omp critical(part3)
                    {
                            for(size_t idx = 0;idx < (size_t)pwave->pbasis;idx++) 
                                vexx_global[idx] += scale * std::real(p[idx]) * psi_j[idx];
                    }
                }
                if(!omp_tid) MPI_Test(&req, &flag, &mrstatus);
            }
            delete RT1;

            // We wait for communication from previous row to finish and then copy it into place
            MPI_Wait(&req, &mrstatus);
            if(i)
            {
                if(i < nstates_occ)
                {
                    double t1 = 0.0;
                    double *old_vexx = &vexx[(size_t)(i-1) * (size_t)pbasis];
                    for(int idx=0;idx < pbasis;idx++) t1 += (atbuf[idx] - old_vexx[idx])*(atbuf[idx] - old_vexx[idx]);
                    tvexx_RMS += t1;
                }
                memcpy(&vexx[(size_t)(i-1) * (size_t)pbasis], atbuf, pbasis * sizeof(double));
            }
            // Remap so we can use MPI_Reduce_scatter
            Remap(vexx_global, arbuf);

            // Zero out vexx_global so it can be used for accumulation in the next iteration of the loop.
            std::fill(vexx_global, vexx_global + pwave->pbasis, 0.0);
            MPI_Ireduce_scatter(arbuf, atbuf, recvcounts.data(), MPI_DOUBLE, MPI_SUM, G.comm, &req);
        }

        // Wait for last transfer to finish and then copy data to correct location
        MPI_Wait(&req, &mrstatus);
        memcpy(&vexx[(size_t)(nstates-1) * (size_t)pbasis], atbuf, pbasis * sizeof(double));

        double scale = (double)nstates_occ*(double)G.get_GLOBAL_BASIS(1);
        MPI_Allreduce(MPI_IN_PLACE, &tvexx_RMS, 1, MPI_DOUBLE, MPI_SUM, this->G.comm);
        vexx_RMS[ct.exx_steps] += sqrt(tvexx_RMS / scale);
        ct.vexx_rms = vexx_RMS[ct.exx_steps];

        MPI_Barrier(G.comm);

        // munmap wavefunction array
        munmap(psi_s, length);
        close(serial_fd);

        delete [] vexx_global;
        MPI_Free_mem(atbuf);
        MPI_Free_mem(arbuf);
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

template <> int Exxbase<double>::VxxIntChol(std::vector<double> &mat, std::vector<double> &CholVec, int cmax, int nst_occ)
{
    
    double tol = 1.0e-5;
    std::vector<double> m_diag, m_appr, m_nu0, delta;
    int nst2 = nst_occ * nst_occ;
    int nst2_perpe = (nst2 + pct.grid_npes-1) /pct.grid_npes;
    int nu, ione = 1;
    m_diag.resize(nst2_perpe);
    m_nu0.resize(nst2_perpe);
    m_appr.assign(nst2_perpe, 0.0);
    delta.assign(nst2_perpe * pct.grid_npes, 0.0);

    for(int i =0; i<nst2; i++) 
    {
        int i_dist = i - nst2_perpe * pct.gridpe;
        if(i_dist >= 0 && i_dist < nst2_perpe)
        {
            m_diag[i_dist] = mat[i*nst2_perpe+i_dist];
            delta[i] = mat[i*nst2_perpe+i_dist];
        }
    }
    
    MPI_Allreduce(MPI_IN_PLACE, delta.data(), nst2, MPI_DOUBLE, MPI_SUM, pct.grid_comm);

    nu = 0;
    double delta_max = 0.0;
    for(int i =0; i<nst2; i++) 
    {
        if(delta_max < std::abs(delta[i])) 
        {
            delta_max = std::abs(delta[i]);
            nu = i;
        }
    }

    delta_max = 1.0/std::sqrt(delta_max);

    dcopy(&nst2_perpe, &mat.data()[nu * nst2_perpe], &ione, CholVec.data(), &ione);
    dscal(&nst2_perpe, &delta_max, CholVec.data(), &ione);

    int ic;
    double *CholVecNu = new double[cmax * nst_occ];
    for(ic= 0; ic < cmax * nst_occ - 1; ic++)
    {
        for(int i = 0; i < nst2_perpe * pct.grid_npes; i++) delta[i] = 0.0;

        for(int i = 0; i < nst2_perpe; i++)
        {
            m_appr[i] += CholVec[ic*nst2_perpe + i] * MyConj(CholVec[ic*nst2_perpe+i]);
            delta[i + pct.gridpe * nst2_perpe] = m_appr[i] - m_diag[i];
        }

        MPI_Allreduce(MPI_IN_PLACE, delta.data(), nst2, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
        delta_max = 0.0;
        for(int i =0; i<nst2; i++) 
        {
            //            printf("\n  aaa %d %d %f ", ic, i, m_diag[i]-m_appr[i]);
            if(delta_max < std::abs(delta[i]))
            {
                delta_max = std::abs(delta[i]);
                nu = i;
            }
        }

        if(ct.verbose && pct.gridpe == 0) printf("\n delta_max %d %d %f\n", ic, nu, delta_max);
        if(delta_max < tol) break;
        delta_max = std::sqrt(delta_max);

        int root = nu/nst2_perpe;
        if(pct.gridpe == root)
        {
            for(int ics = 0; ics < ic+1; ics++)
                CholVecNu[ics] = CholVec[ics * nst2_perpe + nu % nst2_perpe];
        }
        
        MPI_Bcast(CholVecNu, ic+1, MPI_DOUBLE, root, pct.grid_comm);

        for(int i =0; i<nst2_perpe; i++) m_nu0[i] = 0.0;
        for(int ics = 0; ics < ic+1; ics++)
        {
            for(int i =0; i<nst2_perpe; i++)
            {
                m_nu0[i] += MyConj(CholVecNu[ics]) * CholVec[ics * nst2_perpe + i];
            }
        } 

        for(int i =0; i<nst2_perpe; i++)
        {

            CholVec[(ic+1) * nst2_perpe + i] = (mat[nu * nst2_perpe  + i] - m_nu0[i])/delta_max;
        }

    }



    CholVec.erase(CholVec.begin()+(ic+1)*nst2_perpe, CholVec.end());

    return ic+1;

}
template <> int Exxbase<std::complex<double>>::VxxIntChol(
        std::vector<std::complex<double>> &ExxI, std::vector<std::complex<double>> &CholVec, 
        int cmax, int nstates_occ)
{
    return 0;
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
        double *psi_k = (double *)&psi_s[k*pwave->pbasis];
        double *psi_l = (double *)&psi_s[l*pwave->pbasis];
        for(size_t idx=0;idx < pwave->pbasis;idx++) kl_pair[ie*pwave->pbasis + idx] = psi_k[idx]*psi_l[idx];
    }

    delete RT0;
    // Now compute integrals for (i, j, k, l)

    int ij_length = ij_end - ij_start;
    int kl_length = kl_end - kl_start;
    // Now matrix multiply to produce a block of (1,jblocks, 1, nstates_occ) results
    RT0 = new RmgTimer("5-Functional: Exx: gemm");
    alpha = L.get_omega() / ((double)(G.get_NX_GRID(1) * G.get_NY_GRID(1) * G.get_NZ_GRID(1)));
    RmgGemm(trans_a, trans_n, ij_length, kl_length, pwave->pbasis, alpha, ij_pair, pwave->pbasis, kl_pair, pwave->pbasis, beta, Exxints, ij_length);
    delete RT0;
    int pairsize = ij_length * kl_length;
    RT0 = new RmgTimer("5-Functional: Exx: reduce");

    MPI_Allreduce(Exxints, Summedints, pairsize, MPI_DOUBLE, MPI_SUM, LG->comm);
    MPI_Barrier(LG->comm);
    delete RT0;

    RT0 = new RmgTimer("5-Functional: Exx: print");
    size_t nst2 = nstates_occ * nstates_occ;
    int nst2_perpe = (nst2 + pct.grid_npes-1) /pct.grid_npes;
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

            if(ct.ExxIntChol && mode == EXX_DIST_FFT)
            {

                int ij = i * nstates_occ + j;
                int ji = j * nstates_occ + i;
                int ij_dist = i * nstates_occ + j - nst2_perpe * pct.gridpe;
                int ji_dist = j * nstates_occ + i - nst2_perpe * pct.gridpe;

                int kl = k * nstates_occ + l;
                int lk = l * nstates_occ + k;
                int kl_dist = k * nstates_occ + l - nst2_perpe * pct.gridpe;
                int lk_dist = l * nstates_occ + k - nst2_perpe * pct.gridpe;

                if(ij_dist >= 0 && ij_dist < nst2_perpe)
                {
                    ExxInt[lk * nst2_perpe + ij_dist] = Summedints[ie * ij_length + ic];
                    ExxInt[kl * nst2_perpe + ij_dist] = Summedints[ie * ij_length + ic];
                }

                if(ji_dist >= 0 && ji_dist < nst2_perpe)
                {
                    ExxInt[lk * nst2_perpe + ji_dist] = Summedints[ie * ij_length + ic];
                    ExxInt[kl * nst2_perpe + ji_dist] = Summedints[ie * ij_length + ic];
                }

                if(kl_dist >= 0 && kl_dist < nst2_perpe)
                {
                    ExxInt[ij * nst2_perpe + kl_dist] = Summedints[ie * ij_length + ic];
                    ExxInt[ji * nst2_perpe + kl_dist] = Summedints[ie * ij_length + ic];
                }

                if(lk_dist >= 0 && lk_dist < nst2_perpe)
                {
                    ExxInt[ij * nst2_perpe + lk_dist] = Summedints[ie * ij_length + ic];
                    ExxInt[ji * nst2_perpe + lk_dist] = Summedints[ie * ij_length + ic];
                }
            }
            else if(ct.ExxIntChol && mode == EXX_LOCAL_FFT)
            {

                int idx = (i * nstates_occ+j) * nstates_occ * nstates_occ + k * nstates_occ + l;
                ExxInt[idx]  = Summedints[ie * ij_length + ic];
                idx = (j * nstates_occ+i) * nstates_occ * nstates_occ + k * nstates_occ + l;
                ExxInt[idx]  = Summedints[ie * ij_length + ic];
                idx = (i * nstates_occ+j) * nstates_occ * nstates_occ + l * nstates_occ + k;
                ExxInt[idx]  = Summedints[ie * ij_length + ic];
                idx = (j * nstates_occ+i) * nstates_occ * nstates_occ + l * nstates_occ + k;
                ExxInt[idx]  = Summedints[ie * ij_length + ic];

                idx = (k * nstates_occ+l) * nstates_occ * nstates_occ + i * nstates_occ + j;
                ExxInt[idx]  = Summedints[ie * ij_length + ic];
                idx = (k * nstates_occ+l) * nstates_occ * nstates_occ + j * nstates_occ + i;
                ExxInt[idx]  = Summedints[ie * ij_length + ic];
                idx = (l * nstates_occ+k) * nstates_occ * nstates_occ + i * nstates_occ + j;
                ExxInt[idx]  = Summedints[ie * ij_length + ic];
                idx = (l * nstates_occ+k) * nstates_occ * nstates_occ + j * nstates_occ + i;
                ExxInt[idx]  = Summedints[ie * ij_length + ic];
            }

            else if(LG->get_rank()==0 && kl >= ij)fprintf(fp, "%14.8e %d %d %d %d\n", Summedints[ie * ij_length + ic], i, j, k, l);
        }
    }

    delete RT0;
}

// This computes exact exchange integrals
// and writes the result into vfile.
template <> void Exxbase<double>::Vexx_integrals(std::string &vfile)
{

    double scale = 1.0 / (double)pwave->global_basis;


    if(ct.ExxIntChol)
    {
        // std::string filename = wavefile + "ExxInt" + "_spin"+std::to_string(pct.spinpe);
        //  exxint_fd = open(filename.c_str(), O_RDWR|O_CREAT|O_TRUNC, (mode_t)0600);
        //  if(exxint_fd < 0)
        //      throw RmgFatalException() << "Error! Could not open " << filename << " . Terminating.\n";

        //  size_t length = nstates_occ * nstates_occ * nstates_occ * nstates_occ *sizeof(double);
        //  ExxInt = (double *)mmap(NULL, length, PROT_READ, MAP_PRIVATE, exxint_fd, 0);

        size_t nst2 = nstates_occ * nstates_occ;
        size_t nst2_perpe = (nst2 + pct.grid_npes -1)/pct.grid_npes;
        size_t length = nst2 * nst2_perpe;
        if(mode == EXX_LOCAL_FFT)
            length = nst2 * nst2;
        ExxInt.resize(length, 0.0);
        length = nstates_occ * ct.exxchol_max * nst2_perpe;
        ExxCholVec.resize(length);

    }
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

    kl_pair = new double[block_size * pwave->pbasis];
    ij_pair = new double[block_size * pwave->pbasis];

    Exxints = new double[block_size * block_size];
    Summedints = new double[block_size * block_size];

    wf_fft = new std::complex<double> [pwave->pbasis];


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
        std::string filename = wavefile + "_spin"+std::to_string(pct.spinpe) + "_kpt0";
        serial_fd = open(filename.c_str(), O_RDONLY, (mode_t)0600);
        if(serial_fd < 0)
            throw RmgFatalException() << "Error! Could not open " << filename << " . Terminating.\n";

        size_t length = nstates * pwave->pbasis *sizeof(double);
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
            double *psi_i = (double *)&psi_s[i*pwave->pbasis];
            double *psi_j = (double *)&psi_s[j*pwave->pbasis];
            fftpair(psi_i, psi_j, wf_fft, gfac);

            // store (i,j) fft pair in ij_pair
            // forward and then backward fft should not have the scale, Wenchang
            for(size_t idx=0;idx < pwave->pbasis;idx++) ij_pair[ic * pwave->pbasis + idx] = scale * std::real(wf_fft[idx]);

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

    if(ct.ExxIntChol)
    {
        int length = ExxInt.size();
        size_t nst2 = nstates_occ * nstates_occ;
        size_t nst2_perpe = (nst2 + pct.grid_npes -1)/pct.grid_npes;
        if(mode == EXX_LOCAL_FFT)
        {
            MPI_Allreduce(MPI_IN_PLACE, ExxInt.data(), length, MPI_DOUBLE, MPI_SUM, pct.grid_comm);

            std::vector<double> ExxIntGlob;
            ExxIntGlob = ExxInt;
            ExxInt.resize(nst2 * nst2_perpe);

            for(size_t i = 0; i < nst2; i++)
            {
                for(size_t j = 0; j < nst2_perpe; j++)
                {
                    size_t j_glob = pct.gridpe * nst2_perpe + j;
                    if( j_glob < nst2)
                        ExxInt[i * nst2_perpe + j] = ExxIntGlob[i* nst2 + j_glob];
                }
            }

        }
        Nchol = VxxIntChol(ExxInt, ExxCholVec, ct.exxchol_max, nstates_occ);
        std::vector<double> eigs;

        eigs.resize(nstates_occ, 0.0);
        for(int st = 0; st < nstates_occ; st++)
        {
            eigs[st] = ct.kp[0].eigs[st];
        }

        std::vector<double> ExxCholVecGlob;
        ExxCholVecGlob.resize(Nchol * nst2, 0.0);
        for(int i = 0; i < Nchol; i++)
        {
            for(size_t j = 0; j < nst2_perpe; j++)
            {
                size_t j_glob = pct.gridpe * nst2_perpe + j;
                ExxCholVecGlob[j_glob * Nchol + i] = ExxCholVec[i * nst2_perpe + j];

            }
        }

        MPI_Allreduce(MPI_IN_PLACE, ExxCholVecGlob.data(), Nchol*nst2, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
        
        if(pct.worldrank == 0)
            WriteForAFQMC(nstates_occ, Nchol, nstates_occ, nstates_occ, eigs, ExxCholVecGlob, Hcore);
    }
}

template <> void Exxbase<std::complex<double>>::Vexx_integrals(std::string &vfile)
{
    printf("Exx mode not programmed yet\n");
}

template <class T> Exxbase<T>::~Exxbase(void)
{

    if(mode == EXX_DIST_FFT) return;

    if(ct.ExxIntChol)
    {
        //close(exxint_fd);
        //size_t length = nstates_occ * nstates_occ * nstates_occ * nstates_occ * sizeof(T);
        //munmap(ExxInt, length);
        ExxInt.clear();
        ExxCholVec.clear();
    }
    close(serial_fd);
    size_t length = nstates * pwave->pbasis * sizeof(T);
    munmap(psi_s, length);
    //std::string filename= wavefile + "_spin"+ std::to_string(pct.spinpe);
    //unlink(filename.c_str());

    GpuFreeManaged(gfac);
    delete LG;
    MPI_Comm_free(&lcomm); 
    delete pwave;
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

    for(int ik = 0; ik < ct.num_kpts_pe; ik++)
    {
        int ik_glob = ik + pct.kstart;
        std::string filename = wavefile + "_spin"+std::to_string(pct.spinpe) + "_kpt" + std::to_string(ik_glob);
        int flag = MPI_File_open(G.comm, filename.c_str(), amode, fileinfo, &mpi_fhand);
        if(flag) 
            throw RmgFatalException() << "Error! Could not open " << filename << " . Terminating.\n";
        MPI_Offset disp = 0;

        T *wfptr;
        MPI_File_set_view(mpi_fhand, disp, wftype, grid_c, "native", MPI_INFO_NULL);
        int dis_dim = G.get_P0_BASIS(1);
        for(int st=0;st < nstates;st++)
        {
            wfptr = &psi[ik * ct.max_states * dis_dim + st * dis_dim];
            MPI_File_read_all(mpi_fhand, wfptr, dis_dim, wftype, &status);
        }
        MPI_Barrier(G.comm);
        MPI_File_close(&mpi_fhand);
    }
    MPI_Barrier(G.comm);
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
            wfptr = &psi[(ik * ct.max_states * dis_dim + st * dis_dim) * ct.noncoll_factor];
            MPI_File_write_all(mpi_fhand, wfptr, dis_dim, wftype, &status);
            if(ct.noncoll)
            {
                wfptr += dis_dim;
                MPI_File_write_all(mpi_fhand, wfptr, dis_dim, wftype, &status);
            }
        }
        MPI_Barrier(G.comm);
        MPI_File_close(&mpi_fhand);
    }
    MPI_Type_free(&grid_c);
    fflush(NULL);
    MPI_Barrier(G.comm);
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

    std::complex<double> *arbuf, *atbuf;
    MPI_Alloc_mem(2*pwave->pbasis*sizeof(double), MPI_INFO_NULL, &arbuf);
    MPI_Alloc_mem(2*pbasis*sizeof(double), MPI_INFO_NULL, &atbuf);

    std::complex<double> *vexx_global = new std::complex<double>[pwave->pbasis]();

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

        // Save a copy of vexx in the upper part of the orbital storage
        std::complex<double> *prev_vexx = &psi[(size_t)ik * (size_t)ct.max_states * (size_t)pbasis + (size_t)nstates * (size_t)pbasis];
        size_t pstop = (size_t)nstates * (size_t)pbasis;
        std::complex<double> *cur_vexx = &vexx[(size_t)ik * (size_t)ct.run_states * (size_t)pbasis];
        for(size_t idx=0;idx < pstop;idx++)
        {
            prev_vexx[idx] = vexx[(size_t)ik * (size_t)ct.run_states * (size_t)pbasis + idx];
            cur_vexx[idx] = 0.0;
        }


        MPI_Barrier(G.comm);

        for(int iq = 0; iq < ct.klist.num_k_all; iq++)
        {

            int ikindex = ct.klist.k_map_index[iq];
            int isym = ct.klist.k_map_symm[iq];
            int isyma = std::abs(isym) -1;


            kq[0] = ct.kp[ik_glob].kpt[0] - ct.klist.k_all_xtal[iq][0];
            kq[1] = ct.kp[ik_glob].kpt[1] - ct.klist.k_all_xtal[iq][1];
            kq[2] = ct.kp[ik_glob].kpt[2] - ct.klist.k_all_xtal[iq][2];
            //   while(kq[0] > 0.5) kq[0] -= 1.0;
            //   while(kq[1] > 0.5) kq[1] -= 1.0;
            //   while(kq[2] > 0.5) kq[2] -= 1.0;
            //   while(kq[0] <-0.5) kq[0] += 1.0;
            //   while(kq[1] <-0.5) kq[1] += 1.0;
            //   while(kq[2] <-0.5) kq[2] += 1.0;
            double v0, v1, v2;

            v0 = kq[0] *Rmg_L.b0[0] + kq[1] *Rmg_L.b1[0] + kq[2] *Rmg_L.b2[0];
            v1 = kq[0] *Rmg_L.b0[1] + kq[1] *Rmg_L.b1[1] + kq[2] *Rmg_L.b2[1];
            v2 = kq[0] *Rmg_L.b0[2] + kq[1] *Rmg_L.b1[2] + kq[2] *Rmg_L.b2[2];

            kq[0] = v0 * twoPI;
            kq[1] = v1 * twoPI;
            kq[2] = v2 * twoPI;


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

                        symm_ijk(&Rmg_Symm->sym_rotate[isyma *9], &Rmg_Symm->ftau_wave[isyma*3], ix, iy, iz, ixx, iyy, izz, nx_grid, ny_grid, nz_grid);

                        if(isym >= 0)
                        {
                            for(int st = 0; st < nstates_occ; st++)
                            {
                                psi_q[st * nbasis + ix * ny_grid * nz_grid + iy * nz_grid + iz]
                                    = (psi_q_map[st * nbasis + ixx * ny_grid * nz_grid + iyy * nz_grid + izz]);
                            }
                        }
                        else
                        {
                            for(int st = 0; st < nstates_occ; st++)
                            {
                                psi_q[st * nbasis + ix * ny_grid * nz_grid + iy * nz_grid + iz]
                                    = std::conj(psi_q_map[st * nbasis + ixx * ny_grid * nz_grid + iyy * nz_grid + izz]);
                            }
                        }
                    }
                }
            }


            MPI_Request req=MPI_REQUEST_NULL;
            MPI_Status mrstatus;
            int flag=0;

            std::fill(vexx_global, vexx_global + pwave->pbasis, 0.0);
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
                    if(my_rank == (j % npes))
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
                                vexx_global[idx] += scale * p[idx] * psi_j[idx] / (double)ct.klist.num_k_all;
                        }
                    }
                    if(!omp_tid) MPI_Test(&req, &flag, &mrstatus);
                }
                delete RT1;

                // We wait for communication from previous row to finish and then copy it into place
                MPI_Wait(&req, &mrstatus);
                if(i)
                {
                    for(size_t idx=0;idx < (size_t)pbasis;idx++) 
                    {
                        vexx[(size_t)ik * (size_t)ct.run_states * (size_t)pbasis + (size_t)(i-1) * (size_t)pbasis + idx] += atbuf[idx];
                    } 
                }
                // Remap so we can use MPI_Reduce_scatter
                Remap(vexx_global, arbuf);

                // Zero out vexx_global so it can be used for accumulation in the next iteration of the loop.
                std::fill(vexx_global, vexx_global + pwave->pbasis, 0.0);
                MPI_Ireduce_scatter(arbuf, atbuf, recvcounts.data(), MPI_DOUBLE_COMPLEX, MPI_SUM, G.comm, &req);

            }

            // Wait for last transfer to finish and then copy data to correct location
            MPI_Wait(&req, &mrstatus);
            for(size_t idx=0;idx < (size_t)pbasis;idx++) 
            {
                vexx[(size_t)ik * (size_t)ct.run_states * (size_t)pbasis + (size_t)(nstates-1) * (size_t)pbasis + idx] += atbuf[idx];
            } 


            munmap(psi_q_map, length);
            close(serial_fd_q);
        }

        MPI_Barrier(G.comm);

        // munmap wavefunction array
        munmap(psi_s, length);
        close(serial_fd);

        double tvexx_RMS = 0.0;
        size_t length = (size_t)nstates_occ * (size_t)pbasis;
        cur_vexx = &vexx[(size_t)ik * (size_t)ct.run_states * (size_t)pbasis];
        for(size_t idx=0;idx < length;idx++) tvexx_RMS += std::real((prev_vexx[idx] - cur_vexx[idx])*std::conj(prev_vexx[idx] - cur_vexx[idx]));
        MPI_Allreduce(MPI_IN_PLACE, &tvexx_RMS, 1, MPI_DOUBLE, MPI_SUM, this->G.comm);
        double vscale = (double)nstates_occ*(double)G.get_GLOBAL_BASIS(1);

        vexx_RMS[ct.exx_steps] += sqrt(tvexx_RMS / vscale) * ct.kp[ik_glob].kweight;

    } // end loop over k-points

    double t1 = vexx_RMS[ct.exx_steps];
    MPI_Allreduce(MPI_IN_PLACE, &t1, 1, MPI_DOUBLE, MPI_SUM, pct.kpsub_comm);
    vexx_RMS[ct.exx_steps] = t1;
    ct.vexx_rms = t1;

    delete [] vexx_global;

    MPI_Free_mem(atbuf);
    MPI_Free_mem(arbuf);

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
    if(ct.exxdiv_treatment == EXX_DIV_NONE) return;

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
                    if (ct.gamma_extrapolation)
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

    if(!ct.gamma_extrapolation)
    {
        if(scr_type == ERFC_SCREENING)
            exxdiv += 1.0/4.0/erfc_scrlen/erfc_scrlen;
        else if( std::abs(yukawa) > eps )
            exxdiv += 1.0/yukawa;
        else
            exxdiv += -alpha;
    }
    exxdiv *= fourPI /ct.klist.num_k_all;



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

    exxdiv *= ct.klist.num_k_all;
    if(ct.verbose && pct.gridpe == 0) 
    {
        printf("\n exxdiv = %f %f %f", exxdiv, aa, alpha);
        printf("\n erfc_scrlen = %f", erfc_scrlen);
    }
}


template void Exxbase<double>::Remap(double *, double *);
template void Exxbase<std::complex<double>>::Remap(std::complex<double> *, std::complex<double> *);
template <class T> void Exxbase<T>::Remap(T *inbuf, T *rbuf)
{
    // Remaps so we can use MPI_Reduce_scatter rather than MPI_Allreduce
    int npes = G.get_NPES();
    int gdimy = G.get_NY_GRID(1);
    int gdimz = G.get_NZ_GRID(1);

#pragma omp parallel for
    for(size_t rank=0;rank < (size_t)npes;rank++)
    {
        size_t dimx_r = (size_t)dimsx[rank]; 
        size_t dimy_r = (size_t)dimsy[rank]; 
        size_t dimz_r = (size_t)dimsz[rank]; 
        size_t offset_r = recvoffsets[rank];
        size_t xoffset_r = (size_t)xoffsets[rank];
        size_t yoffset_r = (size_t)yoffsets[rank];
        size_t zoffset_r = (size_t)zoffsets[rank];
        for(size_t ix=0;ix < dimx_r;ix++)
        {
            for(size_t iy=0;iy < dimy_r;iy++)
            {
                for(size_t iz=0;iz < dimz_r;iz++)
                {
                    rbuf[offset_r + ix*dimy_r*dimz_r + iy*dimz_r + iz] = 
                        inbuf[(ix+xoffset_r)*(size_t)gdimy*(size_t)gdimz + (iy+(size_t)yoffset_r)*(size_t)gdimz + iz + (size_t)zoffset_r];
                }
            }
        }
    }
}

template void Exxbase<double>::SetHcore(double *Hij, int lda);
template void Exxbase<std::complex<double>>::SetHcore(std::complex<double> *Hij, int lda);
template <class T> void Exxbase<T>::SetHcore(T *Hij, int lda)
{
    Hcore.resize(nstates_occ * nstates_occ);
    for(int i = 0; i < nstates_occ; i++)
        for(int j = 0; j < nstates_occ; j++)
            Hcore[i * nstates_occ + j] = Hij[i* lda + j];

}
