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
#include "Exx_on.h"
#include "prototypes_on.h"
#include "blas.h"

// This class implements exact exchange for ON module using the Exxbase clase
// DelocalOrbital: is delocalized and store K_kl * Phi_l(r), 
// K_kl is the charge density matrix, C_ki * f(i) C_il =  mat_X from DaigScalapack.cpp
// K_kl dimension: num_tot * num_thispe for LocalOrbital


template Exx_on<double>::Exx_on(BaseGrid &, BaseGrid &, Lattice &, const std::string &, LocalObject<double> &phi, double *, int);
template Exx_on<std::complex<double>>::Exx_on(BaseGrid &, BaseGrid &, Lattice &, const std::string &, 
    LocalObject<std::complex<double>> &phi, double *, int );

template Exx_on<double>::~Exx_on(void);
template Exx_on<std::complex<double>>::~Exx_on(void);
template <class T> Exx_on<T>::~Exx_on()
{
};

template <class T> Exx_on<T>::Exx_on (
          BaseGrid &G_in,
          BaseGrid &G_h_in,
          Lattice &L_in,
          const std::string &wavefile_in,
          LocalObject<T> &Phi_in, double *occ_in, 
          int mode_in) : G(G_in), G_h(G_h_in), L(L_in), wavefile(wavefile_in), Phi(Phi_in), occ(occ_in), mode(mode_in)
{
    RmgTimer RT0("5-Functional: Exx init");
    DePhi = (T *)RmgMallocHost(Phi.num_tot * Phi.pbasis * sizeof(T));
    Exxb = new Exxbase<T>(G, G_h, L, wavefile, Phi.num_tot, occ, DePhi, ct.exx_mode);
    Omega_j = (T *)RmgMallocHost(Phi.num_tot * Phi.pbasis * sizeof(T));
    PreOrbital = (T *)RmgMallocHost(Phi.num_thispe * Phi.pbasis * sizeof(T));
    Xij_mat =  (T *)RmgMallocHost(Phi.num_tot * Phi.num_tot * sizeof(T));
    Xij_mat_Sinv =  (T *)RmgMallocHost(Phi.num_tot * Phi.num_tot * sizeof(T));
    for(int i = 0; i < Phi.num_tot * Phi.num_tot; i++) Xij_mat[i] = 0.0;
    for(int i = 0; i < Phi.num_tot * Phi.pbasis; i++) Omega_j[i] = 0.0;
}

template <> void Exx_on<double>::Omega(double *rho_matrix, bool use_float_fft)
{
    double one = 1.0, zero = 0.0;
    int num_tot = Phi.num_tot;
    int num_thispe = Phi.num_thispe;
    int pbasis = Phi.pbasis;
    RmgTimer RT0("5-Functional: Exx potential");

    RmgGemm("N", "N", pbasis, num_tot, num_thispe, one, Phi.storage_cpu, pbasis,
            rho_matrix, num_thispe, zero, DePhi, pbasis);


    double scale = - 1.0 / (double)Exxb->pwave->global_basis;

    // Clear vexx
    for(int idx=0;idx < Phi.num_tot*pbasis;idx++) Omega_j[idx] = 0.0;


    double *zeropsi = new double[pbasis]();
    std::vector<std::complex<double> *> pvec;
    std::vector<std::complex<float> *> wvec;
    pvec.resize(ct.OMP_THREADS_PER_NODE);
    wvec.resize(ct.OMP_THREADS_PER_NODE);
    for(int tid=0;tid < ct.OMP_THREADS_PER_NODE;tid++)
    {
        pvec[tid] = (std::complex<double> *)fftw_malloc(sizeof(std::complex<double>) * Exxb->pwave->pbasis);
        wvec[tid] = (std::complex<float> *)fftw_malloc(sizeof(std::complex<float>) * Exxb->pwave->pbasis);
    }
    
    if(mode != EXX_DIST_FFT)
    {
        throw RmgFatalException() << "only dist_fft mode for Hybrid now" << "\n";
    }
    // Loop over v_kj, here k, and j are all orbitals
    // v_kj = sum_r DePhi_k(r) * Phi_j(r)
    for(int j=0;j < num_tot;j++)
    {   
        int j_local = Phi.index_global_to_proj[j];
        double *phi_j, *Ome_j;
        if(j_local < 0)
        {
            phi_j = zeropsi;
        }
        else
        {
            phi_j = &Phi.storage_cpu[j_local * pbasis];
        }
        Ome_j = &Omega_j[j * pbasis];
#pragma omp parallel for schedule(dynamic)
        for(int k=0; k < num_tot;k++)
        {
            double *psi_k = (double *)&DePhi[k*pbasis];
            int k_local = Phi.index_global_to_proj[k];
            double *phi_k;
            if(k_local < 0)
                phi_k = zeropsi;
            else
                phi_k = &Phi.storage_cpu[k_local * pbasis];

            int omp_tid = omp_get_thread_num();
            std::complex<double> *p = pvec[omp_tid];
            std::complex<float> *w = wvec[omp_tid];

            RmgTimer RT1("5-Functional: Exx potential fft");

            if(use_float_fft)
            {
                for(int idx = 0; idx < pbasis; idx++) w[idx] = std::complex<float>(psi_k[idx] * phi_j[idx], 0.0);
                Exxb->pwave->FftForward(w, w);
                for(int ig=0;ig < pbasis;ig++) w[ig] *= Exxb->gfac[ig];
                Exxb->pwave->FftInverse(w, w);
                for(int ig=0;ig < pbasis;ig++) p[ig] = w[ig];
            }
            else
            {
                for(int idx = 0; idx < pbasis; idx++) p[idx] = std::complex<double>(psi_k[idx] * phi_j[idx], 0.0);
                Exxb->pwave->FftForward(p, p);
                for(int ig=0;ig < pbasis;ig++) p[ig] *= Exxb->gfac[ig];
                Exxb->pwave->FftInverse(p, p);
            }

#pragma omp critical(part1)
            {               
                if(k_local >= 0)
                    for(int idx = 0;idx < pbasis;idx++)Ome_j[idx] += scale * std::real(p[idx]) * phi_k[idx];
            }


        }
    }

    for(int tid=0;tid < ct.OMP_THREADS_PER_NODE;tid++)
    {
        fftw_free(wvec[tid]);
        fftw_free(pvec[tid]);
    }
    delete [] zeropsi;
}

template <> void Exx_on<std::complex<double>>::Omega(std::complex<double> *rho_matrix, bool use_fft_float)
{
}

template <> void Exx_on<std::complex<double>>::Xij(std::complex<double> *Sij_inverse, LocalObject<std::complex<double>> &Phi){};
template <> void Exx_on<double>::Xij(double *Sij_inverse, LocalObject<double> &Phi)
{
    
    double t1 = Rmg_G->get_NX_GRID(Phi.density);
    t1 *= Rmg_G->get_NY_GRID(Phi.density);
    t1 *= Rmg_G->get_NZ_GRID(Phi.density);

    double vol = Rmg_L.get_omega() /t1;

    int na = Phi.num_thispe ;
    int ntot = Phi.num_tot; 

    int pbasis = Phi.pbasis;

    double zero = 0.0, one = 1.0;

    double *mat_tem = (double *)RmgMallocHost(ntot*ntot*sizeof(double));
    RmgGemm("C", "N", na, ntot, pbasis, vol, Phi.storage_cpu, pbasis,
               (double *)Omega_j, pbasis, zero, mat_tem, na);
    for(int idx = 0; idx < na * ntot; idx++) Xij_mat[idx] = 0.0;

    for(int st1 = 0; st1 < ntot; st1++)
    {
        int st1_local = Phi.index_global_to_proj[st1];
        if (st1_local < 0 ) continue;

        for(int st2 = 0; st2 < ntot; st2++)
        {
            Xij_mat[st2 * ntot + st1] = mat_tem[st2 * na + st1_local];
        }

    }

    size_t size = ntot * ntot;
    MPI_Allreduce(MPI_IN_PLACE, Xij_mat, size, MPI_DOUBLE, MPI_SUM, Phi.comm);

    RmgGemm("N", "N", ntot, ntot, ntot, one, Sij_inverse, ntot,
            Xij_mat, ntot, zero, mat_tem, ntot);
    RmgGemm("N", "N", ntot, ntot, ntot, one, mat_tem, ntot,
            Sij_inverse, ntot, zero, Xij_mat_Sinv, ntot);

    RmgFreeHost(mat_tem);

};

template <> double Exx_on<std::complex<double>>::Exxenergy(std::complex<double> * )
{return 0.0;};
template <> double Exx_on<double>::Exxenergy(double *rho_mat_global)
{
    double exx = 0.0;
    for(int idx = 0; idx < Phi.num_tot * Phi.num_tot; idx++) exx += ct.exx_fraction * rho_mat_global[idx] * this->Xij_mat[idx];
    // for spin polarized, rho_mat_global use the true occupation, for spin non polarized, we use occ[] *0.5 
    // so the energy should be scale correctly.
    MPI_Allreduce(MPI_IN_PLACE, &exx, 1, MPI_DOUBLE, MPI_SUM, pct.spin_comm);
    if(ct.nspin == 2) exx *= 0.5;
    return exx;
}
template <> void Exx_on<std::complex<double>>::HijExx(std::complex<double> *, LocalObject<std::complex<double>> &Phi){};
template <> void Exx_on<double>::HijExx(double *Hij_glob, LocalObject<double> &Phi)
{
    // calculte the overlap <PreOrbital| current Phi>

    double t1 = Rmg_G->get_NX_GRID(Phi.density);
    t1 *= Rmg_G->get_NY_GRID(Phi.density);
    t1 *= Rmg_G->get_NZ_GRID(Phi.density);

    double vol = Rmg_L.get_omega() /t1;

    int na = Phi.num_thispe ;
    int ntot = Phi.num_tot ;

    int pbasis = Phi.pbasis;

    double zero = 0.0, one = 1.0;

    double *mat_glob = (double *)RmgMallocHost(ntot*ntot*sizeof(double));
    double *mat_tem = (double *)RmgMallocHost(ntot*ntot*sizeof(double));
    RmgGemm("C", "N", na, na, pbasis, vol, PreOrbital, pbasis,
            Phi.storage_cpu, pbasis, zero, mat_tem, na);

    mat_local_to_glob(mat_tem, mat_glob, Phi, Phi, 0, ntot, 0, ntot, true);

    RmgGemm("N", "N", ntot, ntot, ntot, one, Xij_mat_Sinv, ntot,
            mat_glob, ntot, zero, mat_tem, ntot);

    RmgGemm("T", "N", ntot, ntot, ntot, ct.exx_fraction, mat_glob, ntot,
            mat_tem, ntot, one, Hij_glob, ntot);


    RmgFreeHost(mat_tem);
    RmgFreeHost(mat_glob);

}


template <> void Exx_on<std::complex<double>>::OmegaSinv(std::complex<double> *, LocalObject<std::complex<double>> &Phi){};
template <> void Exx_on<double>::OmegaSinv(double *Sij_inverse, LocalObject<double> &Phi)
{
    int ntot = Phi.num_thispe ;

    int pbasis = Phi.pbasis;

    double zero = 0.0, one = 1.0;

    double *phi_tem = (double *)RmgMallocHost(ntot*pbasis*sizeof(double));

    RmgGemm("N", "N", pbasis, ntot, ntot, one, Omega_j, pbasis,
            Sij_inverse, ntot, zero, phi_tem, pbasis);

    size_t size = ntot * pbasis * sizeof(double);
    memcpy(Omega_j, phi_tem, size);

    RmgFreeHost(phi_tem);

}


template <> void Exx_on<std::complex<double>>::OmegaRes(std::complex<double> *, LocalObject<std::complex<double>> &Phi){};
template <> void Exx_on<double>::OmegaRes(double *res, LocalObject<double> &Phi)
{
    double t1 = Rmg_G->get_NX_GRID(Phi.density);
    t1 *= Rmg_G->get_NY_GRID(Phi.density);
    t1 *= Rmg_G->get_NZ_GRID(Phi.density);

    double vol = Rmg_L.get_omega() /t1;

    int na = Phi.num_thispe ;
    int ntot = Phi.num_tot ;

    int pbasis = Phi.pbasis;

    double zero = 0.0, one = 1.0;

    double *mat_glob = (double *)RmgMallocHost(ntot*ntot*sizeof(double));
    double *mat_local = (double *)RmgMallocHost(ntot*na*sizeof(double));
    RmgGemm("C", "N", na, na, pbasis, vol, PreOrbital, pbasis,
            Phi.storage_cpu, pbasis, zero, mat_local, na);

    mat_local_to_glob(mat_local, mat_glob, Phi, Phi, 0, ntot, 0, ntot, true);
    for(int st1 = 0; st1 < ntot; st1++)
    for(int st2 = 0; st2 < na; st2++)
    {
        int st2_glob = Phi.index_proj_to_global[st2];
        mat_local[st2 * ntot + st1] = mat_glob[st2_glob * ntot + st1];
    }

    RmgGemm("N", "N", pbasis, na, ntot, ct.exx_fraction, Omega_j, pbasis,
            mat_local, ntot, one, res, pbasis);

    RmgFreeHost(mat_local);
    RmgFreeHost(mat_glob);

}

template <> void Exx_on<double>::Omega_rmg(double *Cij_local, double *Cij_global, bool use_float_fft)
{
    double one = 1.0, zero = 0.0;
    int num_tot = Phi.num_tot;
    int num_thispe = Phi.num_thispe;
    int pbasis = Phi.pbasis;
    RmgTimer RT0("5-Functional: Exx potential");

    // DePhi: waveufnctions: Phi * C
    RmgGemm("N", "N", pbasis, num_tot, num_thispe, one, Phi.storage_cpu, pbasis,
            Cij_local, num_thispe, zero, DePhi, pbasis);


    // Clear vexx
    for(int idx=0;idx < Phi.num_tot*pbasis;idx++) Omega_j[idx] = 0.0;
    this->Exxb->Vexx(Omega_j, true);
    //this->Exxb->Vexx(Omega_j, use_float_fft);

    RmgGemm("T", "N", num_tot, num_tot, pbasis, one, Omega_j, pbasis,
            DePhi, pbasis, zero, Cij_local, num_tot);
    size_t size = num_tot * num_tot;
    MPI_Allreduce(MPI_IN_PLACE, Cij_local, size, MPI_DOUBLE, MPI_SUM, Phi.comm);
    

    int info;
    int lwork = num_tot * num_tot;
    int *ipiv = new int[num_tot];
    double *work = new double[lwork];

    dgetrf(&num_tot, &num_tot, Cij_global, &num_tot, ipiv, &info);
    if (info != 0)
    {
        rmg_printf ("error in dgetrf in Exx_on::Omega_rmg  with INFO = %d \n", info);
        fflush (NULL);
        exit (0);
    }
    dgetri(&num_tot, Cij_global, &num_tot, ipiv, (double *)work, &lwork, &info);
    if (info != 0)
    {
        rmg_printf ("error in dgetri in Exx_on::Omega_rmg  with INFO = %d \n", info);
        fflush (NULL);
        exit (0);
    }

    RmgGemm("N", "N", pbasis, num_tot, num_tot, one, Omega_j, pbasis,
            Cij_global, num_tot, zero, DePhi, pbasis);


    size = num_tot * pbasis * sizeof(double);
    memcpy(Omega_j, DePhi, size);
    
    delete [] ipiv; 
    delete [] work;



}

template <> void Exx_on<std::complex<double>>::Omega_rmg(std::complex<double> *rho_matrix, std::complex<double> *cij_glob,  bool use_fft_float)
{
}
