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
#include "blas.h"
#include "Wannier.h"


void transpose(std::complex<double> *m, int w, int h);
template Wannier<double>::Wannier(BaseGrid &, Lattice &, const std::string &, int, int, int, double, double, double *);
template Wannier<std::complex<double>>::Wannier(BaseGrid &, Lattice &, const std::string &, int, int,int, double, double, std::complex<double> *);

template Wannier<double>::~Wannier(void);
template Wannier<std::complex<double>>::~Wannier(void);
template <class T> Wannier<T>::~Wannier ()
{
}

template <class T> Wannier<T>::Wannier (
        BaseGrid &G_in,
        Lattice &L_in,
        const std::string &wavefile_in,
        int nstates_in,
        int nwannier_in,
        int scdm_in,
        double scdm_mu_in,
        double scdm_sigma_in,
        T *psi_in) : 
    G(G_in), L(L_in), wavefile(wavefile_in), nstates(nstates_in), n_wannier(nwannier_in), scdm(scdm_in), 
    scdm_mu(scdm_mu_in), scdm_sigma(scdm_sigma_in), psi(psi_in)
{
    RmgTimer RT0("5-Wannier: total");
    if(ct.kpoint_mesh[0] <= 0 || ct.kpoint_mesh[1] <= 0 || ct.kpoint_mesh[2] <=0)
    {
        throw RmgFatalException() << "kpoint mesh must be set up  \n";
    }
    if(ct.kpoint_is_shift[0] != 0 || ct.kpoint_is_shift[1] != 0 || ct.kpoint_is_shift[2] !=0)
    {
        throw RmgFatalException() << "kpoint must include gamma point, kpoint_is_shift=0, 0, 0  \n";
    }

    ngrid = G.get_NX_GRID(1) * G.get_NY_GRID(1) * G.get_NZ_GRID(1);
    ngrid_noncoll = ngrid * ct.noncoll_factor;

    if(scdm == ISOLATED_ENTANGLEMENT) nstates = n_wannier;
    std::vector<double> occs;
    occs.resize(nstates, 1.0);
    Exx = new Exxbase<T>(G, G, L, wavefile, nstates, occs.data(), psi, EXX_LOCAL_FFT);
    RmgTimer *RT1 = new RmgTimer("5-Wannier: writesingle file");
    Exx->WriteWfsToSingleFile();
    MPI_Barrier(MPI_COMM_WORLD);
    WriteWinEig();
    delete RT1;
}

template <> void Wannier<double>::SetAmn()
{
        throw RmgFatalException() << "scdm need more than one gamma point  \n";
}
template <> void Wannier<std::complex<double>>::SetAmn()
{
    double tol = 1.0e-5;
    int ik_gamma = -1;
    
    for(int kpt = 0; kpt < ct.num_kpts; kpt++)
    {
        if( std::abs(ct.kp[kpt].kpt[0]) + std::abs(ct.kp[kpt].kpt[0]) + std::abs(ct.kp[kpt].kpt[0]) < tol)
        {
            ik_gamma = kpt;
            break;
        }
 
    }
    
    if(ik_gamma < 0) 
        throw RmgFatalException() << "cannot find gamma kpoint \n";

    int *piv = new int[ngrid_noncoll]();
    //size_t length = (size_t)nstates * (size_t)ngrid_noncoll * sizeof(std::complex<double>);
    size_t length = (size_t)nstates * (size_t)ngrid_noncoll;
    psi_s = new std::complex<double>[length];
    if(pct.gridpe == 0)
    {
        std::string filename = wavefile + "_spin"+std::to_string(pct.spinpe) + "_kpt" + std::to_string(ik_gamma);
        RmgTimer *RT1 = new RmgTimer("5-Wannier: mmap");
        serial_fd = open(filename.c_str(), O_RDONLY, (mode_t)0600);
        if(serial_fd < 0)
            throw RmgFatalException() << "Error! Could not open " << filename << " . Wannier Terminating.\n";
        delete RT1;

        read(serial_fd, psi_s, length * sizeof(std::complex<double>));
     //   psi_s = (std::complex<double> *)mmap(NULL, length, PROT_READ, MAP_PRIVATE, serial_fd, 0);

        for(int st = 0; st < nstates; st++)
        {
            
            double focc = 1.0;
            //double eigs = ct.kp[ik_gamma].eigs[st];
            //double tem = (eigs - scdm_mu)/scdm_sigma;
            double eigs = 0.0;
            double tem = (eigs - scdm_mu)/scdm_sigma;
            if(scdm == ISOLATED_ENTANGLEMENT) 
            {
                std::cout<< scdm << " sss " << std::endl;
                focc = 1.0;
            }
            else if(scdm == GAU_ENTANGLEMENT) focc = std::exp(-tem * tem);
            else if(scdm == ERFC_ENTANGLEMENT) focc = 0.5 * erfc(tem);
            else
                throw RmgFatalException() << "scdm = " << scdm << ISOLATED_ENTANGLEMENT<< "  wrong value \n";
            
            for(int idx = 0; idx < ngrid_noncoll; idx++)
            {
                psi_s[st * ngrid_noncoll + idx] *= focc;
                printf("\n aaa %d %f %f  %f %f", idx, psi_s[st * ngrid_noncoll + idx], psi[st * ngrid_noncoll + idx]);
            }
        }

        transpose(psi_s, ngrid_noncoll, nstates);


        std::complex<double> *tau = new std::complex<double>[2*nstates];
        double *rwork = new double[2 * ngrid_noncoll];
        int Lwork = -1, info;
        std::complex<double> Lwork_tmp;
        zgeqp3(&nstates, &ngrid_noncoll, psi_s, &nstates, piv, tau, &Lwork_tmp, &Lwork, rwork,&info);

        Lwork = (int)(std::real(Lwork_tmp)) + 1;

        std::complex<double> *cwork = new std::complex<double>[Lwork];

        zgeqp3(&nstates, &ngrid_noncoll, psi_s, &nstates, piv, tau, cwork, &Lwork, rwork,&info);

        if(info != 0) throw RmgFatalException() << "Error! in zgeqp3 at" << __FILE__ << __LINE__ << "  Wannier Terminating.\n";
        for(int i = 0; i < n_wannier; i++) printf("\n piv %d   %d ", piv[i], info);
        delete [] cwork;
        delete [] rwork;
        delete [] tau;
    //    munmap(psi_s, length);
        close(serial_fd);
    }
    MPI_Bcast(piv, ngrid_noncoll, MPI_INT, 0, pct.grid_comm);

    std::cout << " bbb " << std::endl;
    int nx_grid = G.get_NX_GRID(1);
    int ny_grid = G.get_NY_GRID(1);
    int nz_grid = G.get_NZ_GRID(1);

    int num_q = ct.klist.num_k_all;
    std::complex<double> *psi_wan = new std::complex<double>[n_wannier * nstates];
    std::complex<double> *Amn = new std::complex<double>[num_q * n_wannier * nstates]();

    for(int iq = pct.gridpe; iq < num_q; iq+= pct.grid_npes)
    {
        int ik = ct.klist.k_map_index[iq];
        int isym = ct.klist.k_map_symm[iq];
        int isyma = std::abs(isym);

        int sym_rot_inv[9], ftau_inv[3], det, sr[3][3];

    std::cout << " iq= " << iq << std::endl;
        for(int i= 0; i < 3; i++)
        {
            for(int j= 0; j < 3; j++)
            {
                sr[i][j] = Rmg_Symm->sym_rotate[isyma *9 +i*3 + j];
            }
        }
        det = sr[0][0] * ( sr[1][1] * sr[2][2] - sr[2][1] * sr[1][2] ) -
            sr[0][1] * ( sr[1][0] * sr[2][2] - sr[2][0] * sr[1][2] ) +
            sr[0][2] * ( sr[1][0] * sr[2][1] - sr[2][0] * sr[1][1] );

        for(int i= 0; i < 3; i++)
        {
            for(int j= 0; j < 3; j++)
            {
                sym_rot_inv[i*3+j] = sr[(i+1)%3 ][(j+1)%3] * sr[(i+2)%3 ][(j+2)%3];
                sym_rot_inv[i*3+j] -= sr[(i+1)%3 ][(j+2)%3] * sr[(i+2)%3 ][(j+1)%3];
                sym_rot_inv[i*3+j] /= det;

            }
        }

        for(int i= 0; i < 3; i++)
        {
            ftau_inv[i] = 0;
            for(int j= 0; j < 3; j++)
            {
                ftau_inv[i] -= sym_rot_inv[i*3+j] * Rmg_Symm->ftau_wave[isyma*3 + j];
            }
        }
        std::string filename = wavefile + "_spin"+std::to_string(pct.spinpe) + "_kpt" + std::to_string(ik);
        RmgTimer *RT1 = new RmgTimer("5-Wannier: mmap");
        serial_fd = open(filename.c_str(), O_RDONLY, (mode_t)0600);
        if(serial_fd < 0)
            throw RmgFatalException() << "Error! Could not open " << filename << " . Wannier Terminating.\n";
        delete RT1;

        //psi_s = (std::complex<double> *)mmap(NULL, length, PROT_READ, MAP_PRIVATE, serial_fd, 0);
        read(serial_fd, psi_s, length * sizeof(std::complex<double>));

        int ix, iy, iz, ixx, iyy, izz;
        for(int iw = 0; iw < n_wannier; iw++)
        {
            int grid_idx = piv[iw]-1;
            int ispin = 0, ispin_rev=0;
            if(grid_idx >= ngrid) 
            {
                ispin = 1;
                grid_idx -= ngrid;
            }

            if(ct.noncoll) ispin_rev = (ispin + 1) %2;
            ix = grid_idx/ny_grid/nz_grid;
            iy = (grid_idx/nz_grid ) % ny_grid;
            iz = grid_idx % nz_grid;
            symm_ijk(sym_rot_inv, ftau_inv, ix, iy, iz, ixx, iyy, izz, nx_grid, ny_grid, nz_grid);

            printf("\n ddd %d %d %d %d %d %d   %d", ix, iy, iz, ixx, iyy, izz, grid_idx);
            for(int st = 0; st < nstates; st++)
            {
                double focc = 1.0;
                double eigs = ct.kp[ik_gamma].eigs[st];
                double tem = (eigs - scdm_mu)/scdm_sigma;
                if(scdm == ISOLATED_ENTANGLEMENT ) focc = 1.0;
                else if(scdm == GAU_ENTANGLEMENT) focc = std::exp(-tem * tem);
                else if(scdm == ERFC_ENTANGLEMENT) focc = 0.5 * erfc(tem);
                else
                    throw RmgFatalException() << "scdm = " << scdm << "  wrong value \n";
                if(isym >= 0) 
                {
                    psi_wan[iw * nstates + st] = psi_s[st * ngrid_noncoll + ispin * ngrid 
                        + ixx * ny_grid * nz_grid + iyy * nz_grid + izz] * focc;

                }
                else
                {
                    psi_wan[iw * nstates + st] = std::conj(psi_s[st * ngrid_noncoll + ispin_rev * ngrid 
                            + ixx * ny_grid * nz_grid + iyy * nz_grid + izz])* focc;
                }
            }

        }

        double *sigma = new double[n_wannier];
        double *rwork = new double[5*n_wannier];
        std::complex<double> *Umat = new std::complex<double>[nstates * n_wannier];
        std::complex<double> *VTmat = new std::complex<double>[n_wannier * n_wannier];

        int info;
        int Lwork = -1;
        std::complex<double> Lwork_tmp;
        zgesvd("S", "S", &nstates, &n_wannier, psi_wan, &nstates,
           sigma, Umat, &nstates, VTmat,&n_wannier,&Lwork_tmp, &Lwork, rwork, &info);
        Lwork = (int)(std::real(Lwork_tmp)) + 1;
        std::complex<double> *work = new std::complex<double>[Lwork];

        zgesvd("S", "S", &nstates, &n_wannier, psi_wan, &nstates,
           sigma, Umat, &nstates, VTmat,&n_wannier, work, &Lwork, rwork, &info);
        if(info != 0) throw RmgFatalException() << "Error! in zgesvd at" << __FILE__ << __LINE__ << "  Wannier Terminating.\n";

        std::complex<double> one(1.0), zero(0.0);
        zgemm("N", "N", &nstates, &n_wannier, &n_wannier, &one, Umat, &nstates, VTmat, &n_wannier, &zero, &Amn[iq*nstates*n_wannier], &nstates);

        delete [] sigma;
        delete [] rwork;
        delete [] work;
        delete [] Umat;
        delete [] VTmat;

       // munmap(psi_s, length);
        close(serial_fd);
    }

    
    int count = num_q * nstates * n_wannier;
    MPI_Allreduce(MPI_IN_PLACE, Amn, count, MPI_DOUBLE, MPI_SUM, pct.grid_comm);

    if(pct.imgpe == 0)
    {
        printf("\n Amn done");
        mkdir ("Wannier90_rmg", S_IRWXU);
        time_t tt;

        char *timeptr;
        time (&tt);
        timeptr = ctime (&tt);
        FILE *famn = fopen("Wannier90_rmg/wannier90.amn", "w+");
        fprintf(famn, "Created with scdm on %s", timeptr);
        fprintf(famn, "%8d  %8d  %8d    %10.6f  %10.6f ", nstates, num_q, n_wannier, scdm_mu, scdm_sigma);
        for(int iq = 0; iq <num_q; iq++)
        {
            for(int iw = 0; iw < n_wannier; iw++)
            {
                for(int st = 0; st < nstates; st++)
                {
                    fprintf(famn, "\n  %5d %5d %5d %18.12f %18.12f", st+1, iw+1, iq+1, std::real(Amn[iq * nstates * n_wannier + iw * nstates + st]),
                            std::imag(Amn[iq * nstates * n_wannier + iw * nstates + st]));
                }
            }
        }
    }


}

void transpose(std::complex<double> *m, int w, int h)
{
    int start, next, i;
    std::complex<double> tmp;

    for (start = 0; start <= w * h - 1; start++) {
        next = start;
        i = 0;
        do {	i++;
            next = (next % h) * w + next / h;
        } while (next > start);
        if (next < start || i == 1) continue;

        tmp = m[next = start];
        do {
            i = (next % h) * w + next / h;
            m[next] = (i == start) ? tmp : m[i];
            next = i;
        } while (next > start);
    }
}
template <> void Wannier<double>::SetMmn()
{
        throw RmgFatalException() << "scdm need more than one gamma point  \n";
}
template <> void Wannier<std::complex<double>>::SetMmn()
{

    double tol = 1.0e-5;
    int num_q = ct.klist.num_k_all;
    std::vector<std::vector<int>> kneighbor;
    kneighbor.resize(num_q);

    // find neighboring shell of kpoints
    // since the kpoint must be in a uniform grid, each kpoint has the same number of neighbors.
    std::vector<double> kdist_shell;
    double kd[3];
    for(int iq = 1; iq < num_q; iq++)
    {

       //kd[0] =  

    }
    
    
 
}

template void Wannier<double> ::WriteWinEig();
template void Wannier<std::complex<double>> ::WriteWinEig();
template <class T> void Wannier<T>::WriteWinEig()
{
    if(pct.imgpe == 0)
    {
        mkdir ("Wannier90_rmg", S_IRWXU);
        FILE *fwin = fopen("Wannier90_rmg/wannier90.win", "w+");
        fprintf(fwin, "write_hr        = true");
        fprintf(fwin, "\nnum_wann        =  %d", n_wannier);
        fprintf(fwin, "\nnum_iter        = 20");
        fprintf(fwin, "\nauto_projections= true\n");
        fprintf(fwin, "\nbegin atoms_cart");
        fprintf(fwin, "\nbohr");
        for(size_t ion = 0; ion < Atoms.size(); ion++)
        {

            ION &Atom = Atoms[ion];
            SPECIES &Type = Species[Atom.species];

            fprintf (fwin, "\n%-2s   %12.6e  %12.6e  %12.6e",
                    Type.atomic_symbol, Atom.crds[0], Atom.crds[1], Atom.crds[2]);

        }         
        fprintf(fwin, "\nend atoms_cart\n");
        fprintf(fwin, "\nbegin unit_cell_cart");
        fprintf(fwin, "\nbohr");
        fprintf(fwin, "\n%10.6f   %10.6f   %10.6f", L.get_a0(0),L.get_a0(1),L.get_a0(2));
        fprintf(fwin, "\n%10.6f   %10.6f   %10.6f", L.get_a1(0),L.get_a1(1),L.get_a1(2));
        fprintf(fwin, "\n%10.6f   %10.6f   %10.6f", L.get_a2(0),L.get_a2(1),L.get_a2(2));
        fprintf(fwin, "\nend_unit_cell_cart");
        fprintf(fwin, "\n\nmp_grid : %d %d %d", ct.kpoint_mesh[0], ct.kpoint_mesh[1], ct.kpoint_mesh[2]);

        fprintf(fwin, "\n\nbegin kpoints");
        for(int iq = 0; iq <ct.klist.num_k_all; iq++)
        {
            fprintf(fwin, "\n%10.6f   %10.6f   %10.6f", ct.klist.k_all_cart[iq][0], ct.klist.k_all_cart[iq][1], ct.klist.k_all_cart[iq][2] );
        }
        fprintf(fwin, "\nend kpoints");
            
        fclose(fwin);

        FILE *feig= fopen("Wannier90_rmg/wannier90.eig", "w+");
        for(int iq = 0; iq <ct.klist.num_k_all; iq++)
        {
            int ik = ct.klist.k_map_index[iq];
            for(int st = 0; st < nstates; st++)
            {
                fprintf(feig, "%5d %5d  %18.12f\n", st+1, iq+1, ct.kp[ik].eigs[st] *Ha_eV);
            }
        }
        fclose(feig);

    }
}
