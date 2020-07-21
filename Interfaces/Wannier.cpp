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
#include <iostream>
#include <string>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>



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
template Wannier<double>::Wannier(BaseGrid &, Lattice &, const std::string &, int, int, int, double, double, double *, Kpoint<double> **Kptr);
template Wannier<std::complex<double>>::Wannier(BaseGrid &, Lattice &, const std::string &, int, int,int, double, double, std::complex<double>
*, Kpoint<std::complex<double>> **Kptr);

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
        T *psi_in, Kpoint<T> **Kptr) : 
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


    if(!ct.norm_conserving_pp && ct.localize_projectors)
    {
        throw RmgFatalException() << "for ultra soft pseudopotential, set localize_projectors to be false for wannier90 interface  \n";
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
    Read_nnkpts();
    SetAmn();
    SetMmn(Kptr);
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
            
            double vel = L.get_omega() / ((double)ngrid);
            double norm_coeff = 0.0;
            for(int idx = 0; idx < ngrid_noncoll; idx++)
            {
                norm_coeff += std::norm(psi_s[st * ngrid_noncoll + idx]);
            }
            norm_coeff = std::sqrt(norm_coeff * vel);
            for(int idx = 0; idx < ngrid_noncoll; idx++)
            {
                psi_s[st * ngrid_noncoll + idx] *= norm_coeff;
            }

            double focc = 1.0;
            //double eigs = ct.kp[ik_gamma].eigs[st];
            //double tem = (eigs - scdm_mu)/scdm_sigma;
            double eigs = 0.0;
            double tem = (eigs - scdm_mu)/scdm_sigma;
            if(scdm == ISOLATED_ENTANGLEMENT) 
            {
                focc = 1.0;
            }
            else if(scdm == GAU_ENTANGLEMENT) focc = std::exp(-tem * tem);
            else if(scdm == ERFC_ENTANGLEMENT) focc = 0.5 * erfc(tem);
            else
                throw RmgFatalException() << "scdm = " << scdm << ISOLATED_ENTANGLEMENT<< "  wrong value \n";
            
            for(int idx = 0; idx < ngrid_noncoll; idx++)
            {
                psi_s[st * ngrid_noncoll + idx] = focc * (psi_s[st * ngrid_noncoll + idx]);
            }
        }

        transpose(psi_s, ngrid_noncoll, nstates);

        std::complex<double> *tau = new std::complex<double>[2*nstates]();
        double *rwork = new double[2 * ngrid_noncoll]();
        int Lwork = -1, info;
        std::complex<double> Lwork_tmp;
        zgeqp3(&nstates, &ngrid_noncoll, psi_s, &nstates, piv, tau, &Lwork_tmp, &Lwork, rwork,&info);

        Lwork = (int)(std::real(Lwork_tmp)) + 1;

        std::complex<double> *cwork = new std::complex<double>[Lwork]();

        zgeqp3(&nstates, &ngrid_noncoll, psi_s, &nstates, piv, tau, cwork, &Lwork, rwork,&info);

        if(info != 0) throw RmgFatalException() << "Error! in zgeqp3 at" << __FILE__ << __LINE__ << "  Wannier Terminating.\n";
        //for(int i = 0; i < n_wannier; i++) printf("\n piv %d   %d ", piv[i], info);
        delete [] cwork;
        delete [] rwork;
        delete [] tau;
    //    munmap(psi_s, length);
        close(serial_fd);
    }
    MPI_Bcast(piv, ngrid_noncoll, MPI_INT, 0, pct.grid_comm);

    //piv[0] = 437;
    //piv[1] = 3685;
    //piv[2] = 722;
    //piv[3] = 823;

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
        int isyma = std::abs(isym)-1;

        ReadRotatePsi(ik, isym, isyma, wavefile, psi_s);
        for(int st = 0; st < nstates; st++)
        {
            double vel = L.get_omega() / ((double)ngrid);
            double norm_coeff = 0.0;
            for(int idx = 0; idx < ngrid_noncoll; idx++)
            {
                norm_coeff += std::norm(psi_s[st * ngrid_noncoll + idx]);
            }
            norm_coeff = std::sqrt(norm_coeff * vel);
            for(int idx = 0; idx < ngrid_noncoll; idx++)
            {
                psi_s[st * ngrid_noncoll + idx] *= norm_coeff;
            }
        }

        int ix, iy, iz, ixx, iyy, izz;
        std::complex<double> phase_center;
        for(int iw = 0; iw < n_wannier; iw++)
        {
            int grid_idx = piv[iw]-1;
            int ispin = 0;
            if(grid_idx >= ngrid) 
            {
                ispin = 1;
                grid_idx -= ngrid;
            }

            ix = grid_idx/ny_grid/nz_grid;
            iy = (grid_idx/nz_grid ) % ny_grid;
            iz = grid_idx % nz_grid;

            ixx = ix;
            if(ixx > nx_grid/2) ixx -= nx_grid;
            iyy = iy;
            if(iyy > ny_grid/2) iyy -= ny_grid;
            izz = iz;
            if(izz > nz_grid/2) izz -= nz_grid;
            double kr_center =
                ct.klist.k_all_xtal[iq][0] * ixx/(double)nx_grid +
                ct.klist.k_all_xtal[iq][1] * iyy/(double)ny_grid +
                ct.klist.k_all_xtal[iq][2] * izz/(double)nz_grid ;
            phase_center = std::exp(std::complex<double>(0.0, -kr_center * twoPI));


            for(int st = 0; st < nstates; st++)
            {
                double focc = 1.0;
                double eigs = ct.kp[ik].eigs[st];
                double tem = (eigs - scdm_mu)/scdm_sigma;
                if(scdm == ISOLATED_ENTANGLEMENT ) focc = 1.0;
                else if(scdm == GAU_ENTANGLEMENT) focc = std::exp(-tem * tem);
                else if(scdm == ERFC_ENTANGLEMENT) focc = 0.5 * erfc(tem);
                else
                {
                    throw RmgFatalException() << "scdm = " << scdm << "  wrong value \n";
                }
                psi_wan[iw * nstates + st] = std::conj(psi_s[st * ngrid_noncoll + ispin * ngrid 
                        + ix * ny_grid * nz_grid + iy * nz_grid + iz]) * focc * phase_center;
            }
        }


        double *sigma = new double[n_wannier]();
        double *rwork = new double[5*n_wannier]();
        std::complex<double> *Umat = new std::complex<double>[nstates * n_wannier]();
        std::complex<double> *VTmat = new std::complex<double>[n_wannier * n_wannier]();

        int info;
        int Lwork = -1;
        std::complex<double> Lwork_tmp;
        zgesvd("S", "S", &nstates, &n_wannier, psi_wan, &nstates,
                sigma, Umat, &nstates, VTmat,&n_wannier,&Lwork_tmp, &Lwork, rwork, &info);
        Lwork = (int)(std::real(Lwork_tmp)) + 1;
        std::complex<double> *work = new std::complex<double>[Lwork]();

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
    MPI_Allreduce(MPI_IN_PLACE, Amn, count, MPI_DOUBLE_COMPLEX, MPI_SUM, pct.grid_comm);

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
        fprintf(fwin, "\nbands_plot       = true\n");

        fprintf(fwin, "bands_plot_format= xmgrace\n");
        fprintf(fwin, "bands_num_points: 100\n");
        fprintf(fwin, "begin kpoint_path\n");
        fprintf(fwin, "   G 0.0 0.0 0.0 L 0.0 0.0 1.0\n");
        fprintf(fwin, "  L 0.0 0.0 1.0 N 0.0 1.0 1.0\n");
        fprintf(fwin, "end kpoint_path\n");

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
            fprintf(fwin, "\n%10.6f   %10.6f   %10.6f", ct.klist.k_all_xtal[iq][0], ct.klist.k_all_xtal[iq][1], ct.klist.k_all_xtal[iq][2] );
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


template void Wannier<double> ::Read_nnkpts();
template void Wannier<std::complex<double>> ::Read_nnkpts();
template <class T> void Wannier<T>::Read_nnkpts()
{
    std::ifstream fnnk("Wannier90_rmg/wannier90.nnkp");
    if(!fnnk.is_open())
        throw RmgFatalException() << "generte wannier90.nnkpts first by running wannier90.x -pp \n";

    std::string oneline;
    // find the first nnkpts 
    while(getline(fnnk, oneline) )
    {
        if(oneline.find("nnkpts") != std::string::npos ) break;
    }
    getline(fnnk, oneline);

    boost::trim(oneline);
    ct.klist.num_k_nn = std::stoi(oneline);

    ct.klist.k_neighbors.resize(boost::extents[ct.klist.num_k_all][ct.klist.num_k_nn][4]);

    std::vector<std::string> kn_info;
    if(pct.imgpe == 0 && ct.verbose) 
        printf("\n kpoint neighbors from nnkpts \n");

    std::string whitespace_delims = " \n\t";
    for(int i = 0; i < ct.klist.num_k_all; i++)
    {
        for(int j = 0; j < ct.klist.num_k_nn; j++)
        {
            getline(fnnk, oneline);

            boost::trim(oneline);
            boost::algorithm::split( kn_info, oneline, boost::is_any_of(whitespace_delims), boost::token_compress_on );
            if(std::stoi(kn_info[0]) != i+1 || kn_info.size() != 5) 
            {
                std::cout << " nnkpts problem i+1 = " << i+1 << "  ik = " << kn_info[0] << std::endl;
                throw RmgFatalException() << "nnkpts in wannier90.nnkpts file has problem  \n";
            }
            ct.klist.k_neighbors[i][j][0] = std::stoi(kn_info[1]) -1;
            ct.klist.k_neighbors[i][j][1] = std::stoi(kn_info[2]);
            ct.klist.k_neighbors[i][j][2] = std::stoi(kn_info[3]);
            ct.klist.k_neighbors[i][j][3] = std::stoi(kn_info[4]);
            if(pct.imgpe == 0 && ct.verbose) 
                printf("%5d  %5d  %5d  %5d  %5d \n", std::stoi(kn_info[0])-1,std::stoi(kn_info[1])-1,
                        std::stoi(kn_info[2]),std::stoi(kn_info[3]),std::stoi(kn_info[4]));


        }
    }

}
template  void Wannier<double>::ReadRotatePsi(int ikindex, int isym, int isyma, std::string wavefile, double *psi_k);
template  void Wannier<std::complex<double>>::ReadRotatePsi(int ikindex, int isym, int isyma, std::string wavefile, std::complex<double> *psi_k);
template <class T> void Wannier<T>::ReadRotatePsi(int ikindex, int isym, int isyma, std::string wavefile, T *psi_k)
{

    int nx_grid = G.get_NX_GRID(1);
    int ny_grid = G.get_NY_GRID(1);
    int nz_grid = G.get_NZ_GRID(1);
    int nbasis = nx_grid * ny_grid * nz_grid;
    int nbasis_noncoll = nx_grid * ny_grid * nz_grid * ct.noncoll_factor;
    size_t length = nstates * nbasis_noncoll * sizeof(T);
    T *psi_map;

    std::string filename = wavefile + "_spin"+std::to_string(pct.spinpe) + "_kpt" + std::to_string(ikindex);
    int serial_fd = open(filename.c_str(), O_RDONLY, (mode_t)0600);
    if(serial_fd < 0)
        throw RmgFatalException() << "Error! Could not open " << filename << " . Terminating.\n";

    psi_map = (T *)mmap(NULL, length, PROT_READ, MAP_PRIVATE, serial_fd, 0);


    int ixx, iyy, izz;
    for (int ix = 0; ix < nx_grid; ix++) {
        for (int iy = 0; iy < ny_grid; iy++) {
            for (int iz = 0; iz < nz_grid; iz++) {

                symm_ijk(&Rmg_Symm->sym_rotate[isyma *9], &Rmg_Symm->ftau_wave[isyma*3], ix, iy, iz, ixx, iyy, izz, nx_grid, ny_grid, nz_grid);

                if(ct.noncoll) 
                    throw RmgFatalException() << "wave function rotation need to be done correctly  \n";

                if(isym >= 0)
                {
                    for(int st = 0; st < nstates; st++)
                    {
                        psi_k[st * nbasis + ix * ny_grid * nz_grid + iy * nz_grid + iz]
                            = (psi_map[st * nbasis + ixx * ny_grid * nz_grid + iyy * nz_grid + izz]);
                    }
                }
                else
                {
                    for(int st = 0; st < nstates; st++)
                    {
                        psi_k[st * nbasis + ix * ny_grid * nz_grid + iy * nz_grid + iz]
                            = MyConj(psi_map[st * nbasis + ixx * ny_grid * nz_grid + iyy * nz_grid + izz]);
                    }
                }
            }
        }
    }

    munmap(psi_map, length);
    close(serial_fd);

}


template  void Wannier<double>::SetMmn(Kpoint<double> **Kptr);
template  void Wannier<std::complex<double>>::SetMmn(Kpoint<std::complex<double>> **Kptr);
template <class T> void Wannier<T>::SetMmn(Kpoint<T> **Kptr)
{

    int num_q = ct.klist.num_k_all;
    int num_kn = ct.klist.num_k_nn;

    int nx_grid = G.get_NX_GRID(1);
    int ny_grid = G.get_NY_GRID(1);
    int nz_grid = G.get_NZ_GRID(1);
    int nbasis = nx_grid * ny_grid * nz_grid;
    int nbasis_noncoll = nbasis * ct.noncoll_factor;
    size_t length = nstates * nbasis_noncoll * sizeof(T);
    T *psi_k = (T *)GpuMallocManaged(length);
    T *psi_q = (T *)GpuMallocManaged(length);
    T *Mmn = (T *)GpuMallocManaged(num_q * num_kn * nstates * nstates * sizeof(T));

    double vel = L.get_omega() / ((double)nbasis);
    T alpha(vel);
    T beta = 0.0;

    T kphase(1.0);
    double kr;
    //double *kphase_R = (double *)&kphase;
    std::complex<double> *kphase_C = (std::complex<double> *)&kphase;
    for(int ik = 0; ik < num_q; ik++)
    {

        int ik_irr = ct.klist.k_map_index[ik];
        int isym = ct.klist.k_map_symm[ik];
        int isyma = std::abs(isym) -1;

        ReadRotatePsi(ik_irr, isym, isyma, wavefile, psi_k);

        for(int ikn = 0; ikn < num_kn; ikn++)
        {
            int ikn_index = ct.klist.k_neighbors[ik][ikn][0];
            int ikn_irr = ct.klist.k_map_index[ikn_index];
            int isym_kn = ct.klist.k_map_symm[ikn_index];
            int isyma_kn = std::abs(isym_kn) -1;
            ReadRotatePsi(ikn_irr, isym_kn, isyma_kn, wavefile, psi_q);


            //  phase factor when kneight need a periodic BZ folding
            for (int iz = 0; iz < nz_grid; iz++) {
                for (int iy = 0; iy < ny_grid; iy++) {
                    for (int ix = 0; ix < nx_grid; ix++) {
                        kr = ct.klist.k_neighbors[ik][ikn][1] * ix/(double)nx_grid
                            +ct.klist.k_neighbors[ik][ikn][2] * iy/(double)ny_grid
                            +ct.klist.k_neighbors[ik][ikn][3] * iz/(double)nz_grid;
                        *kphase_C = std::exp( std::complex<double>(0.0, -kr * twoPI));
                        for(int st = 0; st < nstates; st++)
                        {
                            psi_q[st * nbasis_noncoll + ix * ny_grid * nz_grid + iy * nz_grid + iz] *= kphase;
                            if(ct.noncoll)
                                psi_q[st * nbasis_noncoll + nbasis + ix * ny_grid * nz_grid + iy * nz_grid + iz] *= kphase;
                        }

                    }
                }
            }

            RmgGemm("C", "N", nstates, nstates, nbasis_noncoll, alpha, psi_k, nbasis_noncoll, psi_q, nbasis_noncoll,
                    beta, &Mmn[(ik*num_kn+ikn)*nstates*nstates], nstates);

            if(!ct.norm_conserving_pp || ct.noncoll)
                Mmn_us(ik_irr, isym, isyma, ikn_irr, isym_kn, isyma_kn, &Mmn[(ik*num_kn+ikn)*nstates*nstates], Kptr);

        }


    }


    if(pct.imgpe == 0)
    {
        time_t tt;

        char *timeptr;
        time (&tt);
        timeptr = ctime (&tt);
        FILE *fmmn = fopen("Wannier90_rmg/wannier90.mmn", "w+");
        fprintf(fmmn, "Created with RMG on %s", timeptr);
        fprintf(fmmn, "%8d  %8d  %8d    ", nstates, num_q, num_kn);
        for(int iq = 0; iq <num_q; iq++)
        {
            for(int ikn = 0; ikn < num_kn; ikn++)
            {
                fprintf(fmmn,"\n %5d  %5d  %5d  %5d  %5d", iq+1, ct.klist.k_neighbors[iq][ikn][0]+1,
                        ct.klist.k_neighbors[iq][ikn][1],ct.klist.k_neighbors[iq][ikn][2],ct.klist.k_neighbors[iq][ikn][3]);
                std::complex<double> *Mmn_C = (std::complex<double> *)&Mmn[(iq*num_kn+ikn)*nstates*nstates];
                for(int st1 = 0; st1 < nstates; st1++)
                    for(int st2 = 0; st2 < nstates; st2++)
                    {
                        fprintf(fmmn, "\n  %18.12f %18.12f", std::real(Mmn_C[ st1* nstates + st2]),
                                std::imag(Mmn_C[ st1* nstates + st2]));
                    }
            }
        }
        printf("\n Mmn done");
        fclose(fmmn);
    }


    GpuFreeManaged(psi_k);
    GpuFreeManaged(psi_q);
    GpuFreeManaged(Mmn);

}
template  void Wannier<double>::Mmn_us(int, int, int, int, int, int, double *, Kpoint<double> **Kptr);
template  void Wannier<std::complex<double>>::Mmn_us(int, int, int, int, int, int, std::complex<double> *, Kpoint<std::complex<double>> **Kptr);
template <class T> void Wannier<T>::Mmn_us(int ik_irr, int isym, int isyma, int kn_irr, int isym_kn, int isyma_kn, T *Mmn_onekpair, Kpoint<T> **Kptr)
{
    
    size_t num_nonloc_ions = Kptr[ik_irr]->BetaProjector->get_num_nonloc_ions();
    if(Atoms.size() != num_nonloc_ions)
    {
        throw RmgFatalException() << "set localize_projectors = false to use Wannier90 interface: Terminating.\n";
    }
    
    int *nonloc_ions_list = Kptr[ik_irr]->BetaProjector->get_nonloc_ions_list();
    int num_tot_proj = Kptr[ik_irr]->BetaProjector->get_num_tot_proj();
    int pstride = Kptr[ik_irr]->BetaProjector->get_pstride();

    size_t alloc = (size_t)num_tot_proj * (size_t)nstates * ct.noncoll_factor * sizeof(T);
    T *sint_ik_irr = (T *)GpuMallocManaged(alloc);
    T *sint_kn_irr = (T *)GpuMallocManaged(alloc);
    T *sint_ik = (T *)GpuMallocManaged(alloc);
    T *sint_kn = (T *)GpuMallocManaged(alloc);
    for(int i = 0; i < num_tot_proj * nstates * ct.noncoll_factor; i++) 
    {
        sint_ik[i] = 0.0;
        sint_kn[i] = 0.0;
    }

//  from <beta|psi> of irreducible k, rotate to get symmetry-related <beta|psu> |Sk 
    boost::multi_array_ref<T, 4> sint_ik_irr_4d{sint_ik_irr, boost::extents[nstates][ct.noncoll_factor][Atoms.size()][pstride]};
    boost::multi_array_ref<T, 4> sint_kn_irr_4d{sint_kn_irr, boost::extents[nstates][ct.noncoll_factor][Atoms.size()][pstride]};
    boost::multi_array_ref<T, 4> sint_ik_4d{sint_ik, boost::extents[nstates][ct.noncoll_factor][Atoms.size()][pstride]};
    boost::multi_array_ref<T, 4> sint_kn_4d{sint_kn, boost::extents[nstates][ct.noncoll_factor][Atoms.size()][pstride]};

    // Repack the sint array
    for(int istate = 0; istate < nstates; istate++)
    {
        size_t sindex = istate * num_nonloc_ions * pstride * ct.noncoll_factor;
        for(int ispin = 0; ispin < ct.noncoll_factor; ispin++)
        { 
            for (size_t ion = 0; ion < num_nonloc_ions; ion++)
            {
                int proj_index = ion * pstride;
                T *psint_ik = &Kptr[ik_irr]->newsint_local[proj_index + sindex + ispin * num_nonloc_ions* pstride];
                T *psint_kn = &Kptr[kn_irr]->newsint_local[proj_index + sindex + ispin * num_nonloc_ions* pstride];
                int gion = nonloc_ions_list[ion];
                for (int i = 0; i < pstride; i++)
                {
                    sint_ik_irr_4d[istate][ispin][gion][i] = psint_ik[i];
                    sint_kn_irr_4d[istate][ispin][gion][i] = psint_kn[i];

                }
            }
        }
    }



    for (int ion = 0; ion < ct.num_ions; ion++)
    {
        int ion_ik = Rmg_Symm->sym_atom[isyma * ct.num_ions + ion];
        int ion_kn = Rmg_Symm->sym_atom[isyma_kn * ct.num_ions + ion];
        SPECIES &AtomType = Species[Atoms[ion].species];

        for(int ih = 0; ih < AtomType.nh; ih++)
        {
            int l_ih = AtomType.nhtol[ih]; 
            int m_ih = AtomType.nhtom[ih];
            int b_ih = AtomType.indv[ih];
            for(int jh = 0; jh < AtomType.nh; jh++)
            {
                int m_jh = AtomType.nhtom[jh];
                int b_jh = AtomType.indv[jh];

                if(b_ih != b_jh) continue;

                for(int is = 0; is < ct.noncoll_factor; is++)
                {
                    int is_ik = is;
                    int is_kn = is;
                    if(Rmg_Symm->time_rev[isyma]) is_ik = (is +1 )%2;
                    if(Rmg_Symm->time_rev[isyma_kn]) is_kn = (is +1)%2;

                    for(int st = 0; st < nstates; st++)
                    {

                        sint_ik_4d[st][is_ik][ion][ih] += Rmg_Symm->rot_ylm[isyma][l_ih][m_ih][m_jh] * sint_ik_irr_4d[st][is][ion_ik][jh];
                        sint_kn_4d[st][is_kn][ion][ih] += Rmg_Symm->rot_ylm[isyma_kn][l_ih][m_ih][m_jh] * sint_kn_irr_4d[st][is][ion_kn][jh];
                    }
                }
            }
        }
    }

    T ZERO_t(0.0);
    T ONE_t(1.0);


    double *qqq;

    int M_cols = (size_t)num_tot_proj * ct.noncoll_factor;
    size_t alloc1 = M_cols * M_cols;
    T *M_qqq = (T *)GpuMallocManaged(sizeof(T) * alloc1);
    std::complex<double> *M_qqq_C = (std::complex<double> *) M_qqq;

    for (size_t i = 0; i < alloc1; i++)
    {
        M_qqq[i] = ZERO_t;
    }


    // set up M_qqq and M_dnm, this can be done outside in the
    // init.c or get_ddd get_qqq, we need to check the order
    int proj_index = 0;
    for (size_t ion = 0; ion < num_nonloc_ions; ion++)
    {

        /*Actual index of the ion under consideration*/
        proj_index = ion * ct.max_nl;
        int gion = nonloc_ions_list[ion];

        SPECIES &AtomType = Species[Atoms[gion].species];
        int nh = AtomType.nh;
        qqq = Atoms[gion].qqq;

        for (int i = 0; i < nh; i++)
        {
            int inh = i * nh;
            for (int j = 0; j < nh; j++)
            {

                if(ct.noncoll)
                {
                    int it0 = proj_index + i;
                    int jt0 = proj_index + j;
                    int it1 = proj_index + i + num_tot_proj;
                    int jt1 = proj_index + j + num_tot_proj;
                    M_qqq_C[it0 * num_tot_proj * 2 + jt0] = Atoms[gion].qqq_so[inh+j + 0 * nh *nh];
                    M_qqq_C[it0 * num_tot_proj * 2 + jt1] = Atoms[gion].qqq_so[inh+j + 1 * nh *nh];
                    M_qqq_C[it1 * num_tot_proj * 2 + jt0] = Atoms[gion].qqq_so[inh+j + 2 * nh *nh];
                    M_qqq_C[it1 * num_tot_proj * 2 + jt1] = Atoms[gion].qqq_so[inh+j + 3 * nh *nh];
                }
                else
                {
                    int idx = (proj_index + i) * num_tot_proj + proj_index + j;
                    M_qqq[idx] = (T)qqq[inh+j];
                }
            }

        }
    }



    int dim_dnm = num_tot_proj * ct.noncoll_factor;

    //M_dnm: dim_dnm * dim_dnm matrxi
    //sint_compack: dim_dnm * num_states == num_tot_proj * ct.noncoll_factor * num_states
    //nwork: dim_dnm * num_states == num_tot_proj * ct.noncoll_factor * num_states
    //  in the first RmgGemm, nwork is a matrix of (dim_dnm) * num_states 
    //  in the second RmgGemm, nwork is a matrix of num_tot_proj * (tot_states) 

    // leading dimension is num_tot_proj * 2 for noncollinear

    char *transn = "n";
    char *transc = "c";
    RmgGemm (transn, transn, dim_dnm, nstates, dim_dnm, 
            ONE_t, M_qqq,  dim_dnm, sint_kn, dim_dnm,
            ZERO_t,  sint_kn_irr, dim_dnm);

    RmgGemm (transc, transn, nstates, nstates, dim_dnm,
            ONE_t, sint_ik,  dim_dnm, sint_kn_irr, dim_dnm,
            ONE_t,  Mmn_onekpair, nstates);


    GpuFreeManaged(M_qqq);
    GpuFreeManaged(sint_ik);
    GpuFreeManaged(sint_kn);
    GpuFreeManaged(sint_ik_irr);
    GpuFreeManaged(sint_kn_irr);

}
