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
#include <boost/math/special_functions/math_fwd.hpp>
#include <boost/math/special_functions/bessel.hpp>




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
#include "Scalapack.h"

double radial_func(int radial_type, double zona, double r);
void InitDelocalizedWeight_onek(int kpt, double kvec[3], Pw &pwave);
void DelocalizedWeight_one(int kpt, double kvec[3], Pw &pwave);        

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
    G(G_in), L(L_in), wavefile(wavefile_in), nstates_tot(nstates_in), n_wannier(nwannier_in), scdm(scdm_in), 
    scdm_mu(scdm_mu_in), scdm_sigma(scdm_sigma_in), psi(psi_in)
{

    nx_grid = G.get_NX_GRID(1);
    ny_grid = G.get_NY_GRID(1);
    nz_grid = G.get_NZ_GRID(1);
    px0_grid = G.get_PX0_GRID(1);
    py0_grid = G.get_PY0_GRID(1);
    pz0_grid = G.get_PZ0_GRID(1);
    px_offset = G.get_PX_OFFSET(1);
    py_offset = G.get_PY_OFFSET(1);
    pz_offset = G.get_PZ_OFFSET(1);
    nbasis = px0_grid * py0_grid * pz0_grid;
    nbasis_noncoll = nbasis * ct.noncoll_factor;
    exclude_bands.resize(nstates_tot);
    for(int ib = 0; ib < nstates_tot; ib++) exclude_bands[ib] = false;

    if(!ct.norm_conserving_pp && ct.localize_projectors)
    {
        throw RmgFatalException() << "for ultra soft pseudopotential, set localize_projectors to be false for wannier90 interface  \n";
    }
    if(ct.is_gamma)
    {
        throw RmgFatalException() << "Wannier90 interface not programmed for gamma point  \n";
    }
    ngrid = G.get_NX_GRID(1) * G.get_NY_GRID(1) * G.get_NZ_GRID(1);
    ngrid_noncoll = ngrid * ct.noncoll_factor;

    MPI_Barrier(MPI_COMM_WORLD);

    Read_nnkpts();
    nstates = nstates_tot - num_exclude_bands;
}

template void Wannier<double> ::AmnMmn(std::string wavefile);
template void Wannier<std::complex<double>> ::AmnMmn(std::string wavefile);
template <class T> void Wannier<T>::AmnMmn(std::string wavefile)
{
    WriteWinEig();
//  setup forward beta for all of kpoints
    double kvec[3];
    if(!ct.norm_conserving_pp)
    {
        RmgTimer *RT1 = new RmgTimer("7-Wannier: init_weight");
        BaseGrid LG(Rmg_G->get_NX_GRID(1), Rmg_G->get_NY_GRID(1), Rmg_G->get_NZ_GRID(1), 1, 1, 1, 0, 1);
        int rank = Rmg_G->get_rank();
        MPI_Comm lcomm;
        MPI_Comm_split(Rmg_G->comm, rank+1, rank, &lcomm);
        LG.set_rank(0, lcomm);
        Pw pwave (LG, Rmg_L, 1, false);

        for(int kpt = pct.gridpe; kpt < ct.klist.num_k_all; kpt+=pct.grid_npes)
        {

            kvec[0] = ct.klist.k_all_cart[kpt][0];
            kvec[1] = ct.klist.k_all_cart[kpt][1];
            kvec[2] = ct.klist.k_all_cart[kpt][2];
            InitDelocalizedWeight_onek(kpt, kvec, pwave);
        }

        for(int kpt = pct.gridpe; kpt < ct.klist.num_k_ext; kpt+=pct.grid_npes)
        {
            kvec[0] = ct.klist.k_ext_cart[kpt][0];
            kvec[1] = ct.klist.k_ext_cart[kpt][1];
            kvec[2] = ct.klist.k_ext_cart[kpt][2];
            InitDelocalizedWeight_onek(kpt+ct.klist.num_k_all, kvec, pwave);

        }

        MPI_Barrier(MPI_COMM_WORLD);
        delete RT1;
        RT1 = new RmgTimer("7-Wannier: weight");
        for(int kpt = pct.gridpe; kpt < ct.klist.num_k_all; kpt+=pct.grid_npes)
        {

            kvec[0] = ct.klist.k_all_cart[kpt][0];
            kvec[1] = ct.klist.k_all_cart[kpt][1];
            kvec[2] = ct.klist.k_all_cart[kpt][2];

            DelocalizedWeight_one(kpt, kvec, pwave);        
        }
        for(int kpt = pct.gridpe; kpt < ct.klist.num_k_ext; kpt+=pct.grid_npes)
        {
            kvec[0] = ct.klist.k_ext_cart[kpt][0];
            kvec[1] = ct.klist.k_ext_cart[kpt][1];
            kvec[2] = ct.klist.k_ext_cart[kpt][2];
            DelocalizedWeight_one(kpt+ct.klist.num_k_all, kvec, pwave);        
        }
        delete RT1;
    }
    RmgTimer *RT1 = new RmgTimer("7-Wannier: Amn");
    if(scdm < 0)
    {
        SetAmn_proj();
    }
    else
    {
        SetAmn_scdm();
    }
    delete RT1;
    RT1 = new RmgTimer("7-Wannier: Mmn");
    SetMmn();
    delete RT1;
    MPI_Barrier(MPI_COMM_WORLD);
}

template <> void Wannier<double>::SetAmn_scdm()
{
    throw RmgFatalException() << "scdm need more than one gamma point  \n";
}
template <> void Wannier<std::complex<double>>::SetAmn_scdm()
{
    double tol = 1.0e-5;
    int ik_gamma = -1;

    RmgTimer *RT0;
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


    int scalapack_groups = 1;
    int last = 1;
    Scalapack  ScaL (scalapack_groups, pct.thisimg, ct.images_per_node, nstates,
            ct.scalapack_block_factor, last, pct.grid_comm);
    int desca[9];
    ScaL.ComputeDesca(nstates, ngrid_noncoll, desca);
    int n_dist, m_dist, nprow, npcol, myrow, mycol;
    m_dist = ScaL.ComputeMdim(nstates);
    n_dist = ScaL.ComputeNdim(ngrid_noncoll);
    nprow = ScaL.GetRows();
    npcol = ScaL.GetCols();
    myrow = ScaL.GetRow();
    mycol = ScaL.GetCol();

    int *piv = new int[ngrid_noncoll]();
    int *piv_dist = new int[n_dist]();
    //size_t length = (size_t)nstates * (size_t)ngrid_noncoll * sizeof(std::complex<double>);

    if(ScaL.Participates())
    {
        size_t length = (size_t)m_dist * (size_t)n_dist;
        psi_s = new std::complex<double>[length]();
        RT0 = new RmgTimer("7-Wannier: Amn: read");
        //for(int ik = 0; ik < ct.num_kpts_pe; ik++)
        std::string filename = wavefile + "_spin"+std::to_string(pct.spinpe) + "_kpt" + std::to_string(ik_gamma);
        serial_fd = open(filename.c_str(), O_RDONLY, (mode_t)0600);
        if(serial_fd < 0)
            throw RmgFatalException() << "Error! Could not open " << filename << " . Wannier Terminating.\n";
        length = nstates_tot * ngrid_noncoll * sizeof(std::complex<double>);
        std::complex<double> *psi_map;

        psi_map = (std::complex<double> *)mmap(NULL, length, PROT_READ, MAP_PRIVATE, serial_fd, 0);
        delete RT0;
        //   psi_s = (std::complex<double> *)mmap(NULL, length, PROT_READ, MAP_PRIVATE, serial_fd, 0);


        RT0 = new RmgTimer("7-Wannier: Amn: norm and scale psi");
        int st_w = -1;
        for(int st = 0; st < nstates_tot; st++)
        {

            if(exclude_bands[st]) continue;
            st_w++;

            if(myrow != ((st_w/desca[4]) % nprow ) ) continue;

            double focc = 1.0;
            double eigs = ct.kp[ik_gamma].eigs[st] * Ha_eV;
            double tem = (eigs - scdm_mu)/scdm_sigma;
            if(scdm == ISOLATED_ENTANGLEMENT) 
            {
                focc = 1.0;
            }
            else if(scdm == GAU_ENTANGLEMENT) focc = std::exp(-tem * tem);
            else if(scdm == ERFC_ENTANGLEMENT) focc = 0.5 * erfc(tem);
            else
                throw RmgFatalException() << "scdm = " << scdm << ISOLATED_ENTANGLEMENT<< "  wrong value \n";
            int st_local = ((st_w/desca[4]) /nprow) * desca[4] + st_w%desca[4];
            for(int idx = 0; idx < n_dist; idx++)
            {
                int idx_g = (idx/desca[5]) * npcol * desca[5] + mycol * desca[5] + idx%desca[5];
                psi_s[ idx * m_dist + st_local] = std::conj(psi_map[st * ngrid_noncoll + idx_g]) * focc;
            }
        }
        munmap(psi_map, length);
        close(serial_fd);

        delete RT0;

        RT0 = new RmgTimer("7-Wannier: Amn: zgeqp3");
        std::complex<double> *tau = new std::complex<double>[nstates]();
        int Lwork = -1, Lrwork = -1, info;
        std::complex<double> Lwork_tmp;
        double Lrwork_tmp;

        int ia = 1, ja = 1;
        pzgeqpf(&nstates, &ngrid_noncoll, psi_s, &ia, &ja, desca, piv_dist, tau, &Lwork_tmp, &Lwork, &Lrwork_tmp, &Lrwork, &info);

        Lwork = (int)(std::real(Lwork_tmp)) + 1;
        Lrwork = (int)(Lrwork_tmp) + 1;

        std::complex<double> *cwork = new std::complex<double>[Lwork]();
        double *rwork = new double[Lrwork]();

        pzgeqpf(&nstates, &ngrid_noncoll, psi_s, &ia, &ja, desca, piv_dist, tau, cwork, &Lwork, rwork, &Lrwork, &info);
        delete RT0;

        if(info != 0) throw RmgFatalException() << "Error! in zgeqp3 at" << __FILE__ << __LINE__ << "  Wannier Terminating.\n";
        //for(int i = 0; i < n_wannier; i++) printf("\n piv %d   %d ", piv[i], info);
        delete [] psi_s;
        delete [] rwork;
        delete [] tau;
    }

    if(myrow == 0)
    {
        for(int idx = 0; idx < n_dist; idx++)
        {
            int idx_g = (idx/desca[5]) * npcol * desca[5] + mycol * desca[5] + idx%desca[5];
            piv[idx_g] = piv_dist[idx];
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, piv, n_wannier, MPI_INT, MPI_SUM, pct.grid_comm);

    if(pct.gridpe == 0 && 1)
    {
        printf("\n wannier center in crystal unit");
        for(int iw = 0; iw < n_wannier; iw++)
        {
            int ix = (piv[iw]-1)/ny_grid/nz_grid;
            int iy = ((piv[iw]-1)/nz_grid) % ny_grid;
            int iz = (piv[iw]-1)%nz_grid;

            int spin = 1;
            if(piv[iw] >= nx_grid * ny_grid * nz_grid) spin = -1;
            printf("\n %d  %d %f  %f  %f", iw, spin, ix/(double)nx_grid, iy/(double)ny_grid, iz/(double)nz_grid);
        }
    }


    int num_q = ct.klist.num_k_all;
    std::complex<double> *psi_wan = new std::complex<double>[n_wannier * nstates];
    std::complex<double> *Amn = new std::complex<double>[num_q * n_wannier * nstates]();

    for(int iq = pct.gridpe; iq < num_q; iq+= pct.grid_npes)
    {
        int ik = ct.klist.k_map_index[iq];
        int isym = ct.klist.k_map_symm[iq];
        int isyma = std::abs(isym)-1;

        RT0 = new RmgTimer("7-Wannier: Amn: read and rotate");
        ReadRotatePsiwan(iq, ik, isym, isyma, wavefile, psi_wan, piv);
        delete RT0;


        RT0 = new RmgTimer("7-Wannier: Amn: zgesvd");
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
        delete RT0;

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

    RT0 = new RmgTimer("7-Wannier: Amn: write");
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

    delete RT0;

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
        FILE *fwin = fopen("Wannier90_rmg/wannier90.win_init", "w+");
        fprintf(fwin, "write_hr        = true");
        fprintf(fwin, "\nnum_wann        =  %d", n_wannier);
        fprintf(fwin, "\nnum_bands        =  %d", nstates);
        fprintf(fwin, "\nnum_iter        = 20");
        fprintf(fwin, "\nauto_projections= true\n");

        fprintf(fwin, "\ndis_win_max       = 19.2d0");
        fprintf(fwin, "\ndis_froz_max      =  9.8d0");

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
        fprintf(fwin, "\n%16.10f   %16.10f   %16.10f", L.get_a0(0),L.get_a0(1),L.get_a0(2));
        fprintf(fwin, "\n%16.10f   %16.10f   %16.10f", L.get_a1(0),L.get_a1(1),L.get_a1(2));
        fprintf(fwin, "\n%16.10f   %16.10f   %16.10f", L.get_a2(0),L.get_a2(1),L.get_a2(2));
        fprintf(fwin, "\nend_unit_cell_cart");
        fprintf(fwin, "\n\nmp_grid : %d %d %d", ct.kpoint_mesh[0], ct.kpoint_mesh[1], ct.kpoint_mesh[2]);

        fprintf(fwin, "\n\nbegin kpoints");
        for(int iq = 0; iq <ct.klist.num_k_all; iq++)
        {
            fprintf(fwin, "\n%16.10f   %16.10f   %16.10f", ct.klist.k_all_xtal[iq][0], ct.klist.k_all_xtal[iq][1], ct.klist.k_all_xtal[iq][2] );
        }
        fprintf(fwin, "\nend kpoints");

        fclose(fwin);

        FILE *feig= fopen("Wannier90_rmg/wannier90.eig", "w+");
        for(int iq = 0; iq <ct.klist.num_k_all; iq++)
        {
            int ik = ct.klist.k_map_index[iq];
            int st_w = 0;
            for(int st = 0; st < nstates_tot; st++)
            {
                if(exclude_bands[st]) continue;
                fprintf(feig, "%5d %5d  %18.12f\n", st_w+1, iq+1, ct.kp[ik].eigs[st] *Ha_eV);
                st_w++;
            }
        }
        fclose(feig);

    }
}


template void Wannier<double> ::Read_nnkpts();
template void Wannier<std::complex<double>> ::Read_nnkpts();
template <class T> void Wannier<T>::Read_nnkpts()
{
    double tol = 1.0e-5;
    std::ifstream fnnk("Wannier90_rmg/wannier90.nnkp");
    if(!fnnk.is_open())
    {
        WriteWinEig();
        throw RmgFatalException() << "generte wannier90.nnkpts first by running wannier90.x -pp \n";
    }

    std::string oneline;
    std::vector<std::string> kn_info;
    std::string whitespace_delims = " \n\t";

    while(getline(fnnk, oneline) )
    {
        if(oneline.find("kpoints") != std::string::npos ) break;
    }
    
    getline(fnnk, oneline);

    boost::trim(oneline);
    int num_k = std::stoi(oneline);
    if(num_k != ct.klist.num_k_all)
    {
        if(pct.gridpe == 0)
        {
            printf("\nkpoint_mesh : %d %d %d\n", ct.kpoint_mesh[0], ct.kpoint_mesh[1], ct.kpoint_mesh[2]);
            printf("\n number kpoints in wannier90.nnkpts = %d\n", num_k);
        }

        WriteWinEig();
        throw RmgFatalException() << "num kpoints in wannier90.nnkpts not match the kmesh \n";
    }

    for(int i = 0; i < ct.klist.num_k_all; i++)
    {
        getline(fnnk, oneline);
        boost::trim(oneline);
        boost::algorithm::split( kn_info, oneline, boost::is_any_of(whitespace_delims), boost::token_compress_on );
        if (kn_info.size() != 3)
        {
            std::cout << "kpoint line = " << oneline << std::endl;
            throw RmgFatalException() << "kpoints line is not 3 numbers \n";
        }
        double k_diff = std::abs(ct.klist.k_all_xtal[i][0] - std::stof(kn_info[0]));
        k_diff += std::abs(ct.klist.k_all_xtal[i][1] - std::stof(kn_info[1]));
        k_diff += std::abs(ct.klist.k_all_xtal[i][2] - std::stof(kn_info[2]));

        if (k_diff > 1.0e-5)
        {
            if(pct.gridpe == 0) 
            {
                printf("\n k from mesh %f %f %f", ct.klist.k_all_xtal[i][0], ct.klist.k_all_xtal[i][1], ct.klist.k_all_xtal[i][2]);
                printf("\n k from nnkpts %f %f %f\n", std::stof(kn_info[0]), std::stof(kn_info[1]), std::stof(kn_info[2]));
                fflush(NULL);
            }
            WriteWinEig();
            throw RmgFatalException() << "copy kpoint from wannier90.win_init to wannier.win and runwannier90.x -pp wannier90.win\n";
        }
    }

    // find the first nnkpts 
    while(getline(fnnk, oneline) )
    {
        if(oneline.find("nnkpts") != std::string::npos ) break;
    }
    getline(fnnk, oneline);

    boost::trim(oneline);
    ct.klist.num_k_nn = std::stoi(oneline);

    //k_neighbors [][][4] store the corresponding kindex from all kpoint (< ct.klist.num_k_all) or from  extra kpoints
    //k_neighbors [][][5]: the dk index
    ct.klist.k_neighbors.resize(boost::extents[ct.klist.num_k_all][ct.klist.num_k_nn][6]);

    if(pct.imgpe == 0 && ct.verbose) 
        printf("\n kpoint neighbors from nnkpts \n");

    int num_extra_k = 0;
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

            ct.klist.k_neighbors[i][j][4] = ct.klist.k_neighbors[i][j][0];
            double xkn[3], dk[3];
            int ikn = ct.klist.k_neighbors[i][j][0];
            xkn[0] = ct.klist.k_all_xtal[ikn][0] + ct.klist.k_neighbors[i][j][1];
            xkn[1] = ct.klist.k_all_xtal[ikn][1] + ct.klist.k_neighbors[i][j][2];
            xkn[2] = ct.klist.k_all_xtal[ikn][2] + ct.klist.k_neighbors[i][j][3];
            ct.klist.k_neighbors[i][j][4] = -1;

            for(int ik = 0; ik < ct.klist.num_k_all; ik++)
            {
                dk[0] = xkn[0] - ct.klist.k_all_xtal[ik][0]; 
                dk[1] = xkn[1] - ct.klist.k_all_xtal[ik][1]; 
                dk[2] = xkn[2] - ct.klist.k_all_xtal[ik][2]; 

                if( (std::abs(dk[0]) + std::abs(dk[1]) + std::abs(dk[2]) ) < tol) 
                {
                    ct.klist.k_neighbors[i][j][4] = ik;
                    break;
                }


            }

            for(int ik_ext = 0; ik_ext < num_extra_k; ik_ext++)
            {
                dk[0] = xkn[0] - ct.klist.k_ext_xtal[ik_ext][0]; 
                dk[1] = xkn[1] - ct.klist.k_ext_xtal[ik_ext][1]; 
                dk[2] = xkn[2] - ct.klist.k_ext_xtal[ik_ext][2]; 

                if( (std::abs(dk[0]) + std::abs(dk[1]) + std::abs(dk[2]) ) < tol) 
                {
                    ct.klist.k_neighbors[i][j][4] = ct.klist.num_k_all + ik_ext;
                    break;
                }


            }

            if(ct.klist.k_neighbors[i][j][4] < 0)
            {
                ct.klist.k_ext_xtal[num_extra_k][0] = xkn[0];
                ct.klist.k_ext_xtal[num_extra_k][1] = xkn[1];
                ct.klist.k_ext_xtal[num_extra_k][2] = xkn[2];
                ct.klist.k_neighbors[i][j][4] = ct.klist.num_k_all + num_extra_k;
                num_extra_k++;
            }

            if(num_extra_k > ct.klist.num_k_ext)
                throw RmgFatalException() << "too many neighbors kpoint  " <<num_extra_k<<"  " << ct.klist.num_k_ext << "\n";

            if(pct.imgpe == 0 && ct.verbose) 
                printf("%5d  %5d  %5d  %5d  %5d \n", std::stoi(kn_info[0])-1,std::stoi(kn_info[1])-1,
                        std::stoi(kn_info[2]),std::stoi(kn_info[3]),std::stoi(kn_info[4]));


        }
    }

    ct.klist.num_k_ext = num_extra_k;
    for (int kpt = 0; kpt < ct.klist.num_k_ext; kpt++) {
        double v1, v2, v3;

        v1 = ct.klist.k_ext_xtal[kpt][0] *Rmg_L.b0[0]
            + ct.klist.k_ext_xtal[kpt][1] *Rmg_L.b1[0] 
            + ct.klist.k_ext_xtal[kpt][2] *Rmg_L.b2[0];
        v2 = ct.klist.k_ext_xtal[kpt][0] *Rmg_L.b0[1]
            + ct.klist.k_ext_xtal[kpt][1] *Rmg_L.b1[1] 
            + ct.klist.k_ext_xtal[kpt][2] *Rmg_L.b2[1];
        v3 = ct.klist.k_ext_xtal[kpt][0] *Rmg_L.b0[2]
            + ct.klist.k_ext_xtal[kpt][1] *Rmg_L.b1[2] 
            + ct.klist.k_ext_xtal[kpt][2] *Rmg_L.b2[2];
        ct.klist.k_ext_cart[kpt][0] = v1 * twoPI;
        ct.klist.k_ext_cart[kpt][1] = v2 * twoPI;
        ct.klist.k_ext_cart[kpt][2] = v3 * twoPI;
    }

    while(getline(fnnk, oneline) )
    {
        if(oneline.find("exclude_bands") != std::string::npos ) break;
    }
    getline(fnnk, oneline);
    boost::trim(oneline);
    num_exclude_bands = std::stoi(oneline);
    int num_left = num_exclude_bands;
    while(num_left > 0)
    {
        std::vector<std::string> band_info;
        getline(fnnk, oneline);
        boost::trim(oneline);
        boost::algorithm::split( band_info, oneline, boost::is_any_of("-"), boost::token_compress_on );
        if(band_info.size() == 1) 
        {
            num_left -= 1; 
            int band_index = std::stoi(band_info[0]);
            exclude_bands[band_index-1] = true;
        }
        else if(band_info.size() == 2) 
        {
            int band_index1 = std::stoi(band_info[0]);
            int band_index2 = std::stoi(band_info[1]);
            num_left -= band_index2 - band_index1 +1;
            for(int ib = band_index1; ib <= band_index2; ib++) exclude_bands[ib-1] = true;
        }
        else
        {
            throw RmgFatalException() << "Error! in exclude_bands of wannier90.nnkp file. Terminating.\n";
        }
    }


    // read wannier guiding function info
    if(scdm < 0)
    {
        fnnk.clear();
        fnnk.seekg(0, std::ios::beg);
        while(getline(fnnk, oneline) )
        {
            if(oneline.find("projections") != std::string::npos ) break;
        }
        if(oneline.find("auto") != std::string::npos ) 
        {
            throw RmgFatalException() << "set scdm >= 0 to use auto projection\n";
        }
        if(ct.noncoll && oneline.find("spinor") == std::string::npos ) 
        {
            throw RmgFatalException() << "set spinors = true in .win file and get .nnkp file\n";
        }

        getline(fnnk, oneline);
        boost::trim(oneline);
        if( std::stoi(oneline) != n_wannier)
        {
            throw RmgFatalException() << "n_wannier " <<n_wannier <<" is not equal to num_proj "<<std::stoi(oneline) << "\n";
        }
        Wan_proj.resize(n_wannier);
        for(int iw = 0; iw < n_wannier; iw++)
        {
            std::vector<std::string> wan_info;
            getline(fnnk, oneline);
            boost::trim(oneline);
            boost::algorithm::split( wan_info, oneline, boost::is_any_of(whitespace_delims), boost::token_compress_on );
            Wan_proj[iw].center_xtal[0] = std::stod(wan_info[0]);
            Wan_proj[iw].center_xtal[1] = std::stod(wan_info[1]);
            Wan_proj[iw].center_xtal[2] = std::stod(wan_info[2]);
            Wan_proj[iw].l = std::stoi(wan_info[3]);
            Wan_proj[iw].m = std::stoi(wan_info[4]) -1;
            Wan_proj[iw].radial_type = std::stoi(wan_info[5]);

            L.to_cartesian(Wan_proj[iw].center_xtal, Wan_proj[iw].center_cart);
            getline(fnnk, oneline);
            boost::trim(oneline);
            boost::algorithm::split( wan_info, oneline, boost::is_any_of(whitespace_delims), boost::token_compress_on );
            Wan_proj[iw].zaxis[0] = std::stod(wan_info[0]);
            Wan_proj[iw].zaxis[1] = std::stod(wan_info[1]);
            Wan_proj[iw].zaxis[2] = std::stod(wan_info[2]);
            Wan_proj[iw].xaxis[0] = std::stod(wan_info[3]);
            Wan_proj[iw].xaxis[1] = std::stod(wan_info[4]);
            Wan_proj[iw].xaxis[2] = std::stod(wan_info[5]);
            Wan_proj[iw].zona = std::stod(wan_info[6]);

            if(pct.gridpe == 0 && ct.verbose)
            {
                printf("\n wan %d:", iw);
                printf(" %f %f %f",Wan_proj[iw].center_xtal[0],Wan_proj[iw].center_xtal[1],Wan_proj[iw].center_xtal[2]);
                printf(" %d %d %d",Wan_proj[iw].l, Wan_proj[iw].m, Wan_proj[iw].radial_type);
            }


            if(ct.noncoll)
            {
                getline(fnnk, oneline);
                boost::trim(oneline);
                boost::algorithm::split( wan_info, oneline, boost::is_any_of(whitespace_delims), boost::token_compress_on );
                Wan_proj[iw].spin = std::stoi(wan_info[0]);
                Wan_proj[iw].spin_dir[0] = std::stod(wan_info[1]);
                Wan_proj[iw].spin_dir[1] = std::stod(wan_info[2]);
                Wan_proj[iw].spin_dir[2] = std::stod(wan_info[3]);
            }
        }

    }


}
template  void Wannier<double>::ReadRotatePsiwan(int iq, int ikindex, int isym, int isyma, std::string wavefile, double *psi_wan, int *piv);
template  void Wannier<std::complex<double>>::ReadRotatePsiwan(int iq, int ikindex, int isym, int isyma, std::string wavefile, std::complex<double>
        *psi_wan, int *piv);
template <class T> void Wannier<T>::ReadRotatePsiwan(int iq, int ikindex, int isym, int isyma, std::string wavefile, T *psi_wan, int *piv)
{

    size_t length = nstates * ngrid_noncoll * sizeof(T);
    T *psi_map;

    std::string filename = wavefile + "_spin"+std::to_string(pct.spinpe) + "_kpt" + std::to_string(ikindex);
    int serial_fd = open(filename.c_str(), O_RDONLY, (mode_t)0600);
    if(serial_fd < 0)
        throw RmgFatalException() << "Error! Could not open " << filename << " . Terminating.\n";

    psi_map = (T *)mmap(NULL, length, PROT_READ, MAP_PRIVATE, serial_fd, 0);

    double *norm_coeff = new double[nstates_tot]();
    for(int st = 0; st < nstates; st++)
    {
        double vel = L.get_omega() / ((double)ngrid);
        for(int idx = 0; idx < ngrid_noncoll; idx++)
        {
            norm_coeff[st] += std::norm(psi_map[st * ngrid_noncoll + idx]);
        }
        norm_coeff[st] = std::sqrt(norm_coeff[st] * vel);
    }

    int ix, iy, iz, ixx, iyy, izz;
    std::complex<double> phase_center;
    std::complex<double> *psi_wan_C = (std::complex<double> *)psi_wan;
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


        symm_ijk(&Rmg_Symm->sym_rotate[isyma *9], &Rmg_Symm->ftau_wave[isyma*3], 
                ix, iy, iz, ixx, iyy, izz, nx_grid, ny_grid, nz_grid);
        int st_w = -1;
        for(int st = 0; st < nstates_tot; st++)
        {
            if(exclude_bands[st]) continue;
            st_w++;
            double focc = 1.0;
            double eigs = ct.kp[ikindex].eigs[st] * Ha_eV;
            double tem = (eigs - scdm_mu)/scdm_sigma;
            if(scdm == ISOLATED_ENTANGLEMENT ) focc = 1.0;
            else if(scdm == GAU_ENTANGLEMENT) focc = std::exp(-tem * tem);
            else if(scdm == ERFC_ENTANGLEMENT) focc = 0.5 * erfc(tem);
            else
            {
                throw RmgFatalException() << "scdm = " << scdm << "  wrong value \n";
            }

            if(ct.noncoll) 
            {
                std::complex<double> up = psi_map[st_w * ngrid_noncoll + ixx * ny_grid * nz_grid + iyy * nz_grid + izz];
                std::complex<double> dn = psi_map[st_w * ngrid_noncoll + ngrid + ixx * ny_grid * nz_grid + iyy * nz_grid + izz];
                std::complex<double> up_rot, dn_rot;
                if(isym <0)
                {
                    up_rot = std::conj(Rmg_Symm->rot_spin_wave[isyma][0][0]) * up + std::conj(Rmg_Symm->rot_spin_wave[isyma][0][1]) * dn;
                    dn_rot = std::conj(Rmg_Symm->rot_spin_wave[isyma][1][0]) * up + std::conj(Rmg_Symm->rot_spin_wave[isyma][1][1]) * dn;
                    if(ispin == 0) psi_wan_C[iw * nstates + st_w] = - std::conj(dn_rot);
                    if(ispin == 1) psi_wan_C[iw * nstates + st_w] = std::conj(up_rot);
                }
                else
                {
                    up_rot = Rmg_Symm->rot_spin_wave[isyma][0][0] * up + Rmg_Symm->rot_spin_wave[isyma][0][1] * dn;
                    dn_rot = Rmg_Symm->rot_spin_wave[isyma][1][0] * up + Rmg_Symm->rot_spin_wave[isyma][1][1] * dn;
                    if(ispin == 0) psi_wan_C[iw * nstates + st_w] = up_rot;
                    if(ispin == 1) psi_wan_C[iw * nstates + st_w] = dn_rot;
                }

            }
            else if(isym > 0)
            {
                psi_wan[iw * nstates + st_w] = 
                    psi_map[st * ngrid + ixx * ny_grid * nz_grid + iyy * nz_grid + izz];
            }
            else
            {
                psi_wan[iw * nstates + st_w] = 
                    MyConj(psi_map[st * ngrid + ixx * ny_grid * nz_grid + iyy * nz_grid + izz]);
            }

            psi_wan_C[iw * nstates + st_w] = MyConj(psi_wan[iw * nstates + st_w]) * focc * phase_center /norm_coeff[st];

        }

    }

    delete [] norm_coeff;
    munmap(psi_map, length);
    close(serial_fd);
}

template  void Wannier<double>::ReadRotatePsi(int ikindex, int isym, int isyma, std::string wavefile, double *psi_k);
template  void Wannier<std::complex<double>>::ReadRotatePsi(int ikindex, int isym, int isyma, std::string wavefile, std::complex<double> *psi_k);
template <class T> void Wannier<T>::ReadRotatePsi(int ikindex, int isym, int isyma, std::string wavefile, T *psi_k)
{

    size_t length = nstates * ngrid_noncoll * sizeof(T);
    T *psi_map;

    std::string filename = wavefile + "_spin"+std::to_string(pct.spinpe) + "_kpt" + std::to_string(ikindex);
    int serial_fd = open(filename.c_str(), O_RDONLY, (mode_t)0600);
    if(serial_fd < 0)
        throw RmgFatalException() << "Error! Could not open " << filename << " . Terminating.\n";

    psi_map = (T *)mmap(NULL, length, PROT_READ, MAP_PRIVATE, serial_fd, 0);


    int ixx, iyy, izz;
    std::complex<double> *psi_k_C = (std::complex<double> *)psi_k;
    for (int ix = 0; ix < px0_grid; ix++) {
        for (int iy = 0; iy < py0_grid; iy++) {
            for (int iz = 0; iz < pz0_grid; iz++) {

                int ix_g = ix + px_offset;
                int iy_g = iy + py_offset;
                int iz_g = iz + pz_offset;

                symm_ijk(&Rmg_Symm->sym_rotate[isyma *9], &Rmg_Symm->ftau_wave[isyma*3], 
                        ix_g, iy_g, iz_g, ixx, iyy, izz, nx_grid, ny_grid, nz_grid);

                for(int st = 0; st < nstates; st++)
                {
                    if(ct.noncoll) 
                    {
                        std::complex<double> up = psi_map[st * ngrid_noncoll + ixx * ny_grid * nz_grid + iyy * nz_grid + izz];
                        std::complex<double> dn = psi_map[st * ngrid_noncoll + ngrid + ixx * ny_grid * nz_grid + iyy * nz_grid + izz];
                        std::complex<double> up_rot, dn_rot;
                        if(isym < 0 || Rmg_Symm->time_rev[isyma])
                        {
                            up_rot = std::conj(Rmg_Symm->rot_spin_wave[isyma][0][0]) * up + std::conj(Rmg_Symm->rot_spin_wave[isyma][0][1]) * dn;
                            dn_rot = std::conj(Rmg_Symm->rot_spin_wave[isyma][1][0]) * up + std::conj(Rmg_Symm->rot_spin_wave[isyma][1][1]) * dn;
                            psi_k_C[st * nbasis_noncoll + ix * py0_grid * pz0_grid + iy * pz0_grid + iz] = - std::conj(dn_rot);
                            psi_k_C[st * nbasis_noncoll + nbasis + ix * py0_grid * pz0_grid + iy * pz0_grid + iz] = std::conj(up_rot);
                        }
                        else
                        {
                            up_rot = (Rmg_Symm->rot_spin_wave[isyma][0][0]) * up + (Rmg_Symm->rot_spin_wave[isyma][0][1]) * dn;
                            dn_rot = (Rmg_Symm->rot_spin_wave[isyma][1][0]) * up + (Rmg_Symm->rot_spin_wave[isyma][1][1]) * dn;
                            psi_k_C[st * nbasis_noncoll + ix * py0_grid * pz0_grid + iy * pz0_grid + iz] = up_rot;
                            psi_k_C[st * nbasis_noncoll + nbasis + ix * py0_grid * pz0_grid + iy * pz0_grid + iz] = dn_rot;
                        }

                    }
                    else if(isym >= 0)
                    {
                        psi_k[st * nbasis + ix * py0_grid * pz0_grid + iy * pz0_grid + iz]
                            = (psi_map[st * ngrid + ixx * ny_grid * nz_grid + iyy * nz_grid + izz]);
                    }
                    else
                    {
                        psi_k[st * nbasis + ix * py0_grid * pz0_grid + iy * pz0_grid + iz]
                            = MyConj(psi_map[st * ngrid + ixx * ny_grid * nz_grid + iyy * nz_grid + izz]);
                    }
                }
            }
        }
    }

    munmap(psi_map, length);
    close(serial_fd);

}


template  void Wannier<double>::SetMmn();
template  void Wannier<std::complex<double>>::SetMmn();
template <class T> void Wannier<T>::SetMmn()
{

    double tol = 1.0e-5;
    int num_q = ct.klist.num_k_all;
    int num_kn = ct.klist.num_k_nn;

    size_t length = nstates * nbasis_noncoll * sizeof(T);
    T *psi_k = (T *)RmgMallocHost(length);
    T *psi_q = (T *)RmgMallocHost(length);
    T *Mmn = (T *)RmgMallocHost(num_q * num_kn * nstates * nstates * sizeof(T));

    for(int idx = 0; idx < num_q * num_kn * nstates * nstates; idx++) Mmn[idx] = 0.0;
    double vel = L.get_omega() / ((double)ngrid);
    T alpha(vel);
    T beta = 0.0;

    T kphase(1.0);
    double kr;

    RmgTimer *RT1 = new RmgTimer("7-Wannier: Mmn: set dk");
    double dk[ct.klist.num_k_nn][3], dk1[3];
    for(int ikn = 0; ikn < ct.klist.num_k_nn; ikn++)
    {
        int ikn_map = ct.klist.k_neighbors[0][ikn][4];
        if(ikn_map < ct.klist.num_k_all)
        {
            dk[ikn][0] = -ct.klist.k_all_xtal[0][0] + ct.klist.k_all_xtal[ikn_map][0]; 
            dk[ikn][1] = -ct.klist.k_all_xtal[0][1] + ct.klist.k_all_xtal[ikn_map][1]; 
            dk[ikn][2] = -ct.klist.k_all_xtal[0][2] + ct.klist.k_all_xtal[ikn_map][2]; 

        }
        else
        {
            dk[ikn][0] = -ct.klist.k_all_xtal[0][0] + ct.klist.k_ext_xtal[ikn_map - ct.klist.num_k_all][0];
            dk[ikn][1] = -ct.klist.k_all_xtal[0][1] + ct.klist.k_ext_xtal[ikn_map - ct.klist.num_k_all][1];
            dk[ikn][2] = -ct.klist.k_all_xtal[0][2] + ct.klist.k_ext_xtal[ikn_map - ct.klist.num_k_all][2];
        }
    }

    for(int ik = 0; ik < ct.klist.num_k_all; ik++)
    {
        for(int ikn = 0; ikn < ct.klist.num_k_nn; ikn++)
        {
            int ikn_map = ct.klist.k_neighbors[ik][ikn][4];
            if(ikn_map < ct.klist.num_k_all)
            {
                dk1[0] = -ct.klist.k_all_xtal[ik][0] + ct.klist.k_all_xtal[ikn_map][0]; 
                dk1[1] = -ct.klist.k_all_xtal[ik][1] + ct.klist.k_all_xtal[ikn_map][1]; 
                dk1[2] = -ct.klist.k_all_xtal[ik][2] + ct.klist.k_all_xtal[ikn_map][2]; 

            }
            else
            {
                dk1[0] = -ct.klist.k_all_xtal[ik][0] + ct.klist.k_ext_xtal[ikn_map - ct.klist.num_k_all][0];
                dk1[1] = -ct.klist.k_all_xtal[ik][1] + ct.klist.k_ext_xtal[ikn_map - ct.klist.num_k_all][1];
                dk1[2] = -ct.klist.k_all_xtal[ik][2] + ct.klist.k_ext_xtal[ikn_map - ct.klist.num_k_all][2];
            }

            ct.klist.k_neighbors[ik][ikn][5] = -1;
            for(int idk = 0; idk < ct.klist.num_k_nn; idk++)
            {
                if(std::abs(dk1[0] - dk[idk][0]) + std::abs(dk1[1] - dk[idk][1]) + std::abs(dk1[2] - dk[idk][2]) < tol )
                {
                    ct.klist.k_neighbors[ik][ikn][5] = idk;
                    break;
                }

            }

            if(ct.klist.k_neighbors[ik][ikn][5] < 0)
            {
                printf("\n dk are different %d %d %f %f %f  ", ik, ikn, dk1[0], dk1[1], dk1[2]);
                throw RmgFatalException() << "Kpoint neighbors has different distance \n" ;
            }
        }
    }

    delete RT1;
    std::complex<double> *qqq_dk, *qqq_dk_so, *qq_dk_one, *qq_dk_so_one;
    qqq_dk = new std::complex<double>[ct.klist.num_k_nn * Atoms.size() * ct.max_nl * ct.max_nl];
    qqq_dk_so = new std::complex<double>[4*ct.klist.num_k_nn * Atoms.size() * ct.max_nl * ct.max_nl];

    RT1 = new RmgTimer("7-Wannier: Mmn: qqq_dk");
    for(int ikn = 0; ikn < ct.klist.num_k_nn; ikn++)
    {
        qq_dk_one = &qqq_dk[ikn * Atoms.size() * ct.max_nl * ct.max_nl];
        qq_dk_so_one = &qqq_dk_so[ikn * 4 * Atoms.size() * ct.max_nl * ct.max_nl];
        get_qqq_dk(dk[ikn], qq_dk_one, qq_dk_so_one);
    }
    delete RT1;
    //double *kphase_R = (double *)&kphase;

    int num_tot_proj = 0;
    for(size_t ion = 0; ion < Atoms.size(); ion++)
    {
        SPECIES &AtomType = Species[Atoms[ion].species];
        num_tot_proj += AtomType.nh;
    }

    T *Nlweight_k = NULL;
    T *Nlweight_q = NULL;
    if(!ct.norm_conserving_pp)
    {
        Nlweight_k = (T *)RmgMallocHost(num_tot_proj * nbasis* sizeof(T));
        Nlweight_q = (T *)RmgMallocHost(num_tot_proj * nbasis* sizeof(T));
    }
    std::complex<double> *kphase_C = (std::complex<double> *)&kphase;
    for(int ik = 0; ik < num_q; ik++)
    {
        int ik_irr = ct.klist.k_map_index[ik];
        int isym = ct.klist.k_map_symm[ik];
        int isyma = std::abs(isym) -1;

        RT1 = new RmgTimer("7-Wannier: Mmn: read and rotate");
        if(isym == 1 ) {
            ReadPsiFromSingleFile(ik_irr, wavefile, psi_k);
        }
        else {
            ReadRotatePsi(ik_irr, isym, isyma, wavefile, psi_k);
        }
        delete RT1;

        if(!ct.norm_conserving_pp)
        {
            RT1 = new RmgTimer("7-Wannier: Mmn: us: read NL");
            std::string filename;
            filename = "PROJECTORS/NLprojectors_kpt" + std::to_string(ik);
            std::complex<double> *Nlweight_C = (std::complex<double> *)Nlweight_k;
            ReadNlweight(filename, num_tot_proj, Nlweight_C);

            delete RT1;
        }

        for(int ikn = 0; ikn < num_kn; ikn++)
        {

            if(!ct.norm_conserving_pp)
            {
                RT1 = new RmgTimer("7-Wannier: Mmn: us: read NL");
                int ikn_map = ct.klist.k_neighbors[ik][ikn][4];
                std::string filename = "PROJECTORS/NLprojectors_kpt" + std::to_string(ikn_map);
                std::complex<double> *Nlweight_C = (std::complex<double> *)Nlweight_q;
                ReadNlweight(filename, num_tot_proj, Nlweight_C);
                delete RT1;
            }

            int idk = ct.klist.k_neighbors[ik][ikn][5];
            qq_dk_one = &qqq_dk[idk * Atoms.size() * ct.max_nl * ct.max_nl];
            qq_dk_so_one = &qqq_dk_so[idk * 4 * Atoms.size() * ct.max_nl * ct.max_nl];

            int ikn_index = ct.klist.k_neighbors[ik][ikn][0];
            int ikn_irr = ct.klist.k_map_index[ikn_index];
            int isym_kn = ct.klist.k_map_symm[ikn_index];
            int isyma_kn = std::abs(isym_kn) -1;

            RT1 = new RmgTimer("7-Wannier: Mmn: read and rotate");
            if(isym_kn == 1 ) {
                ReadPsiFromSingleFile(ikn_irr, wavefile, psi_q);
            }
            else {
                ReadRotatePsi(ikn_irr, isym_kn, isyma_kn, wavefile, psi_q);
            }

            delete RT1;

            RT1 = new RmgTimer("7-Wannier: Mmn: phase");
            //  phase factor when kneight need a periodic BZ folding
            for (int iz = 0; iz < pz0_grid; iz++) {
                for (int iy = 0; iy < py0_grid; iy++) {
                    for (int ix = 0; ix < px0_grid; ix++) {
                        kr = ct.klist.k_neighbors[ik][ikn][1] * (ix+px_offset)/(double)nx_grid
                            +ct.klist.k_neighbors[ik][ikn][2] * (iy+py_offset)/(double)ny_grid
                            +ct.klist.k_neighbors[ik][ikn][3] * (iz+pz_offset)/(double)nz_grid;
                        *kphase_C = std::exp( std::complex<double>(0.0, -kr * twoPI));
                        for(int st = 0; st < nstates; st++)
                        {
                            psi_q[st * nbasis_noncoll + ix * py0_grid * pz0_grid + iy * pz0_grid + iz] *= kphase;
                            if(ct.noncoll)
                                psi_q[st * nbasis_noncoll + nbasis + ix * py0_grid * pz0_grid + iy * pz0_grid + iz] *= kphase;
                        }

                    }
                }
            }
            delete RT1;

            RT1 = new RmgTimer("7-Wannier: Mmn: gemm");
            RmgGemm("C", "N", nstates, nstates, nbasis_noncoll, alpha, psi_k, nbasis_noncoll, psi_q, nbasis_noncoll,
                    beta, &Mmn[(ik*num_kn+ikn)*nstates*nstates], nstates);
            delete RT1;


            RT1 = new RmgTimer("7-Wannier: Mmn: us");
            if(!ct.norm_conserving_pp)
            {
                Mmn_us(ik, ikn, psi_k, nstates, psi_q, nstates, &Mmn[(ik*num_kn+ikn)*nstates*nstates], qq_dk_one, qq_dk_so_one,
                        num_tot_proj, Nlweight_k, Nlweight_q);
            }
            delete RT1;

        }
    }

    RT1 = new RmgTimer("7-Wannier: Mmn: Reduce");
    int count = num_q * num_kn * nstates * nstates;
    MPI_Allreduce(MPI_IN_PLACE, Mmn, count, MPI_DOUBLE_COMPLEX, MPI_SUM, pct.grid_comm);
    delete RT1;


    RT1 = new RmgTimer("7-Wannier: Mmn: write");
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

    delete RT1;

    RmgFreeHost(psi_k);
    RmgFreeHost(psi_q);
    RmgFreeHost(Mmn);
    if(Nlweight_k != NULL) RmgFreeHost(Nlweight_k);
    if(Nlweight_q != NULL) RmgFreeHost(Nlweight_q);

}

template void Wannier<double>::ReadNlweight(std::string filename, int nh, std::complex<double> *Nlweight_oneatom);
template void Wannier<std::complex<double>>::ReadNlweight(std::string filename, int nh, std::complex<double> *Nlweight_oneatom);
template <class T> void Wannier<T>::ReadNlweight(std::string filename, int nh, std::complex<double> *Nlweight_oneatom)
{

    MPI_Datatype wftype = MPI_DOUBLE_COMPLEX;

    int sizes_c[4];
    int subsizes_c[4];
    int starts_c[4];

    sizes_c[0] = nh;
    sizes_c[1] = G.get_NX_GRID(1);
    sizes_c[2] = G.get_NY_GRID(1);
    sizes_c[3] = G.get_NZ_GRID(1);

    subsizes_c[0] = nh;
    subsizes_c[1] = G.get_PX0_GRID(1);
    subsizes_c[2] = G.get_PY0_GRID(1);
    subsizes_c[3] = G.get_PZ0_GRID(1);

    starts_c[0] = 0;
    starts_c[1] = G.get_PX_OFFSET(1);
    starts_c[2] = G.get_PY_OFFSET(1);
    starts_c[3] = G.get_PZ_OFFSET(1);

    int order = MPI_ORDER_C;
    MPI_Info fileinfo;
    MPI_Datatype grid_c;
    MPI_Status status;

    MPI_Type_create_subarray(4, sizes_c, subsizes_c, starts_c, order, wftype, &grid_c);
    MPI_Type_commit(&grid_c);

    MPI_Info_create(&fileinfo);

    int amode = MPI_MODE_RDWR;
    MPI_File mpi_fhand ;


    MPI_File_open(G.comm, filename.c_str(), amode, fileinfo, &mpi_fhand);
    MPI_Offset disp = 0;

    MPI_File_set_view(mpi_fhand, disp, wftype, grid_c, "native", MPI_INFO_NULL);
    int dis_dim = G.get_P0_BASIS(1) * nh;

    MPI_File_read_all(mpi_fhand, Nlweight_oneatom, dis_dim, wftype, &status);
    MPI_Barrier(G.comm);
    MPI_File_close(&mpi_fhand);
    MPI_Type_free(&grid_c);

}
template  void Wannier<double>::Mmn_us(int, int, double *, int, double *, int, double *, std::complex<double> *, 
        std::complex<double> *, int, double*, double *);
template  void Wannier<std::complex<double>>::Mmn_us(int, int, std::complex<double> *, int, std::complex<double> *, int,
        std::complex<double> *, std::complex<double> *, std::complex<double> *, int, std::complex<double> *, std::complex<double> *);
template <class T> void Wannier<T>::Mmn_us(int ik, int ikn, T *psi_k, int num_st_k, T *psi_q, int num_st_q, T *Mmn_onekpair, 
        std::complex<double> *qq_dk_one, std::complex<double> *qq_dk_so_one, int num_tot_proj, T *Nlweight_k, T *Nlweight_q)
{



    RmgTimer *RT1;
    T ZERO_t(0.0);
    T ONE_t(1.0);
    double vel = L.get_omega() / ((double)ngrid);
    T alpha(vel);
    int tot_st_k = num_st_k * ct.noncoll_factor;
    int tot_st_q = num_st_q * ct.noncoll_factor;

    size_t alloc = (size_t)num_tot_proj * (size_t)num_st_k * ct.noncoll_factor * sizeof(T);
    T *sint_ik = (T *)RmgMallocHost(alloc);
    T *sint_tem = (T *)RmgMallocHost(alloc);

    alloc = (size_t)num_tot_proj * (size_t)num_st_q * ct.noncoll_factor * sizeof(T);
    T *sint_kn = (T *)RmgMallocHost(alloc);


    RT1 = new RmgTimer("7-Wannier: Mmn: us: betapsi");
    RmgGemm ("C", "N", num_tot_proj, tot_st_k, nbasis, alpha,
            Nlweight_k, nbasis, psi_k, nbasis, ZERO_t, sint_ik, num_tot_proj);
    size_t count = (size_t)num_tot_proj * (size_t)num_st_k * ct.noncoll_factor;
    MPI_Allreduce(MPI_IN_PLACE, sint_ik, count, MPI_DOUBLE_COMPLEX, MPI_SUM, G.comm);

    delete RT1;


    RT1 = new RmgTimer("7-Wannier: Mmn: us: betapsi");
    RmgGemm ("C", "N", num_tot_proj, tot_st_q, nbasis, alpha,
            Nlweight_q, nbasis, psi_q, nbasis, ZERO_t, sint_kn, num_tot_proj);
    count = (size_t)num_tot_proj * (size_t)num_st_q * ct.noncoll_factor;
    MPI_Allreduce(MPI_IN_PLACE, sint_kn, count, MPI_DOUBLE_COMPLEX, MPI_SUM, G.comm);
    delete RT1;

    RT1 = new RmgTimer("7-Wannier: Mmn: us: set qqq");
    double *qqq;

    int M_cols = (size_t)num_tot_proj * ct.noncoll_factor;
    size_t alloc1 = M_cols * M_cols;
    T *M_qqq = (T *)RmgMallocHost(sizeof(T) * alloc1);
    std::complex<double> *M_qqq_C = (std::complex<double> *) M_qqq;

    for (size_t i = 0; i < alloc1; i++)
    {
        M_qqq[i] = ZERO_t;
    }


    // set up M_qqq and M_dnm, this can be done outside in the
    // init.c or get_ddd get_qqq, we need to check the order
    int proj_index = 0;

    for (size_t ion = 0; ion < Atoms.size(); ion++)
    {

        /*Actual index of the ion under consideration*/

        SPECIES &AtomType = Species[Atoms[ion].species];
        int nh = AtomType.nh;

        qqq = Atoms[ion].qqq;

        for (int i = 0; i < nh; i++)
        {
            int inh = i * nh;
            for (int j = 0; j < nh; j++)
            {

                if(ct.is_gamma)
                {
                    int idx = (proj_index + i) * num_tot_proj + proj_index + j;
                    M_qqq[idx] = (T)qqq[inh+j];
                }
                else if(!ct.noncoll)
                {
                    int idx = (proj_index + i) * num_tot_proj + proj_index + j;
                    M_qqq_C[idx] = qq_dk_one[ion * ct.max_nl * ct.max_nl + inh+j];

                }
                else
                {
                    int it0 = proj_index + i;
                    int jt0 = proj_index + j;
                    int it1 = proj_index + i + num_tot_proj;
                    int jt1 = proj_index + j + num_tot_proj;
                    M_qqq_C[it0 * num_tot_proj * 2 + jt0] = qq_dk_so_one[ion * ct.max_nl * ct.max_nl*4 + inh+j + 0 * nh *nh];
                    M_qqq_C[it0 * num_tot_proj * 2 + jt1] = qq_dk_so_one[ion * ct.max_nl * ct.max_nl*4 + inh+j + 1 * nh *nh];
                    M_qqq_C[it1 * num_tot_proj * 2 + jt0] = qq_dk_so_one[ion * ct.max_nl * ct.max_nl*4 + inh+j + 2 * nh *nh];
                    M_qqq_C[it1 * num_tot_proj * 2 + jt1] = qq_dk_so_one[ion * ct.max_nl * ct.max_nl*4 + inh+j + 3 * nh *nh];
                }
            }

        }

        proj_index += nh;
    }

    delete RT1;

    RT1 = new RmgTimer("7-Wannier: Mmn: us: gemm");
    int dim_dnm = num_tot_proj * ct.noncoll_factor;

    //M_dnm: dim_dnm * dim_dnm matrxi
    //sint_compack: dim_dnm * num_states == num_tot_proj * ct.noncoll_factor * num_states
    //nwork: dim_dnm * num_states == num_tot_proj * ct.noncoll_factor * num_states
    //  in the first RmgGemm, nwork is a matrix of (dim_dnm) * num_states 
    //  in the second RmgGemm, nwork is a matrix of num_tot_proj * (tot_states) 

    // leading dimension is num_tot_proj * 2 for noncollinear

    char *transn = "n";
    char *transc = "c";
    RmgGemm (transn, transn, dim_dnm, num_st_q, dim_dnm, 
            ONE_t, M_qqq,  dim_dnm, sint_kn, dim_dnm,
            ZERO_t,  sint_tem, dim_dnm);

    alpha = 1.0/pct.grid_npes;
    RmgGemm (transc, transn, num_st_k, num_st_q, dim_dnm,
            alpha, sint_ik,  dim_dnm, sint_tem, dim_dnm,
            ONE_t,  Mmn_onekpair, num_st_k);

    delete RT1;

    RmgFreeHost(M_qqq);
    RmgFreeHost(sint_ik);
    RmgFreeHost(sint_kn);
    RmgFreeHost(sint_tem);

}

double radial_func(int radial_type, double zona, double r)
{
    double val = 0.0;
    if(radial_type == 1)  val =  2.0 * std::pow(zona, 3.0/2.0) * std::exp(-zona * r);
    if(radial_type == 2)  val =  2.0/std::sqrt(8.0) * std::pow(zona, 3.0/2.0) * 
        (2.0 - zona * r) * std::exp(-0.5 * zona * r);
    if(radial_type == 3)  val =  std::sqrt(4.0/27.0) * std::pow(zona, 3.0/2.0) * 
        (1.0 - 2.0/3.0 * zona * r + 2.0 * (zona * r) * (zona*r)/27.0) * std::exp(-zona * r/3.0);

    return val;
}

void InitRadialfunc_Gspace(int lval, int radial_type, double zona, std::complex<double> *radialfunc_g, int gnum, double delta_g)
{
    using boost::math::policies::policy;
    using boost::math::policies::promote_double;
    typedef policy<promote_double<false> > bessel_policy;
    double xmin=-6.0, rmax = 10.0, dx = 0.025;
    int rg_points = std::round( (std::log(rmax) - xmin)/dx + 1.0);
    double *func_r = new double[rg_points];
    double *r = new double[rg_points];
    double *rab = new double[rg_points];
    double *work1 = new double[rg_points];

    zona = zona / 0.529177;
    for(int i = 0; i < rg_points; i++)
    {
        r[i] = std::exp(xmin + i * dx)/zona;
        rab[i] = r[i] * dx;
        func_r[i] = radial_func(radial_type, zona, r[i]) * r[i] * r[i];

    }

    double alpha = (double)lval + 0.5;

    for(int ig = 0;ig < gnum;ig++)
    {
        double gval = ig * delta_g;

        for (int idx = 0; idx < rg_points; idx++)
        {
            double jarg = r[idx] * gval;
            if(jarg < 1.0e-10) jarg = 1.0e-10;  // take care of sigularity
            work1[idx] = func_r[idx] * sqrt(PI/(2.0*jarg)) * boost::math::cyl_bessel_j(alpha, jarg, bessel_policy());

        }                   /* end for */
        radialfunc_g[ig] =  radint1 (work1, r, rab, rg_points);
    }
    delete [] func_r;
    delete [] r;
    delete [] rab;
    delete [] work1;

}

template  void Wannier<double>::InitGuideFunc();
template  void Wannier<std::complex<double>>::InitGuideFunc();
template <class T> void Wannier<T>::InitGuideFunc()
{
    double tol = 1.0e-5;
    std::pair<int, double> onezona_type;
    zona_index.resize(Wan_proj.size());
    int tot_LM;
    lmax = 0;
    RmgTimer *RT1 = new RmgTimer("7-Wannier: Amn: gf: zona");
    for(size_t iw = 1; iw< Wan_proj.size(); iw++)
    {
        lmax = std::max(lmax, Wan_proj[iw].l);
        // l = -4, -5, spd hybridization
        //  l = -1, -2, -3, sp1, sp2, sp3 hybridization
        if(Wan_proj[iw].l <=-4) lmax = std::max(lmax, 2);
        else if(Wan_proj[iw].l < 0) lmax = std::max(lmax, 1);

        bool zona_found = false;
        for(int izona = 0; izona < (int)zona_list.size(); izona++)
        {
            if( std::abs(Wan_proj[iw].zona - zona_list[izona].second) < tol && Wan_proj[iw].radial_type == zona_list[izona].first )
            {
                zona_index[iw] = izona;
                zona_found = true;
                break;
            }

        }
        if(!zona_found)
        {
            onezona_type = std::make_pair(Wan_proj[iw].radial_type, Wan_proj[iw].zona);
            zona_list.push_back(onezona_type);
            zona_index[iw] = zona_list.size() -1;
        }

    }

    delete RT1;
    RT1 = new RmgTimer("7-Wannier: Amn: gf: csph");
    // setup rotation of axis and hybridization of differetn Ylms. such sp3 ...
    tot_LM = (lmax + 1) * (lmax+1);
    double_2d_array r_rand;
    csph.resize(boost::extents[Wan_proj.size()][tot_LM]);
    r_rand.resize(boost::extents[tot_LM * tot_LM][3]);
    std::srand(224);
    for(int i = 0; i < tot_LM; i++)
    {
        r_rand[i][0] = (2.0 *std::rand())/RAND_MAX -1.0;
        r_rand[i][1] = (2.0 *std::rand())/RAND_MAX -1.0;
        r_rand[i][2] = (2.0 *std::rand())/RAND_MAX -1.0;
    }
    double *ylm_array = new double[tot_LM * tot_LM];
    double *ylm_invert = new double[tot_LM * tot_LM];
    double *ylm_wan = new double[tot_LM];

    for (int i = 0; i < tot_LM; i++)
    {
        double r[3];
        r[0] = r_rand[i][0];
        r[1] = r_rand[i][1];
        r[2] = r_rand[i][2];
        for(int L = 0; L <= lmax; L++)
            for(int M = 0; M < 2*L+1; M++)
            {
                int lm = L * L + M;
                ylm_array[i * tot_LM + lm] = Ylm(L, M, r);
            }
    }

    for(int i = 0; i < tot_LM * tot_LM; i++) ylm_invert[i] = 0.0;
    for(int i = 0; i < tot_LM; i++) ylm_invert[i * tot_LM + i] = 1.0;

    int info, ipvt[tot_LM];
    dgesv (&tot_LM, &tot_LM, ylm_array, &tot_LM, ipvt, ylm_invert, &tot_LM, &info);

    double um[3][3], yaxis[3];
    double bs2 = 1./sqrt(2.0);
    double bs3=1.0/sqrt(3.0);
    double bs6 = 1.0/sqrt(6.0);
    double bs12 = 1.0/sqrt(12.0);

    // mapping m value in wannier to rmg 
    int map_m[4][7];
    map_m[0][0] = 0;
    map_m[1][0] = 2;
    map_m[1][1] = 0;
    map_m[1][2] = 1;
    map_m[2][0] = 2;
    map_m[2][1] = 1;
    map_m[2][2] = 3;
    map_m[2][3] = 4;
    map_m[2][4] = 0;
    map_m[3][0] = 0;
    map_m[3][1] = 1;
    map_m[3][2] = 2;
    map_m[3][3] = 3;
    map_m[3][4] = 4;
    map_m[3][5] = 5;
    map_m[3][6] = 6;

    for(int iw = 0; iw < (int)Wan_proj.size(); iw++)
    {
        yaxis[0] = Wan_proj[iw].zaxis[1] * Wan_proj[iw].xaxis[2] - Wan_proj[iw].zaxis[2] * Wan_proj[iw].xaxis[1];
        yaxis[1] = Wan_proj[iw].zaxis[2] * Wan_proj[iw].xaxis[0] - Wan_proj[iw].zaxis[0] * Wan_proj[iw].xaxis[2];
        yaxis[2] = Wan_proj[iw].zaxis[0] * Wan_proj[iw].xaxis[1] - Wan_proj[iw].zaxis[1] * Wan_proj[iw].xaxis[0];
        for(int i = 0; i < 3; i++)
        {
            um[0][i] = Wan_proj[iw].xaxis[i];
            um[1][i] = yaxis[i];
            um[2][i] = Wan_proj[iw].zaxis[i];
        }

        double r_rot[3];
        for(int ir = 0; ir < tot_LM; ir++)
        {
            for(int i = 0; i < 3; i++)
            {
                r_rot[i] = 0.0;
                for(int j = 0; j < 3; j++)
                {
                    r_rot[i] += um[i][j] * r_rand[ir][j]; 
                }

            }

            if(Wan_proj[iw].l >= 0) 
            {
                int m_rmg = map_m[Wan_proj[iw].l][Wan_proj[iw].m];
                ylm_wan[ir] = Ylm(Wan_proj[iw].l, m_rmg, r_rot);
            }
            if(Wan_proj[iw].l == -1)  // sp1 hybridization
            {
                if(Wan_proj[iw].m == 0) ylm_wan[ir] = bs2 * (Ylm(0, 0, r_rot) + Ylm(1, 0, r_rot));
                if(Wan_proj[iw].m == 1) ylm_wan[ir] = bs2 * (Ylm(0, 0, r_rot) - Ylm(1, 0, r_rot));
            }
            if(Wan_proj[iw].l == -2) // sp2
            {
                if(Wan_proj[iw].m == 0) ylm_wan[ir] = bs3 * Ylm(0,0, r_rot) - bs6 * Ylm(1,0, r_rot) + bs2 * Ylm(1,1,r_rot);
                if(Wan_proj[iw].m == 1) ylm_wan[ir] = bs3 * Ylm(0,0, r_rot) - bs6 * Ylm(1,0, r_rot) - bs2 * Ylm(1,1,r_rot);
                if(Wan_proj[iw].m == 2) ylm_wan[ir] = bs3 * Ylm(0,0, r_rot) +2.0*bs6 * Ylm(1,0, r_rot);
            }
            if(Wan_proj[iw].l == -3) // sp3
            {
                if(Wan_proj[iw].m == 0) ylm_wan[ir] = 0.5 * (Ylm(0,0, r_rot)+Ylm(1,0, r_rot)+Ylm(1,1,r_rot)+Ylm(1,2,r_rot));
                if(Wan_proj[iw].m == 1) ylm_wan[ir] = 0.5 * (Ylm(0,0, r_rot)+Ylm(1,0, r_rot)-Ylm(1,1,r_rot)-Ylm(1,2,r_rot));
                if(Wan_proj[iw].m == 2) ylm_wan[ir] = 0.5 * (Ylm(0,0, r_rot)-Ylm(1,0, r_rot)+Ylm(1,1,r_rot)-Ylm(1,2,r_rot));
                if(Wan_proj[iw].m == 3) ylm_wan[ir] = 0.5 * (Ylm(0,0, r_rot)-Ylm(1,0, r_rot)-Ylm(1,1,r_rot)+Ylm(1,2,r_rot));
            }

            if(Wan_proj[iw].l == -4) // sp3d
            {
                if(Wan_proj[iw].m == 0) ylm_wan[ir] = bs3 * Ylm(0,0, r_rot) - bs6 * Ylm(1,0, r_rot) + bs2 * Ylm(1,1,r_rot);
                if(Wan_proj[iw].m == 1) ylm_wan[ir] = bs3 * Ylm(0,0, r_rot) - bs6 * Ylm(1,0, r_rot) - bs2 * Ylm(1,1,r_rot);
                if(Wan_proj[iw].m == 2) ylm_wan[ir] = bs3 * Ylm(0,0, r_rot) +2.0*bs6 * Ylm(1,0, r_rot);
                if(Wan_proj[iw].m == 3) ylm_wan[ir] = bs2 * (Ylm(1,2,r_rot) + Ylm(2,2, r_rot));
                if(Wan_proj[iw].m == 4) ylm_wan[ir] = bs2 * (-Ylm(1,2,r_rot) + Ylm(2,2, r_rot));
            }
            if(Wan_proj[iw].l == -5) // sp3d
            {
                if(Wan_proj[iw].m == 0) ylm_wan[ir] = bs6 * Ylm(0,0, r_rot) - bs2 * Ylm(1,0, r_rot) - bs12 * Ylm(2,2,r_rot) + 0.5 * Ylm(2,4, r_rot);
                if(Wan_proj[iw].m == 1) ylm_wan[ir] = bs6 * Ylm(0,0, r_rot) + bs2 * Ylm(1,0, r_rot) - bs12 * Ylm(2,2,r_rot) + 0.5 * Ylm(2,4, r_rot);
                if(Wan_proj[iw].m == 2) ylm_wan[ir] = bs6 * Ylm(0,0, r_rot) - bs2 * Ylm(1,1, r_rot) - bs12 * Ylm(2,2,r_rot) - 0.5 * Ylm(2,4, r_rot);
                if(Wan_proj[iw].m == 3) ylm_wan[ir] = bs6 * Ylm(0,0, r_rot) + bs2 * Ylm(1,1, r_rot) - bs12 * Ylm(2,2,r_rot) - 0.5 * Ylm(2,4, r_rot);
                if(Wan_proj[iw].m == 4) ylm_wan[ir] = bs6 * Ylm(0,0, r_rot) - bs2 * Ylm(1,2, r_rot) + bs3 * Ylm(2,2,r_rot);
                if(Wan_proj[iw].m == 5) ylm_wan[ir] = bs6 * Ylm(0,0, r_rot) + bs2 * Ylm(1,2, r_rot) + bs3 * Ylm(2,2,r_rot);
            }
        }

        for(int lm = 0; lm < tot_LM; lm++)
        {
            csph[iw][lm] = 0.0;
            for(int ir = 0; ir < tot_LM; ir++)
            {
                csph[iw][lm] += ylm_invert[lm * tot_LM + ir] * ylm_wan[ir];
            }
        }

    }

    delete RT1;
    RT1 = new RmgTimer("7-Wannier: Amn: gf: radial_q");
    double tpiba = 2.0 * PI / Rmg_L.celldm[0];
    double tpiba2 = tpiba * tpiba;
    double gcut = sqrt(coarse_pwaves->gcut*tpiba2);

    delta_g = gcut/(gnum-1);

    radialfunc_g.resize(boost::extents[zona_list.size()][lmax+1][gnum]); 

    for(int izona = 0; izona < (int)zona_list.size(); izona++)
    {
        for(int l = 0; l <= lmax; l++)
        {
            InitRadialfunc_Gspace(l, zona_list[izona].first, zona_list[izona].second,  &radialfunc_g[izona][l][0], gnum, delta_g);

        }
    }

    delete RT1;
}

template  void Wannier<double>::GuideFunc(int kpt, double *guidefunc);
template  void Wannier<std::complex<double>>::GuideFunc(int kpt, std::complex<double> *guidefunc);
template <class T> void Wannier<T>::GuideFunc(int kpt, T *guidefunc)
{

    double tol = 1.0e-5;
    double tpiba = 2.0 * PI / Rmg_L.celldm[0];
    double ax[3];

    double kvec[3];
    std::complex<double> *guidefunc_g = new std::complex<double>[nbasis];

    RmgTimer *RT1;
    kvec[0] = ct.klist.k_all_cart[kpt][0];
    kvec[1] = ct.klist.k_all_cart[kpt][1];
    kvec[2] = ct.klist.k_all_cart[kpt][2];


    for(int idx = 0; idx < n_wannier * nbasis_noncoll; idx++) guidefunc[idx] = 0.0;
    std::complex<double> I_t(0.0, 1.0);
    std::complex<double> *gf = new std::complex<double>[nbasis];
    for(int iw = 0; iw < (int)Wan_proj.size(); iw++)
    {
        for(int ig = 0; ig < nbasis; ig++) guidefunc_g[ig] = 0.0;
        int izona = zona_index[iw];

        RT1 = new RmgTimer("7-Wannier: Amn: gf: guide_q");
        std::complex<double> phase;
        for(int l = 0; l <= lmax; l++)
        {
            for(int m = 0; m < 2*l + 1; m++)
            {
                int lm = l*l + m;
                if(std::abs(csph[iw][lm]) < tol ) continue;

                for(int ig = 0;ig < nbasis;ig++)
                {
                    ax[0] = coarse_pwaves->g[ig].a[0] * tpiba;
                    ax[1] = coarse_pwaves->g[ig].a[1] * tpiba;
                    ax[2] = coarse_pwaves->g[ig].a[2] * tpiba;

                    ax[0] += kvec[0];
                    ax[1] += kvec[1];
                    ax[2] += kvec[2];
                    double kr = ax[0] * Wan_proj[iw].center_cart[0] + ax[1] * Wan_proj[iw].center_cart[1] + ax[2] * Wan_proj[iw].center_cart[2]; 
                    double gval = std::sqrt(ax[0] * ax[0] + ax[1] * ax[1] + ax[2] * ax[2]);
                    int g_index = (int)(gval/delta_g);
                    if(g_index >= gnum-1) continue;
                    double frac = gval/delta_g - g_index;
                    std::complex<double> radial_part = radialfunc_g[izona][l][g_index] * frac + radialfunc_g[izona][l][g_index+1] * (1.0-frac) ;
                    phase = std::exp( std::complex<double>(0.0, kr));
                    guidefunc_g[ig] += std::pow(-I_t, l) * radial_part * Ylm(l, m, ax) * csph[iw][lm] * phase;
                }

            }
        }

        delete RT1;

        RT1 = new RmgTimer("7-Wannier: Amn: gf: fft");
        coarse_pwaves->FftInverse(guidefunc_g, gf);
        delete RT1;

        std::complex<double> *gf_C = (std::complex<double> *)&guidefunc[iw*nbasis_noncoll];
        double *gf_R = (double *)&guidefunc[iw*nbasis_noncoll];
        if(ct.noncoll)
        {
            std::complex<double> frac_up = 1.0, frac_dn = 1.0;
            if(Wan_proj[iw].spin == 1) 
            {
                frac_up = Wan_proj[iw].spin_dir[2];
                frac_dn = std::complex<double>(Wan_proj[iw].spin_dir[0], Wan_proj[iw].spin_dir[1]);
            }
            else
            {
                frac_dn = Wan_proj[iw].spin_dir[2];
                frac_up = std::complex<double>(Wan_proj[iw].spin_dir[0], Wan_proj[iw].spin_dir[1]);
            }

            for(int idx = 0; idx < nbasis; idx++) 
            {
                gf_C[idx] += frac_up * gf[idx];
                gf_C[nbasis + idx] += frac_dn * gf[idx];
            }
        }
        else if(ct.is_gamma)
        {
            for(int idx = 0; idx < nbasis; idx++) gf_R[idx] = std::real(gf[idx]);
        }
        else
        {
            for(int idx = 0; idx < nbasis; idx++) gf_C[idx] = gf[idx];
        }

        double norm_coeff = 0.0;

        RT1 = new RmgTimer("7-Wannier: Amn: gf: norm");
        for(int idx = 0; idx < nbasis_noncoll; idx++) norm_coeff += std::norm(guidefunc[iw*nbasis_noncoll + idx]);
        double vel = L.get_omega() / ((double)ngrid);
        MPI_Allreduce(MPI_IN_PLACE, &norm_coeff, 1, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
        norm_coeff = 1.0/std::sqrt(norm_coeff * vel);
        for(int idx = 0; idx < nbasis_noncoll; idx++) guidefunc[iw*nbasis_noncoll + idx] *= norm_coeff;
        delete RT1;

    }
    delete [] gf;
    delete [] guidefunc_g;

}

template void Wannier<double>::SetAmn_proj();
template void Wannier<std::complex<double>>::SetAmn_proj();
template <class T> void Wannier<T>::SetAmn_proj()
{

    int num_q = ct.klist.num_k_all;

    RmgTimer *RT1;
    size_t length = nstates * nbasis_noncoll * sizeof(T);
    T *psi_k = (T *)RmgMallocHost(length);
    T *Amn = (T *)RmgMallocHost(num_q * n_wannier * nstates * sizeof(T));

    for(int idx = 0; idx < num_q * n_wannier * nstates; idx++) Amn[idx] = 0.0;
    double vel = L.get_omega() / ((double)ngrid);
    T alpha(vel);
    T beta = 0.0;

    InitGuideFunc();
    T *guidefunc = new T[n_wannier * nbasis_noncoll];
    int num_tot_proj = 0;
    for(size_t ion = 0; ion < Atoms.size(); ion++)
    {
        SPECIES &AtomType = Species[Atoms[ion].species];
        num_tot_proj += AtomType.nh;
    }
    for(int ik = 0; ik < num_q; ik++)
    {
        RT1 = new RmgTimer("7-Wannier: Amn: gf");
        GuideFunc(ik, guidefunc);
        delete RT1;
        int ik_irr = ct.klist.k_map_index[ik];
        int isym = ct.klist.k_map_symm[ik];
        int isyma = std::abs(isym) -1;

        RT1 = new RmgTimer("7-Wannier: Amn: read and rotate");
        if(isym == 1 ) {
            ReadPsiFromSingleFile(ik_irr, wavefile, psi_k);
        }
        else {
            ReadRotatePsi(ik_irr, isym, isyma, wavefile, psi_k);
        }
        delete RT1;

        RT1 = new RmgTimer("7-Wannier: Amn: gemm");
        RmgGemm("C", "N", nstates, n_wannier, nbasis_noncoll, alpha, psi_k, nbasis_noncoll, guidefunc, nbasis_noncoll,
                beta, &Amn[ik*n_wannier*nstates], nstates);
        delete RT1;

        RT1 = new RmgTimer("7-Wannier: Amn: us");
        if(!ct.norm_conserving_pp)
        {
            T *Nlweight_k = (T *)RmgMallocHost(num_tot_proj * nbasis* sizeof(T));
            RmgTimer *RT2 = new RmgTimer("7-Wannier: Amn: us: read NL");
            std::string filename;
            filename = "PROJECTORS/NLprojectors_kpt" + std::to_string(ik);
            std::complex<double> *Nlweight_C = (std::complex<double> *)Nlweight_k;
            ReadNlweight(filename, num_tot_proj, Nlweight_C);

            delete RT2;

            std::complex<double> *qq_dk_one, *qq_dk_so_one;
            qq_dk_one = new std::complex<double>[Atoms.size() * ct.max_nl * ct.max_nl];
            qq_dk_so_one = new std::complex<double>[4*Atoms.size() * ct.max_nl * ct.max_nl];
            for(size_t ion = 0; ion < Atoms.size(); ion++)
            {
                int nh = Species[Atoms[ion].species].nh;
                for(int ih = 0; ih < nh*nh; ih++) 
                {
                    qq_dk_one[ion * ct.max_nl * ct.max_nl + ih] = Atoms[ion].qqq[ih];
                }
                if(ct.noncoll)
                {
                    for(int ih = 0; ih < 4*nh*nh; ih++) 
                    {
                        qq_dk_so_one[ion * ct.max_nl * ct.max_nl * 4 + ih] = Atoms[ion].qqq_so[ih];
                    }
                }
            }
            int ikn = -1;
            Mmn_us(ik, ikn, psi_k, nstates, guidefunc, n_wannier, &Amn[ik*nstates*n_wannier], qq_dk_one, qq_dk_so_one,
                    num_tot_proj, Nlweight_k, Nlweight_k);
            delete [] qq_dk_one;
            delete [] qq_dk_so_one;

            RmgFreeHost(Nlweight_k);

        }
        delete RT1;

    }


    RT1 = new RmgTimer("7-Wannier: Amn: Reduce");
    int count = num_q * n_wannier * nstates;
    MPI_Allreduce(MPI_IN_PLACE, Amn, count, MPI_DOUBLE_COMPLEX, MPI_SUM, pct.grid_comm);
    delete RT1;

    RT1 = new RmgTimer("7-Wannier: Amn: write");
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
        fclose(famn);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    delete RT1;

}


template  void Wannier<double>::ReadPsiFromSingleFile(int ikindex, std::string wavefile, double *psi_k);
template  void Wannier<std::complex<double>>::ReadPsiFromSingleFile(int ikindex, std::string wavefile, std::complex<double> *psi_k);
template <class T> void Wannier<T>::ReadPsiFromSingleFile(int ikindex, std::string wavefile, T *psi_k)
{
    // Write the domain distributed wavefunction array and map it to psi_s
    MPI_Datatype wftype = MPI_DOUBLE;
    if(typeid(T) == typeid(std::complex<double>)) wftype = MPI_DOUBLE_COMPLEX;

    int sizes_c[4];
    int subsizes_c[4];
    int starts_c[4];

    sizes_c[0] = nstates;
    sizes_c[1] = G.get_NX_GRID(1);
    sizes_c[2] = G.get_NY_GRID(1);
    sizes_c[3] = G.get_NZ_GRID(1);

    subsizes_c[0] = nstates;
    subsizes_c[1] = G.get_PX0_GRID(1);
    subsizes_c[2] = G.get_PY0_GRID(1);
    subsizes_c[3] = G.get_PZ0_GRID(1);

    starts_c[0] = 0;
    starts_c[1] = G.get_PX_OFFSET(1);
    starts_c[2] = G.get_PY_OFFSET(1);
    starts_c[3] = G.get_PZ_OFFSET(1);

    int order = MPI_ORDER_C;
    MPI_Info fileinfo;
    MPI_Datatype grid_c;
    MPI_Status status;

    MPI_Type_create_subarray(4, sizes_c, subsizes_c, starts_c, order, wftype, &grid_c);
    MPI_Type_commit(&grid_c);

    MPI_Info_create(&fileinfo);

    int amode = MPI_MODE_RDWR;
    MPI_File mpi_fhand ;

    MPI_Barrier(G.comm);

    std::string filename = wavefile + "_spin"+std::to_string(pct.spinpe) + "_kpt" + std::to_string(ikindex);
    MPI_File_open(G.comm, filename.c_str(), amode, fileinfo, &mpi_fhand);
    MPI_Offset disp = 0;

    int dis_dim = G.get_P0_BASIS(1);
    MPI_File_set_view(mpi_fhand, disp, wftype, grid_c, "native", MPI_INFO_NULL);

    dis_dim = dis_dim * nstates;
    MPI_File_read_all(mpi_fhand, psi_k, dis_dim, wftype, &status);
    if(ct.noncoll) {
        MPI_File_read_all(mpi_fhand, &psi_k[dis_dim], dis_dim, wftype, &status);
    }

    MPI_Barrier(G.comm);
    MPI_File_close(&mpi_fhand);

    MPI_Type_free(&grid_c);
    fflush(NULL);
    MPI_Barrier(G.comm);

}

