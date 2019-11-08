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


#include "const.h"
#include "Exxbase.h"
#include "RmgTimer.h"
#include "RmgException.h"
#include "RmgGemm.h"
#include "transition.h"
#include "rmgtypedefs.h"
#include "pe_control.h"

// This class implements exact exchange for delocalized orbitals.
// The wavefunctions are stored in a single file and are not domain
// decomposed. The file name is given by wavefile_in which is
// mmapped to an array. We only access the array in read-only mode.


template Exxbase<double>::Exxbase(BaseGrid &, Lattice &, const std::string &, int, double *, double *, int);
template Exxbase<std::complex<double>>::Exxbase(BaseGrid &, Lattice &, const std::string &, int, double *, std::complex<double> *, int );

template Exxbase<double>::~Exxbase(void);
template Exxbase<std::complex<double>>::~Exxbase(void);

template <class T> Exxbase<T>::Exxbase (
          BaseGrid &G_in,
          Lattice &L_in,
          const std::string &wavefile_in,
          int nstates_in,
          double *occ_in,
          T *psi_in, int mode_in) : G(G_in), L(L_in), wavefile(wavefile_in), nstates(nstates_in), init_occ(occ_in), psi(psi_in), mode(mode_in)
{
    RmgTimer RT0("5-Functional: Exx init");

    for(int st=0;st < nstates;st++) occ.push_back(init_occ[st]);
    tpiba = 2.0 * PI / L.celldm[0];
    tpiba2 = tpiba * tpiba;
    alpha = L.get_omega() / ((double)(G.get_NX_GRID(1) * G.get_NY_GRID(1) * G.get_NZ_GRID(1)));

    if(mode == EXX_DIST_FFT) 
    {
        pbasis = G.get_P0_BASIS(1);
        pwave = coarse_pwaves;
        psi_s = psi;
        LG = &G_in;
    }
    else
    {
        pbasis = G.get_NX_GRID(1) * G.get_NY_GRID(1) * G.get_NZ_GRID(1);
        LG = new BaseGrid(G.get_NX_GRID(1), G.get_NY_GRID(1), G.get_NZ_GRID(1), 1, 1, 1, 0, 1);
        int rank = G.get_rank();
        MPI_Comm_split(G.comm, rank+1, rank, &lcomm);
        LG->set_rank(0, lcomm);
        pwave = new Pw(*LG, L, 1, false);


    }

    setup_gfac();

}

template <> void Exxbase<double>::fftpair(double *psi_i, double *psi_j, std::complex<double> *p)
{
    std::complex<double> ZERO_t(0.0, 0.0);
    for(int idx=0;idx < pbasis;idx++) p[idx] = std::complex<double>(psi_i[idx] * psi_j[idx], 0.0);
    pwave->FftForward(p, p);
    for(int ig=0;ig < pbasis;ig++) p[ig] *= gfac[ig];
    pwave->FftInverse(p, p);
}

template <> void Exxbase<std::complex<double>>::fftpair(std::complex<double> *psi_i, std::complex<double> *psi_j, std::complex<double> *p)
{
}

// This implements different ways of handling the divergence at G=0
template <> void Exxbase<double>::setup_gfac(void)
{
    gfac = new double[pbasis];

    const std::string &dftname = Functional::get_dft_name_rmg();
    double scrlen = Functional::get_screening_parameter_rmg();

    if((dftname == "gaupbe") || (scr_type == GAU_SCREENING))
    {
        double gau_scrlen = Functional::get_gau_parameter_rmg();
        if (gau_scrlen < 1.0e-5) gau_scrlen = 0.15;
        double a0 = pow(PI / gau_scrlen, 1.5);
        for(int ig=0;ig < pbasis;ig++)
        {
            if(pwave->gmask[ig])
                gfac[ig] = a0 * exp(-tpiba2*pwave->gmags[ig] / 4.0 / gau_scrlen);
            else
                gfac[ig] = 0.0;
        }
    }
    else if(scr_type == ERF_SCREENING)
    {
        double a0 = 4.0 * PI;
        for(int ig=0;ig < pbasis;ig++)
        {
            if((pwave->gmags[ig] > 1.0e-6) && pwave->gmask[ig])
                gfac[ig] = a0 * exp(-tpiba2*pwave->gmags[ig] / 4.0 / (scrlen*scrlen)) / (tpiba2*pwave->gmags[ig]);
            else
                gfac[ig] = 0.0;
        }
    }
    else if(scr_type == ERFC_SCREENING)
    {
        double a0 = 4.0 * PI;
        for(int ig=0;ig < pbasis;ig++)
        {
            if((pwave->gmags[ig] > 1.0e-6) && pwave->gmask[ig])
                gfac[ig] = a0 * (1.0 - exp(-tpiba2*pwave->gmags[ig] / 4.0 / (scrlen*scrlen))) / (tpiba2*pwave->gmags[ig]);
            else
                gfac[ig] = 0.0;
        }
    }
    else  // Default is probably wrong
    {
        for(int ig=0;ig < pbasis;ig++)
        {
            if((pwave->gmags[ig] > 1.0e-6) && pwave->gmask[ig])
                gfac[ig] = 1.0/(pwave->gmags[ig] *tpiba2);
            else
                gfac[ig] = 0.0;
        }
    }

}

template <> void Exxbase<std::complex<double>>::setup_gfac(void)
{
}

// These compute the action of the exact exchange operator on all wavefunctions
// and writes the result into vfile.
template <> void Exxbase<double>::Vexx(double *vexx)
{
    RmgTimer RT0("5-Functional: Exx potential");
    double scale = - 1.0 / (double)pwave->global_basis;

    double ZERO_t(0.0), ONE_T(1.0);
    char *trans_a = "t";
    char *trans_n = "n";

    // Clear vexx
    for(int idx=0;idx < nstates*pbasis;idx++) vexx[idx] = 0.0;

    int nstates_occ = 0;
    for(int st=0;st < nstates;st++) if(occ[st] > 1.0e-6) nstates_occ++;
    MPI_Allreduce(MPI_IN_PLACE, &nstates_occ, 1, MPI_INT, MPI_MAX, G.comm);

    if(mode == EXX_DIST_FFT)
    {
        std::complex<double> *p = (std::complex<double> *)fftw_malloc(sizeof(std::complex<double>) * pbasis);
        // Loop over fft pairs and compute Kij(r) 
        for(int i=0;i < nstates_occ;i++)
        {
            double *psi_i = (double *)&psi_s[i*pbasis];
            for(int j=0;j < nstates_occ;j++)
            {   
                double *psi_j = (double *)&psi_s[j*pbasis];
                RmgTimer RT1("5-Functional: Exx potential fft");
                fftpair(psi_i, psi_j, p);
                //if(G.get_rank()==0)printf("TTTT  %d  %d  %e\n",i,j,std::real(p[1]));
                for(int idx = 0;idx < pbasis;idx++)vexx[i*pbasis +idx] += scale * std::real(p[idx]) * psi_s[j*pbasis + idx];
            }
        }

        fftw_free(p);
    }
    else
    {
        rmg_error_handler (__FILE__,__LINE__,"Exx potential mode not programmed yet. Terminating.");
    }

}

template <> double Exxbase<double>::Exxenergy(double *vexx)
{
    double energy = 0.0;
    for(int st=0;st < nstates;st++)
    {
        for(int i=0;i < pbasis;i++) energy += occ[st]*vexx[st*pbasis + i]*psi[st*pbasis + i];
    }
    energy = 0.5*energy * this->L.get_omega() / (double)this->G.get_GLOBAL_BASIS(1);
    MPI_Allreduce(MPI_IN_PLACE, &energy, 1, MPI_DOUBLE, MPI_SUM, this->G.comm);

    return energy;
}

template <> void Exxbase<std::complex<double>>::Vexx(std::complex<double> *vexx)
{
    RmgTimer RT0("5-Functional: Exx potential");
    rmg_error_handler (__FILE__,__LINE__,"Exx potential not programmed for non-gamma yet. Terminating.");
}

template <> double Exxbase<std::complex<double>>::Exxenergy(std::complex<double> *vexx)
{
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


    char *buf;
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
        if (ij_end > wf_pairs.size()) ij_end = wf_pairs.size();
        RmgTimer *RT1 = new RmgTimer("5-Functional: Exx: ij+fft");
        for(int ij = ij_start; ij < ij_end; ij++)
        {

            int ic = ij - ij_start;
            int i = wf_pairs[ij].first;
            int j = wf_pairs[ij].second;
            double *psi_i = (double *)&psi_s[i*pbasis];
            double *psi_j = (double *)&psi_s[j*pbasis];
            fftpair(psi_i, psi_j, wf_fft);

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
            if (kl_end > wf_pairs.size()) kl_end = wf_pairs.size();
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

    int sizes_c[3], sizes_f[3];
    int subsizes_c[3], subsizes_f[3];
    int starts_c[3], starts_f[3];

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

    int sizes_c[3], sizes_f[3];
    int subsizes_c[3], subsizes_f[3];
    int starts_c[3], starts_f[3];

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

    std::string filename = wavefile + "_spin"+std::to_string(pct.spinpe);
    MPI_File_open(G.comm, filename.c_str(), amode, fileinfo, &mpi_fhand);
    MPI_Offset disp = 0;

    T *wfptr;
    MPI_File_set_view(mpi_fhand, disp, wftype, grid_c, "native", MPI_INFO_NULL);
    int dis_dim = G.get_P0_BASIS(1);
    for(int st=0;st < nstates;st++)
    {
        wfptr = &psi[st * dis_dim];
        MPI_File_write_all(mpi_fhand, wfptr, dis_dim, MPI_DOUBLE, &status);
    }
    MPI_Barrier(G.comm);
    MPI_File_close(&mpi_fhand);
    MPI_Type_free(&grid_c);
    fflush(NULL);
    MPI_Barrier(G.comm);
}

