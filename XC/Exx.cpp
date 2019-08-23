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


#include "Exx.h"
#include "RmgTimer.h"
#include "RmgException.h"

// This class implements exact exchange for delocalized orbitals.
// The wavefunctions are stored in a single file and are not domain
// decomposed. The file name is given by wavefile_in which is
// mmapped to an array. We only access the array in read-only mode.


template Exx<double>::Exx(BaseGrid &, Lattice &, std::string &, int, double *);
template Exx<std::complex<double>>::Exx(BaseGrid &, Lattice &, std::string &, int, std::complex<double> *);

template <class T> Exx<T>::Exx (BaseGrid &G_in,
          Lattice &L_in,
          std::string &wavefile_in,
          int nstates_in,
          T *psi_in) : G(G_in), L(L_in), wavefile(wavefile_in), nstates(nstates_in), psi(psi_in)
{
    RmgTimer RT0("5-Functional: Exact Exchange");

    pbasis = G.get_P0_BASIS(G.default_FG_RATIO);
    N = (size_t)G.get_NX_GRID(G.default_FG_RATIO) *
              (size_t)G.get_NY_GRID(G.default_FG_RATIO) *
              (size_t)G.get_NZ_GRID(G.default_FG_RATIO);

    int ratio = G.default_FG_RATIO;

    // Write the domain distributed wavefunction array and map it to psi_s
    MPI_Datatype wftype = MPI_DOUBLE;
    if(typeid(T) == typeid(std::complex<double>)) wftype = MPI_DOUBLE_COMPLEX;

    int sizes_c[3], sizes_f[3];
    int subsizes_c[3], subsizes_f[3];
    int starts_c[3], starts_f[3];

    sizes_c[0] = G.get_NX_GRID(1);
    sizes_c[1] = G.get_NY_GRID(1);
    sizes_c[2] = G.get_NZ_GRID(1);

    sizes_f[0] = G.get_NX_GRID(ratio);
    sizes_f[1] = G.get_NY_GRID(ratio);
    sizes_f[2] = G.get_NZ_GRID(ratio);

    subsizes_c[0] = G.get_PX0_GRID(1);
    subsizes_c[1] = G.get_PY0_GRID(1);
    subsizes_c[2] = G.get_PZ0_GRID(1);

    subsizes_f[0] = G.get_PX0_GRID(ratio);
    subsizes_f[1] = G.get_PY0_GRID(ratio);
    subsizes_f[2] = G.get_PZ0_GRID(ratio);

    starts_c[0] = G.get_PX_OFFSET(1);
    starts_c[1] = G.get_PY_OFFSET(1);
    starts_c[2] = G.get_PZ_OFFSET(1);

    starts_f[0] = G.get_PX_OFFSET(ratio);
    starts_f[1] = G.get_PY_OFFSET(ratio);
    starts_f[2] = G.get_PZ_OFFSET(ratio);

    int order = MPI_ORDER_C;
    MPI_Info fileinfo;
    MPI_Datatype grid_c, grid_f;
    MPI_Status status;


    MPI_Type_create_subarray(3, sizes_c, subsizes_c, starts_c, order, wftype, &grid_c);
    MPI_Type_commit(&grid_c);

    MPI_Type_create_subarray(3, sizes_f, subsizes_f, starts_f, order, MPI_DOUBLE, &grid_f);
    MPI_Type_commit(&grid_f);

    MPI_Info_create(&fileinfo);

    int amode = MPI_MODE_RDWR|MPI_MODE_CREATE;
    MPI_File mpi_fhand ;

    MPI_Barrier(G.comm);

    MPI_File_open(G.comm, wavefile.c_str(), amode, fileinfo, &mpi_fhand);
    MPI_Offset disp = 0;
    T *wfptr = psi;
    MPI_File_set_view(mpi_fhand, disp, wftype, grid_c, "native", MPI_INFO_NULL);
    for(int st=0;st < nstates;st++)
    {
        MPI_File_write_all(mpi_fhand, wfptr, pbasis, MPI_DOUBLE, &status);
        wfptr += pbasis;
    }
    MPI_File_close(&mpi_fhand);
    fflush(NULL);
    MPI_Barrier(G.comm);

    serial_fd = open(wavefile.c_str(), O_RDONLY, (mode_t)0600);
    if(serial_fd < 0)
        throw RmgFatalException() << "Error! Could not open " << wavefile << " . Terminating.\n";

    size_t length = nstates * N;
    psi_s = (T *)mmap(NULL, length, PROT_READ, MAP_PRIVATE, serial_fd, 0);
    MPI_Barrier(G.comm);

    MPI_Type_free(&grid_f);
    MPI_Type_free(&grid_c);

    // Compute st_start and st_count
    st_start = 0;
    st_count = nstates;
    int rank = G.get_rank();
    int NPES = G.get_NPES(); 
    if(NPES > 1)
    {
        st_count = nstates / NPES;
        st_start = st_count * rank;
        int rem = nstates % NPES;
        if(rank < rem)
        {
            st_count++;
            st_start += rank;
        }
        else
        {
            st_start += rem;
        }

    }
}


template <class T> Exx<T>::~Exx(void)
{
    close(serial_fd);
    size_t length = nstates * N;
    munmap(psi_s, length);
    unlink(wavefile.c_str());
}
