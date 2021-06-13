/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
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

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>
#include <complex>
#if !(defined(_WIN32) || defined(_WIN64))
#include <unistd.h>
#else
#include <io.h>
#endif
#include <time.h>
#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "State.h"
#include "Kpoint.h"
#include "transition.h"


template void Write_Wfs_forWannier (int kpt, Kpoint<double> * kptr, std::vector<bool> exclude_bands, std::string wavefile);
template void Write_Wfs_forWannier (int kpt, Kpoint<std::complex<double> > * kptr, std::vector<bool> exclude_bands, std::string wavefile);
    template <typename T>
void Write_Wfs_forWannier (int kpt, Kpoint<T> * kptr, std::vector<bool> exclude_bands, std::string wavefile)
{
    mkdir ("WfsForWannier90", S_IRWXU);
    std::string wf_filename = wavefile + "_spin" + std::to_string(pct.spinpe);
    wf_filename += "_kpt"+  std::to_string(kpt);

    // Write the domain distributed wavefunction array in a single file.
    MPI_Datatype wftype = MPI_DOUBLE;
    if(typeid(T) == typeid(std::complex<double>)) wftype = MPI_DOUBLE_COMPLEX;

    int sizes_c[4];
    int subsizes_c[4];
    int starts_c[4];

    int nstates_tot = kptr->nstates;

    int dis_dim = Rmg_G->get_P0_BASIS(1);
    int dis_dim_noncoll = dis_dim * ct.noncoll_factor;
    int num_wan_st = nstates_tot;
    for(int st = 0; st < nstates_tot; st++) {
        if(exclude_bands[st]) num_wan_st--;
    }

    std::vector<T> psi;
    psi.resize(num_wan_st * dis_dim_noncoll);
    int st_w = -1;
    for(int st = 0; st < nstates_tot; st++) {

        if(exclude_bands[st]) continue;
        st_w++;
        for(int idx = 0; idx < dis_dim_noncoll; idx++) {
            psi[st_w * dis_dim_noncoll + idx] = kptr->orbital_storage[st* dis_dim_noncoll + idx];
        }

    }
    sizes_c[0] = num_wan_st;
    sizes_c[1] = Rmg_G->get_NX_GRID(1);
    sizes_c[2] = Rmg_G->get_NY_GRID(1);
    sizes_c[3] = Rmg_G->get_NZ_GRID(1);

    subsizes_c[0] = num_wan_st;
    subsizes_c[1] = Rmg_G->get_PX0_GRID(1);
    subsizes_c[2] = Rmg_G->get_PY0_GRID(1);
    subsizes_c[3] = Rmg_G->get_PZ0_GRID(1);

    starts_c[0] = 0;
    starts_c[1] = Rmg_G->get_PX_OFFSET(1);
    starts_c[2] = Rmg_G->get_PY_OFFSET(1);
    starts_c[3] = Rmg_G->get_PZ_OFFSET(1);

    int order = MPI_ORDER_C;
    MPI_Info fileinfo;
    MPI_Datatype grid_c;
    MPI_Status status;

    MPI_Type_create_subarray(4, sizes_c, subsizes_c, starts_c, order, wftype, &grid_c);
    MPI_Type_commit(&grid_c);

    MPI_Info_create(&fileinfo);

    int amode = MPI_MODE_RDWR|MPI_MODE_CREATE;
    MPI_File mpi_fhand ;

    MPI_Barrier(Rmg_G->comm);

    MPI_File_open(Rmg_G->comm, wf_filename.c_str(), amode, fileinfo, &mpi_fhand);
    MPI_Offset disp = 0;

    MPI_File_set_view(mpi_fhand, disp, wftype, grid_c, "native", MPI_INFO_NULL);

    dis_dim = dis_dim * num_wan_st;
    MPI_File_write_all(mpi_fhand, psi.data(), dis_dim, wftype, &status);
    if(ct.noncoll) {
        MPI_File_write_all(mpi_fhand, psi.data()+dis_dim, dis_dim, wftype, &status);
    }

    MPI_Barrier(Rmg_G->comm);
    MPI_File_close(&mpi_fhand);

    MPI_Type_free(&grid_c);
    MPI_Barrier(Rmg_G->comm);
}
