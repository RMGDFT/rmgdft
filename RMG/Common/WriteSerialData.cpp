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
#if !(defined(_WIN32) || defined(_WIN64))
    #include <unistd.h>
#else
    #include <io.h>
#endif
#include <complex>
#include "const.h"
#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "State.h"
#include "Kpoint.h"
#include "transition.h"


/* 

  This routine writes the charge density, hartree and exchange correlation potentials and
  the electronic orbitals to serial files where each global grid object represents one file.

*/

template void WriteSerialData (std::string&, double *, double *, double *, Kpoint<double> **);
template void WriteSerialData (std::string&, double *, double *, double *, Kpoint<std::complex<double> > **);


template <typename KpointType>
void WriteSerialData (std::string& name, double * vh, double * rho, double * vxc, Kpoint<KpointType> ** Kptr)
{
    int sizes_c[3], sizes_f[3];
    int subsizes_c[3], subsizes_f[3];
    int starts_c[3], starts_f[3];
    int ratio = Kptr[0]->G->default_FG_RATIO;
    int pbasis = Kptr[0]->G->get_P0_BASIS(1);
    int fpbasis = Kptr[0]->G->get_P0_BASIS(ratio);
    int pex, pey, pez;
    pe2xyz (pct.gridpe, &pex, &pey, &pez);

    MPI_Info fileinfo;
    MPI_Datatype grid_c, grid_f;
    MPI_Status status;
    MPI_Offset disp;

    MPI_Datatype wftype = MPI_DOUBLE;
    if(typeid(KpointType) == typeid(std::complex<double>)) wftype = MPI_DOUBLE_COMPLEX;

    sizes_c[0] = Kptr[0]->G->get_NX_GRID(1);
    sizes_c[1] = Kptr[0]->G->get_NY_GRID(1);
    sizes_c[2] = Kptr[0]->G->get_NZ_GRID(1);

    sizes_f[0] = Kptr[0]->G->get_NX_GRID(ratio);
    sizes_f[1] = Kptr[0]->G->get_NY_GRID(ratio);
    sizes_f[2] = Kptr[0]->G->get_NZ_GRID(ratio);

    subsizes_c[0] = Kptr[0]->G->get_PX0_GRID(1);
    subsizes_c[1] = Kptr[0]->G->get_PY0_GRID(1);
    subsizes_c[2] = Kptr[0]->G->get_PZ0_GRID(1);

    subsizes_f[0] = Kptr[0]->G->get_PX0_GRID(ratio);
    subsizes_f[1] = Kptr[0]->G->get_PY0_GRID(ratio);
    subsizes_f[2] = Kptr[0]->G->get_PZ0_GRID(ratio);

    starts_c[0] = Kptr[0]->G->get_PX_OFFSET(1);
    starts_c[1] = Kptr[0]->G->get_PY_OFFSET(1);
    starts_c[2] = Kptr[0]->G->get_PZ_OFFSET(1);

    starts_f[0] = Kptr[0]->G->get_PX_OFFSET(ratio);
    starts_f[1] = Kptr[0]->G->get_PY_OFFSET(ratio);
    starts_f[2] = Kptr[0]->G->get_PZ_OFFSET(ratio);


    int order = MPI_ORDER_C;
    MPI_Type_create_subarray(3, sizes_c, subsizes_c, starts_c, order, wftype, &grid_c);
    MPI_Type_commit(&grid_c);

    MPI_Type_create_subarray(3, sizes_f, subsizes_f, starts_f, order, MPI_DOUBLE, &grid_f);
    MPI_Type_commit(&grid_f);

    MPI_Info_create(&fileinfo);

    int amode = MPI_MODE_RDWR|MPI_MODE_CREATE;
    MPI_File mpi_fhand ;

    MPI_Barrier(pct.grid_comm);
    if(pct.spinpe == 0) {
        std::string newname = name + ".vh";
        MPI_File_open(pct.grid_comm, newname.c_str(), amode, fileinfo, &mpi_fhand);
        disp=0;
        MPI_File_set_view(mpi_fhand, disp, MPI_DOUBLE, grid_f, "native", MPI_INFO_NULL);
        MPI_File_write_all(mpi_fhand, vh, fpbasis, MPI_DOUBLE, &status);
        MPI_File_close(&mpi_fhand);

        MPI_Barrier(pct.grid_comm);
        newname = name + ".vxc";
        MPI_File_open(pct.grid_comm, newname.c_str(), amode, fileinfo, &mpi_fhand);
        disp=0;
        MPI_File_set_view(mpi_fhand, disp, MPI_DOUBLE, grid_f, "native", MPI_INFO_NULL);
        for(int ispin = 0; ispin < ct.nspin; ispin++) {
            MPI_File_write_all(mpi_fhand, vxc+ispin*fpbasis, fpbasis, MPI_DOUBLE, &status);
        }
        MPI_File_close(&mpi_fhand);

        MPI_Barrier(pct.grid_comm);
        newname = name + ".rho";
        MPI_File_open(pct.grid_comm, newname.c_str(), amode, fileinfo, &mpi_fhand);
        disp=0;
        MPI_File_set_view(mpi_fhand, disp, MPI_DOUBLE, grid_f, "native", MPI_INFO_NULL);
        for(int ispin = 0; ispin < ct.nspin; ispin++) {
            MPI_File_write_all(mpi_fhand, rho+ispin*fpbasis, fpbasis, MPI_DOUBLE, &status);
        }
        MPI_File_close(&mpi_fhand);
    }

    OrbitalHeader H;
    H.nx = sizes_c[0];
    H.ny = sizes_c[1];
    H.nz = sizes_c[2];

    for (int ik = 0; ik < ct.num_kpts_pe; ik++)
    {
        for (int is = 0; is < Kptr[ik]->nstates; is++)
        {
            H.eig = Kptr[ik]->Kstates[is].eig[0];
            H.occ = Kptr[ik]->Kstates[is].occupation[0];

            std::string wfname = name + "_spin" + std::to_string(pct.spinpe) +
                "_kpt" + std::to_string(pct.kstart+ik) +
                "_wf" + std::to_string(is);

            if(pct.gridpe == 0)
            {
                int fhand = open(wfname.c_str(), O_CREAT | O_TRUNC | O_RDWR, S_IREAD | S_IWRITE);
                if (fhand < 0) {
                    rmg_printf("Can't open restart file %s", wfname.c_str());
                    rmg_error_handler(__FILE__, __LINE__, "Terminating.");
                }
                size_t wsize = write (fhand, &H, sizeof(OrbitalHeader));
                if(wsize != sizeof(OrbitalHeader))
                    rmg_error_handler (__FILE__,__LINE__,"error writing");
                close(fhand);
                fflush(NULL);
            }

            MPI_Barrier(pct.grid_comm);
            MPI_File_open(pct.grid_comm, wfname.c_str(), amode, fileinfo, &mpi_fhand);
            disp = sizeof(OrbitalHeader);
            KpointType *wfptr = Kptr[ik]->Kstates[is].psi;
            MPI_File_set_view(mpi_fhand, disp, wftype, grid_c, "native", MPI_INFO_NULL);
            MPI_File_write_all(mpi_fhand, wfptr, pbasis, wftype, &status);
            MPI_File_close(&mpi_fhand);
        }
    }

    MPI_Type_free(&grid_f);
    MPI_Type_free(&grid_c);

} // WriteSerialData

