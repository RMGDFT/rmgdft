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
#include "transition.h"
#include "LocalObject.h"
#include "prototypes_on.h"
#include "blas.h"



/* 

  This routine writes the charge density, hartree and exchange correlation potentials and
  the electronic orbitals to serial files where each global grid object represents one file.

*/

template void WriteWavefunctions (std::string&, LocalObject<double> &LO , double *Cij, BaseGrid &BG, 
        double *eig, double *occ);
//template void WriteWavefunctions (std::string&, LocalObjectt<std::complex<double> > &LO, std::complex<double> &Cij);

template <typename KpointType>
void WriteWavefunctions (std::string& name, LocalObject<KpointType> &Phi, KpointType *Cij_dist, BaseGrid &BG, 
        double *eig, double *occ)
{
    int sizes_c[3];
    int subsizes_c[3];
    int starts_c[3];
    int pbasis = BG.get_P0_BASIS(1);
    int pex, pey, pez;
    pe2xyz (pct.gridpe, &pex, &pey, &pez);

    MPI_Info fileinfo;
    MPI_Datatype grid_c;
    MPI_Status status;
    MPI_Offset disp;

    MPI_Datatype wftype = MPI_DOUBLE;
    if(typeid(KpointType) == typeid(std::complex<double>)) wftype = MPI_DOUBLE_COMPLEX;

    sizes_c[0] = BG.get_NX_GRID(1);
    sizes_c[1] = BG.get_NY_GRID(1);
    sizes_c[2] = BG.get_NZ_GRID(1);

    subsizes_c[0] = BG.get_PX0_GRID(1);
    subsizes_c[1] = BG.get_PY0_GRID(1);
    subsizes_c[2] = BG.get_PZ0_GRID(1);

    starts_c[0] = BG.get_PX_OFFSET(1);
    starts_c[1] = BG.get_PY_OFFSET(1);
    starts_c[2] = BG.get_PZ_OFFSET(1);


    int order = MPI_ORDER_C;
    MPI_Type_create_subarray(3, sizes_c, subsizes_c, starts_c, order, wftype, &grid_c);
    MPI_Type_commit(&grid_c);

    MPI_Info_create(&fileinfo);

    int amode = MPI_MODE_RDWR|MPI_MODE_CREATE;
    MPI_File mpi_fhand ;

    OrbitalHeader H;
    H.nx = sizes_c[0];
    H.ny = sizes_c[1];
    H.nz = sizes_c[2];

//  calculate the wavefuctions from localized orbitals
    KpointType *psi = new KpointType[Phi.num_tot * pbasis]();
    KpointType *Cij_global = new KpointType[Phi.num_tot * Phi.num_tot]();
    KpointType *Cij_local = new KpointType[Phi.num_tot * Phi.num_thispe]();

    mat_dist_to_global(Cij_dist, pct.desca, Cij_global);

    for(int i = 0; i < Phi.num_thispe; i++)
    {
        int i_glob = Phi.index_proj_to_global[i];
        if(i_glob < 0) continue;
        for(int j = 0; j < Phi.num_tot; j++)
            Cij_local[j * Phi.num_thispe + i] = Cij_global[j * Phi.num_tot + i_glob];
    }

    
    double one(1.0), zero(0.0);
    dgemm("N", "N", &pbasis, &Phi.num_tot, &Phi.num_thispe, &one, Phi.storage_proj, &pbasis,
                Cij_local, &Phi.num_thispe, &zero, psi, &pbasis);
    
//    for (int ik = 0; ik < ct.num_kpts_pe; ik++)
    int ik = 0;
    {
        for (int is = 0; is < ct.num_states; is++)
        {
            H.eig = eig[is];
            H.occ = occ[is];

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
            KpointType *wfptr = &psi[is * pbasis];
            MPI_File_set_view(mpi_fhand, disp, wftype, grid_c, "native", MPI_INFO_NULL);
            MPI_File_write_all(mpi_fhand, wfptr, pbasis, MPI_DOUBLE, &status);
            MPI_File_close(&mpi_fhand);
        }
    }

    MPI_Type_free(&grid_c);
    delete [] psi;
    delete [] Cij_global;
    delete [] Cij_local;

} // WriteSerialData

