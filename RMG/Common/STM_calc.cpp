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


#include <float.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h> 
#include <unistd.h>
#include <unordered_map>
#include <csignal>
#include <sstream>
#include <iomanip>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <iterator>


#include "const.h"
#include "RmgTimer.h"
#include "RmgException.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "transition.h"
#include "Kpoint.h"
#include "InputKey.h"
#include "blas.h"
#include "RmgThread.h"
#include "rmgthreads.h"


#include "../Headers/common_prototypes.h"
#include "../Headers/common_prototypes1.h"
#include "prototypes_tddft.h"
#include "Exxbase.h"
#include "Neb.h"
#include "Wannier.h"
#include "GlobalSums.h"


void OutputSTM(std::vector<double> rho_2d, int NX, int NY, std::string filenamne);
template void STM_calc (Kpoint<double> **Kptr, double *rho, std::vector<double> bias_list, std::vector<double> height_list);
template void STM_calc (Kpoint<std::complex<double>> **Kptr, double *rho, std::vector<double> bias_list, std::vector<double> height_list);
template <typename OrbitalType> void STM_calc (Kpoint<OrbitalType> **Kptr, double *rho, std::vector<double> bias_list, std::vector<double>
height_list)
{

    double z_max = -DBL_MAX;
    double z_min = DBL_MAX;
    for(auto atom = Atoms.begin(); atom!= Atoms.end(); ++atom)
    {
        z_max = std::max(z_max, atom->crds[2] * a0_A);
        z_min = std::min(z_min, atom->crds[2] * a0_A);
    }

    double vaccum = Rmg_L.get_zside() * a0_A - (z_max - z_min);
    for (int i = 0; i < int(height_list.size()); i++)
    {
        if (height_list[i] > vaccum/2.0) 
        {
            height_list[i] = vaccum/2.0;
            printf("STM height list is out of range(Angstrom), max_z =%f, vaccum = %f  height = %f\n", z_max, vaccum, height_list[i]);
        }
    }

    // for STM we don't include augmented charge even when ultra-soft psuedopotentials are used.
    int grid =  1;
    int NX = Rmg_G->get_NX_GRID(grid);
    int NY = Rmg_G->get_NY_GRID(grid);
    int NZ = Rmg_G->get_NZ_GRID(grid);

    int PX0 = Rmg_G->get_PX0_GRID(grid);
    int PY0 = Rmg_G->get_PY0_GRID(grid);
    int PZ0 = Rmg_G->get_PZ0_GRID(grid);

    int PX_OFFSET = Rmg_G->get_PX_OFFSET(grid);
    int PY_OFFSET = Rmg_G->get_PY_OFFSET(grid);
    int PZ_OFFSET = Rmg_G->get_PZ_OFFSET(grid);

    double hz = Rmg_L.get_zside() *a0_A/ NZ;
    int pbasis =  Rmg_G->get_P0_BASIS(grid);
    int nstates = ct.num_states;
    double *work = new double[pbasis ];
    double *work_glob = new double[NX*NY*NZ];

    int iz_max = int((z_max + (height_list.back()) )/hz);
    int iz_idx = iz_max - PZ_OFFSET;

    mkdir("STM", S_IRWXU);
    std::vector<double> rho_xy;
    rho_xy.resize(NX*NY);
    ct.efermi = Fill (Kptr, ct.occ_width, ct.nel, ct.occ_mix, ct.num_states, ct.occ_flag, ct.mp_order);
    for(auto bias_ptr = bias_list.begin(); bias_ptr != bias_list.end(); ++bias_ptr)
    {
        double bias = *bias_ptr;


        std::ostringstream streamObj;
        // Set Fixed -Point Notation
        streamObj << std::fixed;
        streamObj << std::setprecision(2);
        //Add double to stream
        streamObj << std::abs(bias);

        std::string bias_string = streamObj.str();
        if(bias < 0.0) bias_string = "m" + bias_string; 
        std::string filename = "STM/STM_bias_" + bias_string + "_spin"+std::to_string(pct.spinpe) ;
        // change the occupation then calculate the charge density
        double Emin = std::min(ct.efermi*Ha_eV, ct.efermi*Ha_eV + bias);
        double Emax = std::max(ct.efermi*Ha_eV, ct.efermi*Ha_eV + bias);
        for(int kpt = 0;kpt < ct.num_kpts_pe;kpt++)
        {
            for(int st = 0; st < ct.num_states; st++)
            {

                double eig = Kptr[kpt]->Kstates[st].eig[0]*Ha_eV  ;
                if( eig < Emin)
                {
                    Kptr[kpt]->Kstates[st].occupation[0] = 
                        std::exp( - (eig-Emin) * (eig-Emin)/ct.gaus_broad/ct.gaus_broad);
                }
                else if(eig > Emax)
                {
                    Kptr[kpt]->Kstates[st].occupation[0] = 
                        std::exp( - (eig-Emax) * (eig-Emin)/ct.gaus_broad/ct.gaus_broad);
                }
                else
                {
                    Kptr[kpt]->Kstates[st].occupation[0] = 1.0;
                }

                //if(pct.gridpe == 0)printf("\n occ %d %f %f %f", st, ct.efermi * Ha_eV, Kptr[kpt]->Kstates[st].eig[0] * Ha_eV, Kptr[kpt]->Kstates[st].occupation[0]);
            }
        }




        for(int idx = 0;idx < pbasis;idx++)
            work[idx] = 0.0;

        for (int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
        {
            for (int istate = 0; istate < nstates; istate++)
            {

                double scale = Kptr[kpt]->Kstates[istate].occupation[0] * Kptr[kpt]->kp.kweight;

                OrbitalType *psi = Kptr[kpt]->Kstates[istate].psi;

                for (int idx = 0; idx < pbasis; idx++)
                {
                    work[idx] += scale * std::norm(psi[idx]);
                    if(ct.noncoll)
                    {
                        work[idx] += scale * std::norm(psi[idx + pbasis]);
                    }
                }                   /* end for */

            }                       /*end for istate */
        }

        for(int idx =0; idx < NX  * NY * NZ; idx++){
            work_glob[idx] = 0.0;
        }

        MPI_Allreduce(MPI_IN_PLACE, (double *)work, pbasis, MPI_DOUBLE, MPI_SUM, pct.kpsub_comm);

        DistributeToGlobal(work, work_glob);


        OutputCubeFile(work_glob, NX, NY, NZ, filename +".cube", true);


        fill(rho_xy.begin(), rho_xy.end(), 0.0);
        if(iz_idx >=0 && iz_idx < PZ0)
        {
            for(int ix = 0; ix < PX0; ix++)
            {
                for(int iy =0; iy < PY0; iy++)
                {
                    int idx = ix * PY0 * PZ0 + iy*PZ0 + iz_idx;
                    rho_xy[(ix+PX_OFFSET) * NY + iy+PY_OFFSET] = work[idx];
                }
            }
        }

        GlobalSums (rho_xy.data(), NX*NY, pct.grid_comm);
        OutputSTM(rho_xy, NX, NY, filename + ".stm");
    }

}

void OutputSTM(std::vector<double> rho_2d, int NX, int NY, std::string filename)
{

    FILE *fhand = fopen(filename.c_str(), "w");

    fprintf(fhand, "#This is a 2D file to be viewed \n");

    fprintf(fhand, "%d %12.6f %12.6f \n", NX, 
            Rmg_L.get_a0(0) /NX * a0_A, 
            Rmg_L.get_a0(1) /NX * a0_A);
    fprintf(fhand, "%d %12.6f %12.6f \n", NY,
            Rmg_L.get_a1(0) /NY * a0_A, 
            Rmg_L.get_a1(1) /NY * a0_A);

    for (int i=0; i<NX; i++) {
        for (int j=0; j<NY; j++) {
            fprintf(fhand, " %g ", rho_2d[i * NY +j]);
            if( (j % 6) == 5) fprintf(fhand, "\n");
        }

        fprintf(fhand, "\n");
    }
    fclose(fhand);
}

