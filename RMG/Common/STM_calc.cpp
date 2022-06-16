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
    int grid =  Rmg_G->default_FG_RATIO;
    int NX = Rmg_G->get_NX_GRID(grid);
    int NY = Rmg_G->get_NY_GRID(grid);
    int NZ = Rmg_G->get_NZ_GRID(grid);

    int PX0 = Rmg_G->get_PX0_GRID(grid);
    int PY0 = Rmg_G->get_PY0_GRID(grid);
    int PZ0 = Rmg_G->get_PZ0_GRID(grid);

    int FPX_OFFSET, FPY_OFFSET, FPZ_OFFSET;
    Rmg_G->find_node_offsets(pct.gridpe,NX, NY, NZ, &FPX_OFFSET, &FPY_OFFSET, &FPZ_OFFSET);
    double hz = Rmg_L.get_zside() *a0_A/ NZ;

    int iz_max = int((z_max + (height_list.back()) )/hz);
    mkdir("STM", S_IRWXU);
    std::vector<double> rho_xy;
    std::vector<double> rho_3d;
    rho_3d.resize(NX * NY * NZ, 0.0);
    rho_xy.resize(NX*NY);
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

        int factor = ct.noncoll_factor * ct.noncoll_factor;

        int ratio = Rmg_G->default_FG_RATIO;
        int FP0_BASIS = Rmg_G->get_P0_BASIS(ratio);

        GetNewRhoPost(Kptr, rho);

        if(!ct.norm_conserving_pp) {
            double *augrho = new double[FP0_BASIS*factor]();
            GetAugRho(Kptr, augrho);
            for(int idx = 0;idx < FP0_BASIS*factor;idx++) rho[idx] += augrho[idx];
            delete [] augrho;
        }

        if (ct.nspin == 2)
            get_rho_oppo (rho,  &rho[FP0_BASIS]);
        if(ct.AFM)
        {
            Rmg_Symm->symmetrize_rho_AFM(rho, &rho[FP0_BASIS]);
        }
        else
        {
            if(Rmg_Symm) Rmg_Symm->symmetrize_grid_object(rho);
            if(ct.noncoll && Rmg_Symm)
                Rmg_Symm->symmetrize_grid_vector(&rho[FP0_BASIS]);
        }

        OutputCubeFile(rho, Rmg_G->default_FG_RATIO, filename +".cube");

        double rho_max = 0.0;
        double rho_min = DBL_MAX;
        for (int i =0; i < int(height_list.size()); i++)
        {
            std::ostringstream HeightObj;
            // Set Fixed -Point Notation
            HeightObj << std::fixed;
            HeightObj << std::setprecision(1);
            //Add double to Height
            HeightObj << std::abs(height_list[i]);
            fill(rho_xy.begin(), rho_xy.end(), 0.0);
            int iz = int((z_max + height_list[i] )/hz);

            for(int ix = 0; ix < PX0; ix++)
            {
                for(int iy = 0; iy < PY0; iy++)
                {
                    if (iz >=FPZ_OFFSET && iz < PZ0 + FPZ_OFFSET)
                    {
                        rho_xy[(ix+FPX_OFFSET) * NY + iy + FPY_OFFSET] = rho[ix * PY0 * PZ0 + iy * PZ0 + iz - FPZ_OFFSET];
                    }
                }
            }

            int length = NX * NY;
            MPI_Allreduce(MPI_IN_PLACE, rho_xy.data(), length, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
            double rho_0  = *std::max_element(rho_xy.begin(), rho_xy.end());
            double rho_1  = *std::min_element(rho_xy.begin(), rho_xy.end());
            rho_max = std::max(rho_max, rho_0);
            rho_min = std::min(rho_min, rho_1);
            if(pct.imgpe == 0)
                OutputSTM(rho_xy, NX, NY, filename + "_height_"+HeightObj.str() + ".stm");
        }

        // constant current mode STM

        double rho_ave = (rho_max + rho_min) * 0.5;
        if(pct.imgpe == 0 ) printf("\n rho_ave %e at bias %f\n", rho_ave, bias);

        std::fill(rho_3d.begin(), rho_3d.end(), 0.0);
        for(int ix = 0; ix < PX0; ix++)
        {
            for(int iy = 0; iy < PY0; iy++)
            {
                for(int iz = 0; iz < PZ0; iz++)
                {
                    rho_3d[(ix+FPX_OFFSET) * NY *NZ  + (iy + FPY_OFFSET) * NZ + iz +FPZ_OFFSET ] = rho[ix * PY0 * PZ0 + iy * PZ0 + iz];
                }
            }
        }

        int length = NX * NY * NZ;
        MPI_Allreduce(MPI_IN_PLACE, rho_3d.data(), length, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
        std::fill(rho_xy.begin(), rho_xy.end(), 0.0);
        for(int ix = pct.gridpe; ix < NX; ix+= pct.grid_npes)
        {
            for(int iy = 0; iy < NY; iy++)
            {
                int iz0 = -1;
                for(int iz = iz_max; iz >=0; iz--)
                {
                    if( rho_3d[ ix *NY *NZ + iy * NZ + iz] > rho_ave)
                    {
                        iz0 = iz;
                        break;
                    }
                }

                if( iz0 == -1)
                {
                    rho_xy[ix*NY + iy] = 0.0;
                }
                else
                {
                    double  rho1 = rho_3d[ ix * NY * NZ + iy *NZ + iz0];
                    double  rho2 = rho_3d[  ix * NY * NZ + iy *NZ + iz0 +1];

                    if(rho2 > rho1)
                    {
                        printf("\n charge density wrong at ix iy iz %d %d %d   %e %e", ix, iy, iz0, rho1, rho2);
                    }
                    rho_xy[ix*NY + iy] = (iz0 + (rho1 - rho_ave)/(rho1 - rho2) ) * hz;
                    //if(rho_xy[ix*NY + iy] > 30.0) printf("\n WARR: %d %d %d %e %e %e", ix, iy, iz0, rho_ave, rho1, rho2);
                }

            }
        }

        length = NX * NY;
        MPI_Allreduce(MPI_IN_PLACE, rho_xy.data(), length, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
        if(pct.imgpe == 0)
            OutputSTM(rho_xy, NX, NY, filename + "_ConsCurrent.stm");

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
