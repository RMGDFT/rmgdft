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
#include "RmgParallelFft.h"

// write rho in G-space to be read in BGW
// the file is for BGW read in fortran with unformatted 
// Written By Wenchang Lu, NCSU, 2016-10

void WriteBGW_Rhog (double *rho, double *rho_oppo)
{
    int amode;
    ION *iptr;
    SPECIES *sp;
    char stitle[32], sdate[32], stime[32];
    time_t tt;
    time(&tt);
    char *timeptr;
    timeptr = ctime(&tt);
    for(int i = 0; i < 32; i++)
    {
        strncpy(&sdate[i], " ", 1);
        strncpy(&stime[i], " ", 1);
        strncpy(&stitle[i], " ", 1);
    }
    strncpy(sdate, timeptr, 9);
    strncpy(stime, &timeptr[9], 9);
    strncpy(stitle, "RHO-Complex", 11);

    int cell_symmetry = 0; //for cubic and orthorhobic 
    int nrecord = 1;
    int nspin = ct.spin_flag + 1;
    double ecutrho, at[9], adot[9], bg[9], bdot[9]; 
    double tpiba, recvol;
    ecutrho = ct.ecutrho;

    at[0] = Rmg_L.get_a0(0);
    at[1] = Rmg_L.get_a0(1)/at[0];
    at[2] = Rmg_L.get_a0(2)/at[0];
    at[3] = Rmg_L.get_a1(0)/at[0];
    at[4] = Rmg_L.get_a1(1)/at[0];
    at[5] = Rmg_L.get_a1(2)/at[0];
    at[6] = Rmg_L.get_a2(0)/at[0];
    at[7] = Rmg_L.get_a2(1)/at[0];
    at[8] = Rmg_L.get_a2(2)/at[0];

    bg[0] = Rmg_L.get_b0(0)*at[0];
    bg[1] = Rmg_L.get_b0(1)*at[0];
    bg[2] = Rmg_L.get_b0(2)*at[0];
    bg[3] = Rmg_L.get_b1(0)*at[0];
    bg[4] = Rmg_L.get_b1(1)*at[0];
    bg[5] = Rmg_L.get_b1(2)*at[0];
    bg[6] = Rmg_L.get_b2(0)*at[0];
    bg[7] = Rmg_L.get_b2(1)*at[0];
    bg[8] = Rmg_L.get_b2(2)*at[0];

    at[0] = 1.0;
    double alat, alat2, omega, tpiba2;
    alat = Rmg_L.get_a0(0);
    alat2 = alat * alat;
    omega = Rmg_L.get_omega();
    tpiba = 2.0 * PI /alat;
    tpiba2 = tpiba * tpiba;
    recvol = 8.0 * PI*PI*PI / omega;

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            adot[i * 3 + j] = 0.0;
            for (int k = 0; k < 3; k++)
                adot[i*3 + j] += at[i *3 + k] * at[k*3 + j]*alat2;
        }

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            bdot[i * 3 + j] = 0.0;
            for (int k = 0; k < 3; k++)
                bdot[i*3 + j] += bg[i *3 + k] * bg[k*3 + j]*tpiba2;
        }

    int num_pw_rho;
    num_pw_rho = 0;
    int ivec[3];
    double gvec[3], gmags;
    int FNX_GRID = Rmg_G->get_NX_GRID(Rmg_G->default_FG_RATIO);
    int FNY_GRID = Rmg_G->get_NY_GRID(Rmg_G->default_FG_RATIO);
    int FNZ_GRID = Rmg_G->get_NZ_GRID(Rmg_G->default_FG_RATIO);
    int FNBASIS = FNX_GRID * FNY_GRID * FNZ_GRID;

    int *g_g = new int[3*FNX_GRID * FNY_GRID * FNZ_GRID];
    for(int ix = 0; ix < FNX_GRID; ix++)
        for(int iy = 0; iy < FNY_GRID; iy++)
            for(int iz = 0; iz < FNZ_GRID; iz++)
            {
                ivec[0] = ix;
                ivec[1] = iy;
                ivec[2] = iz;
                if(ivec[0] > FNX_GRID/2) ivec[0] = ivec[0] - FNX_GRID;
                if(ivec[1] > FNY_GRID/2) ivec[1] = ivec[1] - FNY_GRID;
                if(ivec[2] > FNZ_GRID/2) ivec[2] = ivec[2] - FNZ_GRID;

                gvec[0] = (double)ivec[0] * Rmg_L.b0[0] + (double)ivec[1] * Rmg_L.b1[0] + (double)ivec[2] * Rmg_L.b2[0];
                gvec[0] *= Rmg_L.celldm[0];
                gvec[1] = (double)ivec[0] * Rmg_L.b0[1] + (double)ivec[1] * Rmg_L.b1[1] + (double)ivec[2] * Rmg_L.b2[1];
                gvec[1] *= Rmg_L.celldm[0];
                gvec[2] = (double)ivec[0] * Rmg_L.b0[2] + (double)ivec[1] * Rmg_L.b1[2] + (double)ivec[2] * Rmg_L.b2[2];
                gvec[2] *= Rmg_L.celldm[0];

                gmags = gvec[0] * gvec[0] + gvec[1] * gvec[1] + gvec[2] * gvec[2];

                if(gmags * tpiba2 < ecutrho) 
                {
                    g_g[num_pw_rho * 3 + 0] = ivec[0];
                    g_g[num_pw_rho * 3 + 1] = ivec[1];
                    g_g[num_pw_rho * 3 + 2] = ivec[2];
                    num_pw_rho++;
                }
            }

    int FP0_BASIS = get_FP0_BASIS();
    int FPX0_GRID = get_FPX0_GRID();
    int FPY0_GRID = get_FPY0_GRID();
    int FPZ0_GRID = get_FPZ0_GRID();

    int incx = FPY0_GRID * FPZ0_GRID;
    int incy = FPZ0_GRID;
    int incx1 = FNY_GRID * FNZ_GRID;
    int incy1 = FNZ_GRID;
    int ix, iy, iz, ii, jj, kk, idx1, idx2;

    ii =  get_FPX_OFFSET();
    jj =  get_FPY_OFFSET();
    kk =  get_FPZ_OFFSET();

    std::complex<double> *rhog = new std::complex<double> [num_pw_rho];
    std::complex<double> *rhog_oppo = new std::complex<double> [num_pw_rho];
    std::complex<double> *rhog_dist = new std::complex<double> [FP0_BASIS];
    std::complex<double> *rhog_global = new std::complex<double> [FNBASIS];


    for(int idx = 0; idx < FP0_BASIS; idx++)
        rhog_dist[idx] = std::complex<double>(rho[idx], 0.0);

    PfftForward(rhog_dist, rhog_dist, *fine_pwaves);


    
    for(int ig = 0; ig < FNBASIS; ig++) rhog_global[ig] = 0.0;
    for (ix = 0; ix < FPX0_GRID; ix++)
        for (iy = 0; iy < FPY0_GRID; iy++)
            for (iz = 0; iz < FPZ0_GRID; iz++)
            {
                idx1 = ix * incx + iy * incy + iz;
                idx2 = (ix + ii) * incx1 + (iy + jj) * incy1 + iz + kk;
                rhog_global[idx2] = rhog_dist[idx1];
            }


    for(int ig = 0; ig < num_pw_rho; ig++) rhog[ig] = 0.0;
    for(int ig = 0; ig < num_pw_rho; ig++)
    {
        ix = g_g[ig *3 + 0];
        iy = g_g[ig *3 + 1];
        iz = g_g[ig *3 + 2];

        if(ix < 0 ) ix += FNX_GRID;
        if(iy < 0 ) iy += FNY_GRID;
        if(iz < 0 ) iz += FNZ_GRID;
        idx2 = ix * incx1 + iy  * incy1 + iz;
        rhog[ig] += rhog_global[idx2];

    }

    MPI_Allreduce(MPI_IN_PLACE, rhog, 2*num_pw_rho, MPI_DOUBLE, MPI_SUM, pct.grid_comm);

    double norm_const = 0.0;
    for(int ig = 0; ig < num_pw_rho; ig++)
        if(g_g[ig*3] == 0 && g_g[ig*3+1] ==0 && g_g[ig*3+2] ==0)
            norm_const = std::real(rhog[ig]);




    if(nspin == 2)
    {
        for(int idx = 0; idx < FP0_BASIS; idx++)
            rhog_dist[idx] = std::complex<double>(rho_oppo[idx], 0.0);

        PfftForward(rhog_dist, rhog_dist, *fine_pwaves);

        for(int ig = 0; ig < FNBASIS; ig++) rhog_global[ig] = 0.0;
        for (ix = 0; ix < FPX0_GRID; ix++)
            for (iy = 0; iy < FPY0_GRID; iy++)
                for (iz = 0; iz < FPZ0_GRID; iz++)
                {
                    idx1 = ix * incx + iy * incy + iz;
                    idx2 = (ix + ii) * incx1 + (iy + jj) * incy1 + iz + kk;
                    rhog_global[idx2] = rhog_dist[idx1];
                }


        for(int ig = 0; ig < num_pw_rho; ig++) rhog[ig] = 0.0;
        for(int ig = 0; ig < num_pw_rho; ig++)
        {
            ix = g_g[ig *3 + 0];
            iy = g_g[ig *3 + 1];
            iz = g_g[ig *3 + 2];

            if(ix < 0 ) ix += FNX_GRID;
            if(iy < 0 ) iy += FNY_GRID;
            if(iz < 0 ) iz += FNZ_GRID;
            idx2 = ix * incx1 + iy  * incy1 + iz;
            rhog_oppo[ig] += rhog_global[idx2];

        }

        MPI_Allreduce(MPI_IN_PLACE, rhog_oppo, 2*num_pw_rho, MPI_DOUBLE, MPI_SUM, pct.grid_comm);

        for(int ig = 0; ig < num_pw_rho; ig++)
            if(g_g[ig*3] == 0 && g_g[ig*3+1] ==0 && g_g[ig*3+2] ==0)
                norm_const += std::real(rhog_oppo[ig]);

        for(int ig = 0; ig < num_pw_rho; ig++)
            rhog_oppo[ig] = rhog_oppo[ig] * ct.nel /norm_const;

    }

    for(int ig = 0; ig < num_pw_rho; ig++)
        rhog[ig] = rhog[ig] * ct.nel /norm_const;

    /*Only one processor will write restart file*/
    if (pct.imgpe == 0)
    {

        amode = S_IREAD | S_IWRITE;
        int fhand = open("rhog.complex", O_CREAT | O_TRUNC | O_RDWR, amode);
        int length;
        length = 96;
        write(fhand, &length, sizeof(int));
        write(fhand, stitle, sizeof(char) *32);
        write(fhand, sdate, sizeof(char) *32);
        write(fhand, stime, sizeof(char) *32);
        write(fhand, &length, sizeof(int));

        length = 5 * sizeof(int) + 1 * sizeof(double);
        write(fhand, &length, sizeof(int));

        write(fhand, &nspin, sizeof(int));
        write(fhand, &num_pw_rho, sizeof(int)); // num_ of pw for rho
        write(fhand, &ct.nsym, sizeof(int));  // num of symmetry
        write(fhand, &cell_symmetry, sizeof(int));
        write(fhand, &ct.num_ions, sizeof(int));
        write(fhand, &ecutrho, sizeof(double));
        write(fhand, &length, sizeof(int));

        length = 3 * sizeof(int);
        write(fhand, &length, sizeof(int));
        write(fhand, &FNX_GRID, sizeof(int)); 
        write(fhand, &FNY_GRID, sizeof(int)); 
        write(fhand, &FNZ_GRID, sizeof(int)); 
        write(fhand, &length, sizeof(int));

        length = 20 * sizeof(double);
        write(fhand, &length, sizeof(int));

        write(fhand, &omega, sizeof(double)); 
        write(fhand, &alat, sizeof(double)); 
        write(fhand, at, sizeof(double) * 9);
        write(fhand, adot, sizeof(double) * 9);

        write(fhand, &length, sizeof(int));
        length = 20 * sizeof(double);
        write(fhand, &length, sizeof(int));

        write(fhand, &recvol, sizeof(double)); 
        write(fhand, &tpiba, sizeof(double));
        write(fhand, bg, sizeof(double) * 9);
        write(fhand, bdot, sizeof(double) * 9);
        write(fhand, &length, sizeof(int));

        length = 9*ct.nsym * sizeof(int);
        write(fhand, &length, sizeof(int));

        write(fhand, ct.sym_rotate, sizeof(int) * 9 * ct.nsym);

        write(fhand, &length, sizeof(int));
        length = 3*ct.nsym * sizeof(double);
        write(fhand, &length, sizeof(int));

        write(fhand, ct.sym_trans, sizeof(double) * 3 * ct.nsym);

        write(fhand, &length, sizeof(int));
        length = ct.num_ions *( 3 * sizeof(double) + sizeof(int));
        write(fhand, &length, sizeof(int));

        for(int ion = 0; ion < ct.num_ions; ion++) 
        {
            iptr = &ct.ions[ion];
            sp = &ct.sp[iptr->species];
            double x[3];
            x[0] = iptr->crds[0]/alat;
            x[1] = iptr->crds[1]/alat;
            x[2] = iptr->crds[2]/alat;
            write(fhand, x, sizeof(double)*3);
            write(fhand, &sp->atomic_number, sizeof(int));
        }

        write(fhand, &length, sizeof(int));

        length = sizeof(int);
        write(fhand, &length, sizeof(int));
        write(fhand, &nrecord, sizeof(int));
        write(fhand, &length, sizeof(int));

        length = sizeof(int);
        write(fhand, &length, sizeof(int));
        write(fhand, &num_pw_rho, sizeof(int));
        write(fhand, &length, sizeof(int));

        length = 3*num_pw_rho*sizeof(int);
        write(fhand, &length, sizeof(int));
        write(fhand, g_g, sizeof(int) * 3 * num_pw_rho);
        write(fhand, &length, sizeof(int));

        length = sizeof(int);
        write(fhand, &length, sizeof(int));
        write(fhand, &nrecord, sizeof(int));
        write(fhand, &length, sizeof(int));
        write(fhand, &length, sizeof(int));
        write(fhand, &num_pw_rho, sizeof(int));
        write(fhand, &length, sizeof(int));

        length = sizeof(std::complex<double>) * nspin * num_pw_rho;
        write(fhand, &length, sizeof(int));
        write(fhand, rhog, sizeof(std::complex<double>) * num_pw_rho);
        if(nspin == 2) 
            write(fhand, rhog_oppo, sizeof(std::complex<double>) * num_pw_rho);
        write(fhand, &length, sizeof(int));
        close(fhand);

    }

}

