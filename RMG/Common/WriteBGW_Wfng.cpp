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

//  write wave function in gspace for one kpoitn. filename to wfng.complex_kpt(x)
//  the file will be read in fortran with unformatted in Berkeley GW
//  written by Wenchang Lu, NCSU, 2016-10-12

template  void WriteBGW_Wfng (int kpt, Kpoint<double> * kptr);
template  void WriteBGW_Wfng (int kpt, Kpoint<std::complex<double> > * kptr);
    template <typename KpointType>
void WriteBGW_Wfng (int kpt, Kpoint<KpointType> * kptr)
{
    SPECIES *sp;
    ION *iptr;
    int amode, fhand=0;
    char stitle[32], sdate[32], stime[32];
    double wfng_dk1, wfng_dk2, wfng_dk3;
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
    strncpy(stitle, "WFN-Complex", 11);

    wfng_dk1 = 0.0;
    wfng_dk2 = 0.0;
    wfng_dk3 = 0.0;
    int cell_symmetry = 0; //for cubic and orthorhobic 
    int nrecord = 1;
    int nspin = ct.spin_flag + 1;
    double ecutrho, ecutwfc, at[9], adot[9], bg[9], bdot[9]; 
    double tpiba, recvol;
    ecutrho = ct.ecutrho;
    ecutwfc = ct.ecutwfc;

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


    int NX_GRID = Rmg_G->get_NX_GRID(1);
    int NY_GRID = Rmg_G->get_NY_GRID(1);
    int NZ_GRID = Rmg_G->get_NZ_GRID(1);

    int NBASIS = NX_GRID * NY_GRID * NZ_GRID;
    int *gk_g = new int[3* NBASIS];

    int num_pw_wfc_max;
    int num_pw_wfc;

    num_pw_wfc = 0;
    for(int ix = 0; ix < NX_GRID; ix++)
        for(int iy = 0; iy < NY_GRID; iy++)
            for(int iz = 0; iz < NZ_GRID; iz++)
            {
                ivec[0] = ix;
                ivec[1] = iy;
                ivec[2] = iz;
                if(ivec[0] > NX_GRID/2) ivec[0] = ivec[0] - NX_GRID;
                if(ivec[1] > NY_GRID/2) ivec[1] = ivec[1] - NY_GRID;
                if(ivec[2] > NZ_GRID/2) ivec[2] = ivec[2] - NZ_GRID;

                gvec[0] = (double)ivec[0] * Rmg_L.b0[0] + (double)ivec[1] * Rmg_L.b1[0] + (double)ivec[2] * Rmg_L.b2[0];
                gvec[0] *= Rmg_L.celldm[0];
                gvec[1] = (double)ivec[0] * Rmg_L.b0[1] + (double)ivec[1] * Rmg_L.b1[1] + (double)ivec[2] * Rmg_L.b2[1];
                gvec[1] *= Rmg_L.celldm[0];
                gvec[2] = (double)ivec[0] * Rmg_L.b0[2] + (double)ivec[1] * Rmg_L.b1[2] + (double)ivec[2] * Rmg_L.b2[2];
                gvec[2] *= Rmg_L.celldm[0];

                gvec[0] = gvec[0] * tpiba + ct.kp[kpt].kvec[0];
                gvec[1] = gvec[1] * tpiba + ct.kp[kpt].kvec[1];
                gvec[2] = gvec[2] * tpiba + ct.kp[kpt].kvec[2];

                gmags = gvec[0] * gvec[0] + gvec[1] * gvec[1] + gvec[2] * gvec[2];
                if(gmags  < ecutwfc) 
                {
                    gk_g[num_pw_wfc * 3 + 0] = ivec[0];
                    gk_g[num_pw_wfc * 3 + 1] = ivec[1];
                    gk_g[num_pw_wfc * 3 + 2] = ivec[2];
                    num_pw_wfc++;
                }
            }



    num_pw_wfc_max=num_pw_wfc;

    int nks = ct.num_states;
    int ntot_states = nspin * nks;
    double *occs = new double[ntot_states];
    double *eigs= new double[ntot_states];
    int *ifmin = new int[ntot_states];
    int *ifmax = new int[ntot_states];

    for(int idx = 0; idx < ntot_states; idx++) 
    {
        eigs[idx] = 0.0;
        occs[idx] = 0.0;
    }

    // Fill eigs
    for(int ispin = 0; ispin<nspin; ispin++)
        for(int st = 0;st < ct.num_states;st++) {
            eigs[ispin * nks + st] = kptr->Kstates[st].eig[ispin];
            occs[ispin * nks + st] = kptr->Kstates[st].occupation[ispin];
        }

    MPI_Allreduce(MPI_IN_PLACE, eigs, ntot_states, MPI_DOUBLE, MPI_SUM, pct.spin_comm);
    MPI_Allreduce(MPI_IN_PLACE, occs, ntot_states, MPI_DOUBLE, MPI_SUM, pct.spin_comm);

    if (nspin == 1) for (int i = 0; i < ntot_states; i++) occs[i] *= 0.5;


    //  change the eigvalue unit from Hatree to Rydeberg
    for (int i = 0; i < ntot_states; i++) eigs[i] *= 2.0;

    for (int i = 0; i < nspin; i++) ifmin[i] = 0;
    for (int i = 0; i <  nspin; i++) ifmax[i] = 0;
    for (int i = 0; i <  nspin; i++) 
        for (int st = 0; st < ct.num_states; st++)
        {
            if(occs[i * ct.num_states + st] > 0.5) 
            {
                if(ifmin[i] == 0) ifmin[i] = st +1;
                ifmax[i] = st+1;
            }
        }

    int length;
    int ione = 1;

    /*Only one processor will write restart file*/
    if (pct.gridpe == 0)
    {

        amode = S_IREAD | S_IWRITE;
        std::string filename("wfng.complex_kpt");
        filename = filename + std::to_string(kpt);
        fhand = open((char *)filename.c_str(), O_CREAT | O_TRUNC | O_RDWR, amode);
        length = 96;
        write(fhand, &length, sizeof(int));
        write(fhand, stitle, sizeof(char) *32);
        write(fhand, sdate, sizeof(char) *32);
        write(fhand, stime, sizeof(char) *32);

        write(fhand, &length, sizeof(int));
        length = 8 * sizeof(int) + 2 * sizeof(double);
        write(fhand, &length, sizeof(int));

        write(fhand, &nspin, sizeof(int));
        write(fhand, &num_pw_rho, sizeof(int)); // num_ of pw for rho
        write(fhand, &ct.nsym, sizeof(int));  // num of symmetry
        write(fhand, &cell_symmetry, sizeof(int));
        write(fhand, &ct.num_ions, sizeof(int));
        write(fhand, &ecutrho, sizeof(double));
        write(fhand, &ione, sizeof(int)); // num of kpoint 
        write(fhand, &ct.num_states, sizeof(int)); // number of bands
        write(fhand, &num_pw_wfc_max, sizeof(int)); // max num of plane waves for wavefunction of all kpoint
        write(fhand, &ecutwfc, sizeof(double)); 

        write(fhand, &length, sizeof(int));
        length = 6 * sizeof(int) + 3 * sizeof(double);
        write(fhand, &length, sizeof(int));

        write(fhand, &FNX_GRID, sizeof(int)); 
        write(fhand, &FNY_GRID, sizeof(int)); 
        write(fhand, &FNZ_GRID, sizeof(int)); 
        write(fhand, ct.kpoint_mesh, sizeof(int) *3);

        write(fhand, &wfng_dk1, sizeof(double)); 
        write(fhand, &wfng_dk2, sizeof(double)); 
        write(fhand, &wfng_dk3, sizeof(double)); 

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

        write(fhand, &num_pw_wfc, sizeof(int));
        write(fhand, &length, sizeof(int));

        length = sizeof(double);
        write(fhand, &length, sizeof(int));
        write(fhand, &ct.kp[kpt].kweight, sizeof(double));
        write(fhand, &length, sizeof(int));

        length = 3* sizeof(double);
        write(fhand, &length, sizeof(int));
        write(fhand, ct.kp[kpt].kpt, sizeof(double)*3);
        write(fhand, &length, sizeof(int));

        length =nspin* sizeof(int);
        write(fhand, &length, sizeof(int));
        write(fhand, ifmin, sizeof(int) * nspin);
        write(fhand, &length, sizeof(int));

        length = nspin* sizeof(int);
        write(fhand, &length, sizeof(int));
        write(fhand, ifmax, sizeof(int) * nspin);
        write(fhand, &length, sizeof(int));

        length = nspin *ct.num_states * sizeof(double);
        write(fhand, &length, sizeof(int));
        write(fhand, eigs, sizeof(double) * ct.num_states * nspin); 
        write(fhand, &length, sizeof(int));
        // eigs(ispin * nband * ct.num_kpts + kpt * nband + iband)

        write(fhand, &length, sizeof(int));
        write(fhand, occs, sizeof(double) * ct.num_states * nspin); 
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
        write(fhand, &num_pw_wfc, sizeof(int));
        write(fhand, &length, sizeof(int));

        length = 3*num_pw_wfc * sizeof(int);
        write(fhand, &length, sizeof(int));
        write(fhand, gk_g, sizeof(int)*num_pw_wfc*3);
        write(fhand, &length, sizeof(int));
    }


    int ix, iy, iz, ii, jj, kk;
    int idx2, idx1, incx, incx1, incy, incy1;

    incx = get_PY0_GRID() * get_PZ0_GRID();
    incy = get_PZ0_GRID();
    incx1 = get_NY_GRID() * get_NZ_GRID();
    incy1 = get_NZ_GRID();

    ii =  get_PX_OFFSET();
    jj =  get_PY_OFFSET();
    kk =  get_PZ_OFFSET();


    std::complex<double> *wfng = new std::complex<double> [num_pw_wfc_max];
    int P0_BASIS = Rmg_G->get_P0_BASIS(1);
    std::complex<double> *wfng_dist = new std::complex<double> [P0_BASIS];
    std::complex<double> *wfng_global = new std::complex<double> [NBASIS];


    //  get wfng in gspace for one kpoint
    for(int istate = 0; istate < ct.num_states; istate++)
    {
        for(int ig=0; ig < num_pw_wfc; ig++) wfng[ig] = 0.0;
        for(int ig=0; ig < NBASIS; ig++) wfng_global[ig] = 0.0;
        if(ct.is_gamma)
        {
            for(int idx = 0; idx < P0_BASIS; idx++)
                wfng_dist[idx] = std::complex<double>(std::real(kptr->Kstates[istate].psi[idx]), 0.0);

        }
        else 
        {
            for(int idx = 0; idx < P0_BASIS; idx++)
                wfng_dist[idx] = std::complex<double>(kptr->Kstates[istate].psi[idx]);
        }

        PfftForward(wfng_dist, wfng_dist, *coarse_pwaves);

        for (ix = 0; ix < get_PX0_GRID(); ix++)
            for (iy = 0; iy < get_PY0_GRID(); iy++)
                for (iz = 0; iz < get_PZ0_GRID(); iz++)
                {
                    idx1 = ix * incx + iy * incy + iz;
                    idx2 = (ix + ii) * incx1 + (iy + jj) * incy1 + iz + kk;
                    wfng_global[idx2] = wfng_dist[idx1];
                }


        for(int ig = 0; ig < num_pw_wfc; ig++)
        {
            ix = gk_g[ig *3 + 0];
            iy = gk_g[ig *3 + 1];
            iz = gk_g[ig *3 + 2];

            if(ix < 0 ) ix += NX_GRID;
            if(iy < 0 ) iy += NY_GRID;
            if(iz < 0 ) iz += NZ_GRID;
            idx2 = ix * incx1 + iy  * incy1 + iz;
            wfng[ig] += wfng_global[idx2];

        }

        MPI_Allreduce(MPI_IN_PLACE, wfng, 2*num_pw_wfc, MPI_DOUBLE, MPI_SUM, pct.grid_comm);

        double norm_const = 0.0;
        for(int ig = 0; ig < num_pw_wfc; ig++)
            norm_const += std::norm(wfng[ig]);

        for(int ig = 0; ig < num_pw_wfc; ig++)
            wfng[ig] /= sqrt(norm_const);


        if(pct.gridpe == 0) 
        {

            length = sizeof(int);
            write(fhand, &length, sizeof(int));
            write(fhand, &nrecord, sizeof(int));
            write(fhand, &length, sizeof(int));

            write(fhand, &length, sizeof(int));
            write(fhand, &num_pw_wfc, sizeof(int));
            write(fhand, &length, sizeof(int));

            length = num_pw_wfc * sizeof(std::complex<double>);
            write(fhand, &length, sizeof(int));
            write(fhand, wfng, sizeof(std::complex<double>) * num_pw_wfc);
            write(fhand, &length, sizeof(int));


        }
    }

    if(pct.gridpe == 0) close(fhand);
}





