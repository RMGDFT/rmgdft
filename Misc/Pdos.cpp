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
#include "Pdos.h"
#include "Scalapack.h"

template Pdos<double>::Pdos(BaseGrid &, Lattice &, const std::string , int, double *, Kpoint<double> **Kptr);
template Pdos<std::complex<double>>::Pdos(BaseGrid &, Lattice &, const std::string , int, std::complex<double>
*, Kpoint<std::complex<double>> **Kptr);

template Pdos<double>::~Pdos(void);
template Pdos<std::complex<double>>::~Pdos(void);
template <class T> Pdos<T>::~Pdos ()
{
    std::cout << "exit pdos" << std::endl;
}

template <class T> Pdos<T>::Pdos(
        BaseGrid &G_in,
        Lattice &L_in,
        const std::string wavefile_in,
        int nstates_in,
        T *psi_in, Kpoint<T> **Kptr) : 
    G(G_in), L(L_in), wavefile(wavefile_in), nstates(nstates_in), psi(psi_in)
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
    pbasis = px0_grid * py0_grid * pz0_grid;
    pbasis_noncoll = pbasis * ct.noncoll_factor;

    if(!ct.norm_conserving_pp )
    {
        std::cout << "Pdos not programmed for ultra soft pseudopotential" << std::endl;
        return;
    }
    ngrid = G.get_NX_GRID(1) * G.get_NY_GRID(1) * G.get_NZ_GRID(1);
    ngrid_noncoll = ngrid * ct.noncoll_factor;

    std::vector<double> occs;
    occs.resize(nstates, 1.0);
    Exx = new Exxbase<T>(G, G, L, wavefile, nstates, occs.data(), psi, EXX_LOCAL_FFT);
    RmgTimer *RT1 = new RmgTimer("8-Pdos: writesingle file");
    Exx->WriteWfsToSingleFile();
    MPI_Barrier(MPI_COMM_WORLD);
    delete RT1;

}

template  void Pdos<double>::ReadRotatePsi(int ikindex, int isym, int isyma, std::string wavefile, double *psi_k);
template  void Pdos<std::complex<double>>::ReadRotatePsi(int ikindex, int isym, int isyma, std::string wavefile, std::complex<double> *psi_k);
template <class T> void Pdos<T>::ReadRotatePsi(int ikindex, int isym, int isyma, std::string wavefile, T *psi_k)
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
                            psi_k_C[st * pbasis_noncoll + ix * py0_grid * pz0_grid + iy * pz0_grid + iz] = - std::conj(dn_rot);
                            psi_k_C[st * pbasis_noncoll + pbasis + ix * py0_grid * pz0_grid + iy * pz0_grid + iz] = std::conj(up_rot);
                        }
                        else
                        {
                            up_rot = (Rmg_Symm->rot_spin_wave[isyma][0][0]) * up + (Rmg_Symm->rot_spin_wave[isyma][0][1]) * dn;
                            dn_rot = (Rmg_Symm->rot_spin_wave[isyma][1][0]) * up + (Rmg_Symm->rot_spin_wave[isyma][1][1]) * dn;
                            psi_k_C[st * pbasis_noncoll + ix * py0_grid * pz0_grid + iy * pz0_grid + iz] = up_rot;
                            psi_k_C[st * pbasis_noncoll + pbasis + ix * py0_grid * pz0_grid + iy * pz0_grid + iz] = dn_rot;
                        }

                    }
                    else if(isym >= 0)
                    {
                        psi_k[st * pbasis + ix * py0_grid * pz0_grid + iy * pz0_grid + iz]
                            = (psi_map[st * ngrid + ixx * ny_grid * nz_grid + iyy * nz_grid + izz]);
                    }
                    else
                    {
                        psi_k[st * pbasis + ix * py0_grid * pz0_grid + iy * pz0_grid + iz]
                            = MyConj(psi_map[st * ngrid + ixx * ny_grid * nz_grid + iyy * nz_grid + izz]);
                    }
                }
            }
        }
    }

    munmap(psi_map, length);
    close(serial_fd);

}


template  void Pdos<double>::Pdos_calc(Kpoint<double> **Kptr, std::vector<double> eigs);
template  void Pdos<std::complex<double>>::Pdos_calc(Kpoint<std::complex<double>> **Kptr, std::vector<double> eigs);
template <class T> void Pdos<T>::Pdos_calc(Kpoint<T> **Kptr, std::vector<double> eigs)
{

    double tol = 1.0e-5;
    int num_q = ct.klist.num_k_all;

    size_t length = nstates * pbasis_noncoll * sizeof(T);
    T *psi_k = (T *)RmgMallocHost(length);

//  total rho projected in x,y,z direction
// magnatic rho projected in x, y, z direction
    int factor = ct.noncoll_factor * ct.noncoll_factor;
    double_4d_array rhop_x, rhop_y, rhop_z;
    rhop_x.resize(boost::extents[factor][num_q][nstates][nx_grid]);
    rhop_y.resize(boost::extents[factor][num_q][nstates][ny_grid]);
    rhop_z.resize(boost::extents[factor][num_q][nstates][nz_grid]);
    std::fill(rhop_x.origin(), rhop_x.origin() + factor * num_q * nstates * nx_grid, 0.0);
    std::fill(rhop_y.origin(), rhop_y.origin() + factor * num_q * nstates * ny_grid, 0.0);
    std::fill(rhop_z.origin(), rhop_z.origin() + factor * num_q * nstates * nz_grid, 0.0);

    double_2d_array rho_temp;
    rho_temp.resize(boost::extents[factor][pbasis]);

    for(int ik = 0; ik < num_q; ik++)
    {
        int ik_irr = ct.klist.k_map_index[ik];
        int isym = ct.klist.k_map_symm[ik];
        int isyma = std::abs(isym) -1;

    //    RmgTimer *RT1 = new RmgTimer("8-Pdos: read and rotate");
        ReadRotatePsi(ik_irr, isym, isyma, wavefile, psi_k);
    //    delete RT1;

        for(int st = 0; st < nstates; st++) {

            for(int idx = 0; idx < pbasis; idx++){
                rho_temp[0][idx] = std::norm(psi_k[st*pbasis_noncoll + idx]);
                if(ct.noncoll)
                {
                    std::complex<double> psiud = 2.0 * psi_k[st*pbasis_noncoll + idx] * std::conj(psi_k[st*pbasis_noncoll + idx + pbasis]);
                    double psiup = std::norm(psi_k[st*pbasis_noncoll + idx]);
                    double psidn = std::norm(psi_k[st*pbasis_noncoll + idx + pbasis]);
                    rho_temp[0][idx] = psiup + psidn;
                    rho_temp[1][idx] = std::real(psiud);
                    rho_temp[2][idx] = std::imag(psiud);
                    rho_temp[3][idx] = psiup - psidn;
                }
                
            }

            for(int is = 0; is < factor; is++){

                for(int iy = 0; iy < py0_grid; iy++){
                    for(int iz = 0; iz < pz0_grid; iz++){
                        for(int ix = 0; ix < px0_grid; ix++){
                            int idx = ix * py0_grid * pz0_grid + iy * pz0_grid + iz;
                            rhop_x[is][ik][st][ix+px_offset] += rho_temp[is][idx];
                        }
                    }
                }

                for(int iy = 0; iy < py0_grid; iy++){
                    for(int ix = 0; ix < px0_grid; ix++){
                        for(int iz = 0; iz < pz0_grid; iz++){
                            int idx = ix * py0_grid * pz0_grid + iy * pz0_grid + iz;
                            rhop_y[is][ik][st][iy+py_offset] += rho_temp[is][idx];
                        }
                    }
                }

                for(int ix = 0; ix < px0_grid; ix++){
                    for(int iy = 0; iy < py0_grid; iy++){
                        for(int iz = 0; iz < pz0_grid; iz++){
                            int idx = ix * py0_grid * pz0_grid + iy * pz0_grid + iz;
                            rhop_z[is][ik][st][iz+pz_offset] += rho_temp[is][idx];
                        }
                    }
                }
            }
        }
    }

    RmgTimer *RT1 = new RmgTimer("8-Pdos: Reduce");
    int count = num_q * nstates * factor;
    MPI_Allreduce(MPI_IN_PLACE, rhop_x.origin(), count * nx_grid, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
    MPI_Allreduce(MPI_IN_PLACE, rhop_y.origin(), count * ny_grid, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
    MPI_Allreduce(MPI_IN_PLACE, rhop_z.origin(), count * nz_grid, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
    delete RT1;

    double delta_e = (ct.Emax - ct.Emin)/(ct.E_POINTS-1);
    double gaus_broad = ct.gaus_broad;

    double_3d_array pdos_x, pdos_y, pdos_z;
    pdos_x.resize(boost::extents[factor][ct.E_POINTS][nx_grid]);
    pdos_y.resize(boost::extents[factor][ct.E_POINTS][ny_grid]);
    pdos_z.resize(boost::extents[factor][ct.E_POINTS][nz_grid]);
    std::fill(pdos_x.origin(), pdos_x.origin() + factor * ct.E_POINTS * nx_grid, 0.0);
    std::fill(pdos_y.origin(), pdos_y.origin() + factor * ct.E_POINTS * ny_grid, 0.0);
    std::fill(pdos_z.origin(), pdos_z.origin() + factor * ct.E_POINTS * nz_grid, 0.0);

    for(int ie = pct.imgpe; ie < ct.E_POINTS; ie++) {
        double energy = ct.Emin + ie * delta_e;
        for(int ik = 0; ik < num_q; ik++) {
            int ik_irr = ct.klist.k_map_index[ik];
            for(int st = 0; st < nstates; st++){
                double tem = energy - eigs[ik_irr * nstates + st];
                double dos_one = std::exp(-tem * tem/gaus_broad/gaus_broad) / gaus_broad /std::sqrt(PI);

                if(dos_one < tol) continue;
                for(int is = 0; is < factor; is++) {
                    for(int ix = 0; ix < nx_grid; ix++) {
                        pdos_x[is][ie][ix] += rhop_x[is][ik][st][ix] * dos_one;
                    }
                    for(int iy = 0; iy < ny_grid; iy++) {
                        pdos_y[is][ie][iy] += rhop_y[is][ik][st][iy] * dos_one;
                    }
                    for(int iz = 0; iz < nz_grid; iz++) {
                        pdos_z[is][ie][iz] += rhop_z[is][ik][st][iz] * dos_one;
                    }
                }

            }
        }

        for(int is = 0; is < factor; is++) {
            double tem = 0.0;
            for(int ix = 0; ix < nx_grid; ix++) {
                tem += pdos_x[is][ie][ix];
            }
            double tem1 = 0.0;
            for(int iy = 0; iy < ny_grid; iy++) {
                tem1 += pdos_y[is][ie][iy];
            }

            if(pct.imgpe == 0) printf("\n dos %f %e", energy, tem);
            for(int iz = 0; iz < nz_grid; iz++) {
            }
        }
    }
    count = ct.E_POINTS * factor;
    MPI_Allreduce(MPI_IN_PLACE, pdos_x.origin(), count * nx_grid, MPI_DOUBLE, MPI_SUM, pct.img_comm);
    MPI_Allreduce(MPI_IN_PLACE, pdos_y.origin(), count * ny_grid, MPI_DOUBLE, MPI_SUM, pct.img_comm);
    MPI_Allreduce(MPI_IN_PLACE, pdos_z.origin(), count * nz_grid, MPI_DOUBLE, MPI_SUM, pct.img_comm);

    if(pct.imgpe == 0)
    {
        time_t tt;
        char *timeptr;
        time (&tt);
        timeptr = ctime (&tt);
        mkdir ("PDOS", S_IRWXU);

        double hx = L.get_xside()/nx_grid * a0_A;
        double hy = L.get_yside()/ny_grid * a0_A;
        double hz = L.get_zside()/nz_grid * a0_A;
        for(int is = 0; is < factor; is++) {
            std::string filename = "PDOS/pdos_"+ std::to_string(is) + "_a.dat"; 
            FILE *file = fopen (filename.c_str(), "w");
            fprintf(file, "#Created with RMG on %s\n", timeptr);
            fprintf(file, "#pdos0: total, pdos1,pdos2,pdos3: magnetic rhoa\n");
            fprintf(file, "#*_a, *_b, *_c: projected on lattice vectors a,b, c\n");

            fprintf (file, "#    x[Angstrom]      E[eV]      pdos\n\n");
            for (int ie = 0; ie < ct.E_POINTS; ie++)
            {
                double energy = ct.Emin + ie * delta_e;
                for (int ix = 0; ix < nx_grid; ix++)
                {

                    fprintf (file, " %10.6f %10.6f %12.6e\n",
                            ix * hx, energy, pdos_x[is][ie][ix]);
                }
                fprintf (file, "\n");
            }

            fclose (file);

            filename = "PDOS/pdos_"+ std::to_string(is) + "_b.dat"; 
            file = fopen (filename.c_str(), "w");
            fprintf(file, "#Created with RMG on %s", timeptr);
            fprintf(file, "#pdos0: total, pdos1,pdos2,pdos3: magnetic rhoa");
            fprintf(file, "#*_a, *_b, *_c: projected on lattice vectors a,b, c");

            fprintf (file, "#    x[Angstrom]      E[eV]      pdos\n\n");
            for (int ie = 0; ie < ct.E_POINTS; ie++)
            {
                double energy = ct.Emin + ie * delta_e;
                for (int iy = 0; iy < ny_grid; iy++)
                {

                    fprintf (file, " %10.6f %10.6f %12.6e\n",
                            iy * hy, energy, pdos_y[is][ie][iy]);
                }
                fprintf (file, "\n");
            }

            fclose (file);
            filename = "PDOS/pdos_"+ std::to_string(is) + "_c.dat"; 
            file = fopen (filename.c_str(), "w");
            fprintf(file, "#Created with RMG on %s", timeptr);
            fprintf(file, "#pdos0: total, pdos1,pdos2,pdos3: magnetic rhoa");
            fprintf(file, "#*_a, *_b, *_c: projected on lattice vectors a,b, c");

            fprintf (file, "#    x[Angstrom]      E[eV]      pdos\n\n");
            for (int ie = 0; ie < ct.E_POINTS; ie++)
            {
                double energy = ct.Emin + ie * delta_e;
                for (int iz = 0; iz < nz_grid; iz++)
                {

                    fprintf (file, " %10.6f %10.6f %12.6e\n",
                            iz * hz, energy, pdos_z[is][ie][iz]);
                }
                fprintf (file, "\n");
            }

            fclose (file);

        }
    }

   // RmgFreeHost(psi_k);
}

