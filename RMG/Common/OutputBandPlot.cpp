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
#include <stdlib.h>
#include <stdio.h>
#include "portability.h"
#include "transition.h"
#include "const.h"
#include "State.h"
#include "Kpoint.h"
#include "BaseThread.h"
#include "TradeImages.h"
#include "RmgTimer.h"
#include "RmgThread.h"
#include "RmgException.h"
#include "GlobalSums.h"
#include "rmgthreads.h"
#include "vhartree.h"
#include "packfuncs.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "Kpoint.h"
#include "Subdiag.h"
#include "../Headers/prototypes.h"
#include "Plots.h"

extern "C" int spg_get_ir_reciprocal_mesh(int *grid_address,
                               int map[],
                               const int mesh[3],
                               const int is_shift[3],
                               const int is_time_reversal,
                               const double lattice[3][3],
                               const double *position,
                               const int types[],
                               const int num_atom,
                               const double symprec);
extern "C" void spg_get_tetrahedra_relative_grid_address(int relative_grid_address[24][4][3],
                     const double rec_lattice[3][3]);
extern "C" double spg_get_tetrahedra_integration_weight(const double omega,
                      const double tetrahedra_omegas[24][4],
                      const char function);

static int mat_inverse_matrix_d3(double m[3][3], double a[3][3], const double precision);
static double mat_get_determinant_d3(double a[3][3]);
static int grid_address_to_index(int g[3], int mesh[3]);



template void OutputBandPlot(Kpoint<double> **);
template void OutputBandPlot(Kpoint<std::complex<double> > **);

template <typename KpointType>
void OutputBandPlot(Kpoint<KpointType> ** Kptr)
{

    int num_energy, k, l, q, r, gp;
    double energy, kx, ky, kz;
    double max_eig, min_eig, t1;
    char filename[MAX_PATH], pngfile[MAX_PATH];
    FILE *bs_f, *dos_f; 
    if(pct.gridpe == 0)
    {

        
        snprintf(filename, MAX_PATH, "%s%s", ct.basename, ".bandstructure.xmgr");
        bs_f = fopen (filename, "w");
        if(!bs_f) {
            rmg_printf("Unable to write band plot data.\n");
            return;
        }
    }


    max_eig = -1000.0;
    min_eig = 1000.0;
    for (int is = 0; is < ct.num_states; is++)
    {
        for(int ik = 0; ik < ct.num_kpts_pe; ik++)
        {


            if(pct.gridpe == 0) fprintf (bs_f, "\n %4d  %16.8f ", ik, Kptr[ik]->Kstates[is].eig[0] * Ha_eV);
            max_eig = std::max(max_eig,  Kptr[ik]->Kstates[is].eig[0] * Ha_eV);
            min_eig = std::min(min_eig,  Kptr[ik]->Kstates[is].eig[0] * Ha_eV);


        }

        if(pct.gridpe ==0 ) fprintf (bs_f, "\n &&");
    }

    if(pct.gridpe ==0) fclose (bs_f);

#if PLPLOT_LIBS && 0
    if(pct.gridpe == 0)
    {
        double *x = new double[ct.num_kpts_pe];
        double *y = new double[ct.num_kpts_pe * ct.num_states];
        snprintf(pngfile, MAX_PATH, "%s%s", ct.basename, ".bandstructure.png");
        x[0] = 0.0;
        for(int ik = 1; ik < ct.num_kpts_pe; ik++)
        {
            kx = (Kptr[ik]->kpt[0] - Kptr[ik-1]->kpt[0]);
            ky = (Kptr[ik]->kpt[1] - Kptr[ik-1]->kpt[1]);
            kz = (Kptr[ik]->kpt[2] - Kptr[ik-1]->kpt[2]);
            x[ik] = x[ik-1] + sqrt(kx*kx + ky*ky + kz*kz);
        }

        for (int is = 0; is < ct.num_states; is++)
            for(int ik = 0; ik < ct.num_kpts_pe; ik++)
                y[is * ct.num_kpts_pe + ik] = Kptr[ik]->Kstates[is].eig[0] * Ha_eV;
        MultiLinePlot(pngfile, "", "E(eV)", "Band Strucutre", x, y, ct.num_kpts_pe, ct.num_states);
    }

#endif

#if 0

    modf(max_eig, &t1);
    max_eig = t1;
    modf(min_eig, &t1);
    min_eig = t1 -1;
    num_energy = (max_eig - min_eig)/0.1 +1;

    int *mesh, *is_shift;
    double *tau, *dos, *dos_g;
    int *grid_address, *map, *weights, *ir_weights;
    int is_time_reversal = 1;
    int meshsize, num_kpts;
    int *ir_gp, *gp_ir_index;

    int ion, *ityp;
    double symprec = 1.0e-5;

    double lattice[3][3];
    double t_e[24][4];
    int g_addr[3];
    int i, j, kpt;

    mesh = ct.kpoint_mesh;
    is_shift = ct.kpoint_is_shift;
    lattice[0][0] = get_a0(0);
    lattice[0][1] = get_a0(1);
    lattice[0][2] = get_a0(2);
    lattice[1][0] = get_a1(0);
    lattice[1][1] = get_a1(1);
    lattice[1][2] = get_a1(2);
    lattice[2][0] = get_a2(0);
    lattice[2][1] = get_a2(1);
    lattice[2][2] = get_a2(2);


    meshsize = mesh[0] * mesh[1] * mesh[2];
    dos = new double[num_energy];
    dos_g = new double[num_energy];

    grid_address = new int[meshsize*3];
    map = new int[meshsize];
    weights = new int[meshsize];
    ir_weights = new int[meshsize];
    ir_gp = new int[meshsize];
    gp_ir_index = new int[meshsize];


    tau = new double[ct.num_ions *3];
    ityp = new int[ct.num_ions];


    /* Set up atomic positions and species for fortran routines */
    for (ion = 0; ion < ct.num_ions; ion++)
    {
        to_crystal (&tau[ion*3], ct.ions[ion].crds);
        ityp[ion] = ct.ions[ion].species;
    }





    num_kpts = spg_get_ir_reciprocal_mesh(grid_address, map,mesh, is_shift, 
            is_time_reversal, lattice, tau, ityp, ct.num_ions, symprec);


    for(kpt = 0; kpt < meshsize; kpt++)
        weights[kpt] = 0;
    for(kpt = 0; kpt < meshsize; kpt++)
        weights[map[kpt]]++;

    j = 0;
    for (i = 0; i < meshsize; i++) {
        if (weights[i] != 0) {
            ir_gp[j] = i;
            ir_weights[j] = weights[i];
            gp_ir_index[i] = j;
            j++;
        } else {
            gp_ir_index[i] = gp_ir_index[map[i]];
        }
    }



    int relative_grid_address[24][4][3];
    double rec_lat[3][3], iw;

    mat_inverse_matrix_d3(rec_lat, lattice, 1e-5);

    spg_get_tetrahedra_relative_grid_address(relative_grid_address, rec_lat);

    int num_pe;
    double tem, ebroading = 0.4;
    MPI_Comm_size (pct.grid_comm, &num_pe);

    for( i = 0; i < num_energy; i++) dos[i] = 0.0;
    for( i = 0; i < num_energy; i++) dos_g[i] = 0.0;
    for( i = pct.gridpe; i < num_energy; i += num_pe)
    {
        energy = min_eig + i * 0.1;

        for(j = 0; j < ct.num_kpts_pe; j++)
        {
            for (k = 0; k < ct.num_states; k++) {
                for (l = 0; l < 24; l++) {
                    for (q = 0; q < 4; q++) {
                        for (r = 0; r < 3; r++) {
                            g_addr[r] = grid_address[ir_gp[j] *3 + r] +
                                relative_grid_address[l][q][r];
                        }

                        gp = grid_address_to_index(g_addr, mesh);
                        t_e[l][q] = Kptr[gp_ir_index[gp]]->Kstates[k].eig[0] * Ha_eV;
                    }
                }
                iw = spg_get_tetrahedra_integration_weight(energy, t_e, 'I');
                dos[i] += iw/meshsize * ir_weights[j];
                tem = energy - Kptr[j]->Kstates[k].eig[0] * Ha_eV;
                dos_g[i] += exp(-tem * tem /ebroading/ebroading) /meshsize * ir_weights[j] /ebroading/1.7724;
            }
        }

    }


    GlobalSums (dos, num_energy, pct.grid_comm);
    GlobalSums (dos_g, num_energy, pct.grid_comm);
    if(pct.gridpe == 0)
    {
        snprintf(filename, MAX_PATH, "%s%s", ct.basename, ".dos.xmgr");
        dos_f = fopen (filename, "w");

        for( i = 0; i < num_energy; i++)
        {
            energy = min_eig + i * 0.1;
            fprintf (dos_f, "\n %f  %16.8f  %16.8f", energy, dos[i], dos_g[i]);
        }

        fclose (dos_f);
    }

    delete [] grid_address;
    delete [] map ;
    delete [] weights ;
    delete [] ir_weights ;
    delete [] ir_gp ;
    delete [] gp_ir_index ;


    delete [] tau ;
    delete [] ityp ;
#endif
}


static double mat_get_determinant_d3(double a[3][3])
{
    return a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        + a[0][1] * (a[1][2] * a[2][0] - a[1][0] * a[2][2])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
}

static int mat_inverse_matrix_d3(double m[3][3],
        double a[3][3],
        const double precision)
{
    double det;
    det = mat_get_determinant_d3(a);

    m[0][0] = (a[1][1] * a[2][2] - a[1][2] * a[2][1]) / det;
    m[1][0] = (a[1][2] * a[2][0] - a[1][0] * a[2][2]) / det;
    m[2][0] = (a[1][0] * a[2][1] - a[1][1] * a[2][0]) / det;
    m[0][1] = (a[2][1] * a[0][2] - a[2][2] * a[0][1]) / det;
    m[1][1] = (a[2][2] * a[0][0] - a[2][0] * a[0][2]) / det;
    m[2][1] = (a[2][0] * a[0][1] - a[2][1] * a[0][0]) / det;
    m[0][2] = (a[0][1] * a[1][2] - a[0][2] * a[1][1]) / det;
    m[1][2] = (a[0][2] * a[1][0] - a[0][0] * a[1][2]) / det;
    m[2][2] = (a[0][0] * a[1][1] - a[0][1] * a[1][0]) / det;
    return 1;
}

static int grid_address_to_index(int g[3], int mesh[3])
{
    int i;
    int gm[3];

    for (i = 0; i < 3; i++) {
        gm[i] = g[i] % mesh[i];
        if (gm[i] < 0) {
            gm[i] += mesh[i];
        }
    }
    return (gm[0] + gm[1] * mesh[0] + gm[2] * mesh[0] * mesh[1]);
}

