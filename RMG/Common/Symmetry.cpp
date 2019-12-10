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

#include <cstdint>
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "RmgException.h"
#include "GlobalSums.h"
#include "Symmetry.h"

template void Symmetry::symmetrize_grid_object<double>(double *);
template void Symmetry::symmetrize_grid_object<std::complex<double>>(std::complex<double> *);


void inline symm_ijk(int *srotate, int *strans, int &ix, int &iy, int &iz, int &ixx, int &iyy, int &izz,
        int &nx, int &ny, int &nz)
{

    // rotate ijk by symmetry, and then translation
    int ipoint[3], opoint[3];

    ipoint[0] = ix;
    ipoint[1] = iy;
    ipoint[2] = iz;

    for(int i = 0; i < 3; i++)
    {
        opoint[i] = 0;
        for(int j = 0; j < 3; j++)
            opoint[i] += srotate[i *3 + j] * ipoint[j] ;
        opoint[i] += strans[i];
    }

    ixx = (opoint[0] + nx)%nx;
    while(ixx < 0) ixx += nx;
    while(ixx >= nx) ixx -= nx;
    iyy = (opoint[1] + ny)%ny;
    while(iyy < 0) iyy += ny;
    while(iyy >= ny) iyy -= ny;
    izz = (opoint[2] + nz)%nz;
    while(izz < 0) izz += nz;
    while(izz >= nz) izz -= nz;
}


Symmetry::Symmetry (
          BaseGrid &G_in,
          Lattice &L_in,
          int density) : G(G_in), L(L_in)
{
    symprec = 1.0e-5, angprec = 1.0;
    px_grid = G.get_PX0_GRID(density);
    py_grid = G.get_PY0_GRID(density);
    pz_grid = G.get_PZ0_GRID(density);
    max_pdim = std::max(px_grid, py_grid);
    max_pdim = std::max(max_pdim, pz_grid);

    nx_grid = G.get_NX_GRID(density);
    ny_grid = G.get_NY_GRID(density);
    nz_grid = G.get_NZ_GRID(density);

    xoff = G.get_PX_OFFSET(density);
    yoff = G.get_PY_OFFSET(density);
    zoff = G.get_PZ_OFFSET(density);

    pbasis = px_grid * py_grid * pz_grid;
    nbasis = nx_grid * ny_grid * nz_grid;

    double lattice[9];
    lattice[0*3+0] = L.get_a0(0);
    lattice[1*3+0] = L.get_a0(1);
    lattice[2*3+0] = L.get_a0(2);
    lattice[0*3+1] = L.get_a1(0);
    lattice[1*3+1] = L.get_a1(1);
    lattice[2*3+1] = L.get_a1(2);
    lattice[0*3+2] = L.get_a2(0);
    lattice[1*3+2] = L.get_a2(1);
    lattice[2*3+2] = L.get_a2(2);

    double *tau = new double[4*3 * ct.num_ions];
    int *ityp = new int[4 * ct.num_ions];

    /* Set up atomic positions and species for external routines */
    for (int ion = 0; ion < ct.num_ions; ion++)
    {
        L.to_crystal (&tau[ion * 3], Atoms[ion].crds);
        ityp[ion] = Atoms[ion].species;
    }

    int nsym_atom = spg_get_multiplicity(lattice, tau, ityp, ct.num_ions, symprec, angprec);
    int *sa = new int[9 * nsym_atom];
    std::vector<double> translation(3 * nsym_atom);
    s.resize(9 * nsym_atom);
    ftau.resize(3 * nsym_atom);

    ct.sym_trans = new double[3 * nsym_atom];
    ct.sym_rotate = new int[9 * nsym_atom];

    nsym_atom = spg_get_symmetry(sa, translation.data(),  nsym_atom, lattice, tau, ityp, ct.num_ions, symprec, angprec);

    nsym = 0;
    for(int kpt = 0; kpt < nsym_atom; kpt++)
    {
        double intpart;
        double frac1 = modf(translation[kpt*3 + 0] * nx_grid, &intpart);
        double frac2 = modf(translation[kpt*3 + 1] * ny_grid, &intpart);
        double frac3 = modf(translation[kpt*3 + 2] * nz_grid, &intpart);
        if(frac1 > 0.5) frac1 = 1.0-frac1;
        if(frac2 > 0.5) frac2 = 1.0-frac2;
        if(frac3 > 0.5) frac3 = 1.0-frac3;
        if(frac1 < symprec && frac2 < symprec &&frac3 < symprec)
        {

            if(ct.verbose) printf("\n sym operation after considering real space grid # %d",nsym);
            for(int i = 0; i < 3; i++)
                for(int j = 0; j < 3; j++)
                {
                    s[nsym * 9 + i *3 + j] = sa[kpt * 9 + i *3 + j];
                    ct.sym_rotate[nsym * 9 + i *3 + j] = sa[kpt * 9 + i *3 + j];
                }

            ftau[nsym*3 + 0] = translation[kpt*3 + 0] * nx_grid + symprec;
            ftau[nsym*3 + 1] = translation[kpt*3 + 1] * ny_grid + symprec;
            ftau[nsym*3 + 2] = translation[kpt*3 + 2] * nz_grid + symprec;

            ct.sym_trans[nsym*3+0] = translation[kpt*3+0];
            ct.sym_trans[nsym*3+1] = translation[kpt*3+1];
            ct.sym_trans[nsym*3+2] = translation[kpt*3+2];

            if(ct.verbose)
            {
                for(int i = 0; i < 3; i++)
                {
                    printf("\n      %3d  %3d  %3d", s[nsym * 9 + i *3 + 0],s[nsym * 9 + i *3 + 1],s[nsym * 9 + i *3 + 2]);
                }
                printf("  with translation of (%d %d %d) grids ", ftau[nsym*3 + 0],ftau[nsym*3 + 1],ftau[nsym*3 + 2]);
            }
            nsym++;
        }
        else if(ct.verbose)
        {
            printf("\n translation break a symmetry") ;
            for(int i = 0; i < 3; i++)
            {
                printf("\n      %3d  %3d  %3d", sa[kpt * 9 + i *3 + 0],sa[kpt * 9 + i *3 + 1],sa[kpt * 9 + i *3 + 2]);
            }   
            printf("  with translation of (%f %f %f) grids ", frac1, frac2, frac3);
            
        }
    }   

    // sym index arrays dimensioned to size of smallest possible integer type
    if(max_pdim < 256)
    {
        sym_index_x8.resize(nsym * pbasis);
        sym_index_y8.resize(nsym * pbasis);
        sym_index_z8.resize(nsym * pbasis);
        init_symm_ijk(sym_index_x8, sym_index_y8, sym_index_z8);
    }
    else
    {
        sym_index_x16.resize(nsym * pbasis);
        sym_index_y16.resize(nsym * pbasis);
        sym_index_z16.resize(nsym * pbasis);
        init_symm_ijk(sym_index_x16, sym_index_y16, sym_index_z16);
    }

    delete [] sa;
    delete [] tau;
    delete [] ityp;
}

template <typename T> void Symmetry::init_symm_ijk(std::vector<T> &sym_x_idx, std::vector<T> &sym_y_idx, std::vector<T> &sym_z_idx)
{
    int incx = py_grid * pz_grid;
    int incy = pz_grid;

    int ixx, iyy, izz;
    for(int isy = 0; isy < nsym; isy++)
    {
        for (int ix = 0; ix < px_grid; ix++) {
            for (int iy = 0; iy < py_grid; iy++) {
                for (int iz = 0; iz < pz_grid; iz++) {
                    int ix1 = ix + xoff;
                    int iy1 = iy + yoff;
                    int iz1 = iz + zoff;

                    symm_ijk(&s[isy *9], &ftau[isy*3], ix1, iy1, iz1, ixx, iyy, izz, nx_grid, ny_grid, nz_grid);
                    sym_x_idx[isy * pbasis + ix * incx + iy * incy + iz] = ixx;
                    sym_y_idx[isy * pbasis + ix * incx + iy * incy + iz] = iyy;
                    sym_z_idx[isy * pbasis + ix * incx + iy * incy + iz] = izz;
                }
            }
        }
    }
}

template <typename T>
void Symmetry::symmetrize_grid_object(T *object)
{
    if(max_pdim < 256)
    {
        symmetrize_grid_object_int(object, sym_index_x8, sym_index_y8, sym_index_z8);
    }
    else
    {
        symmetrize_grid_object_int(object, sym_index_x16, sym_index_y16, sym_index_z16);
    }
}

template <typename T, typename U>
void Symmetry::symmetrize_grid_object_int(T *object, const std::vector<U> &sym_x_idx, const std::vector<U> &sym_y_idx, const std::vector<U> &sym_z_idx)
{
    int incx = py_grid * pz_grid;
    int incy = pz_grid;

    int incx1 = 1;
    int incy1 = nx_grid;
    int incz1 = nx_grid * ny_grid;

    // Allocate a global array object and put this processors object into the correct location
    T *da = new T[nbasis]();

    for (int ix = 0; ix < px_grid; ix++) {
        for (int iy = 0; iy < py_grid; iy++) {
            for (int iz = 0; iz < pz_grid; iz++) {
                da[(iz + zoff)*incz1 + (iy + yoff)*incy1 + (ix + xoff)*incx1] = object[ix * incx + iy*incy + iz];
            }
        }
    }

    /* Call global sums to give everyone the full array */
    int length = nbasis;
    if(typeid(T) == typeid(std::complex<double>)) length *= 2;
    GlobalSums ((double *)da, length, pct.grid_comm);

    for(int ix=0;ix < pbasis;ix++) object[ix] = 0.0;

    for(int isy = 0; isy < nsym; isy++)
    {
        for (int ix = 0; ix < px_grid; ix++) {
            for (int iy = 0; iy < py_grid; iy++) {
                for (int iz = 0; iz < pz_grid; iz++) {

                    int ixx = sym_x_idx[isy * pbasis + ix * incx + iy * incy + iz] ;
                    int iyy = sym_y_idx[isy * pbasis + ix * incx + iy * incy + iz] ;
                    int izz = sym_z_idx[isy * pbasis + ix * incx + iy * incy + iz] ;

                    object[ix * incx + iy*incy + iz] += da[izz *incz1 + iyy *incy1 + ixx *incx1];
                }
            }
        }
    }

    double t1 = (double) nsym;
    t1 = 1.0 / t1;
    for(int ix = 0; ix < pbasis; ix++) object[ix] = object[ix] * t1;

    delete [] da;

}


Symmetry::~Symmetry(void)
{
    delete [] ct.sym_trans;
    delete [] ct.sym_rotate;
}
