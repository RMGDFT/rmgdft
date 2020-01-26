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
#include "rmg_error.h"

Symmetry *Rmg_Symm;

template void Symmetry::symmetrize_grid_object<double>(double *);
template void Symmetry::symmetrize_grid_object<std::complex<double>>(std::complex<double> *);


void symm_ijk(int *srotate, int *strans, int &ix, int &iy, int &iz, int &ixx, int &iyy, int &izz,
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

    nx_grid = G.get_NX_GRID(density);
    ny_grid = G.get_NY_GRID(density);
    nz_grid = G.get_NZ_GRID(density);

    max_pdim = std::max(nx_grid, ny_grid);
    max_pdim = std::max(max_pdim, nz_grid);
    int nx_grid_c = G.get_NX_GRID(1);
    int ny_grid_c = G.get_NY_GRID(1);
    int nz_grid_c = G.get_NZ_GRID(1);

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
    ftau_wave.resize(3 * nsym_atom);

    sym_trans.resize(3 * nsym_atom);
    sym_rotate.resize(9 * nsym_atom);

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

            if(ct.verbose && pct.imgpe == 0) printf("\n sym operation after considering real space grid # %d",nsym);
            for(int i = 0; i < 3; i++)
                for(int j = 0; j < 3; j++)
                {
                    s[nsym * 9 + i *3 + j] = sa[kpt * 9 + i *3 + j];
                    sym_rotate[nsym * 9 + i *3 + j] = sa[kpt * 9 + i *3 + j];
                }

            ftau[nsym*3 + 0] = translation[kpt*3 + 0] * nx_grid + symprec;
            ftau[nsym*3 + 1] = translation[kpt*3 + 1] * ny_grid + symprec;
            ftau[nsym*3 + 2] = translation[kpt*3 + 2] * nz_grid + symprec;
            ftau_wave[nsym*3 + 0] = translation[kpt*3 + 0] * nx_grid_c + symprec;
            ftau_wave[nsym*3 + 1] = translation[kpt*3 + 1] * ny_grid_c + symprec;
            ftau_wave[nsym*3 + 2] = translation[kpt*3 + 2] * nz_grid_c + symprec;

            sym_trans[nsym*3+0] = translation[kpt*3+0];
            sym_trans[nsym*3+1] = translation[kpt*3+1];
            sym_trans[nsym*3+2] = translation[kpt*3+2];

            if(ct.verbose && pct.imgpe == 0)
            {
                for(int i = 0; i < 3; i++)
                {
                    printf("\n      %3d  %3d  %3d", s[nsym * 9 + i *3 + 0],s[nsym * 9 + i *3 + 1],s[nsym * 9 + i *3 + 2]);
                }
                printf("  with translation of (%d %d %d) grids ", ftau[nsym*3 + 0],ftau[nsym*3 + 1],ftau[nsym*3 + 2]);
            }
            nsym++;
        }
        else if(ct.verbose && pct.imgpe == 0)
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

    ct.nsym = nsym;
    if(ct.verbose && pct.imgpe == 0) printf("\n number of sym operation before considering real space grid: %d",nsym_atom);
    if(ct.verbose && pct.imgpe == 0) printf("\n number of sym operation  after considering real space grid: %d",nsym);
    assert(nsym >0);
    if(nsym == 1) ct.is_use_symmetry = 0;

    sym_atom.resize(ct.num_ions* nsym);

    //  determine equivenlent ions after symmetry operation
    double xtal[3];
    double ndim[3];
    ndim[0] = (double)nx_grid;
    ndim[1] = (double)ny_grid;
    ndim[2] = (double)nz_grid;
    for(int isym = 0; isym < nsym; isym++)
    {
        for (int ion = 0; ion < ct.num_ions; ion++)
        {

            // xtal is the coordinates of atom operated by symmetry operatio isym.
            for(int i = 0; i < 3; i++)
            {
                xtal[i] = s[isym *9 + i *3 + 0] * Atoms[ion].xtal[0]
                    + s[isym *9 + i *3 + 1] * Atoms[ion].xtal[1]
                    + s[isym *9 + i *3 + 2] * Atoms[ion].xtal[2] +ftau[isym *3 + i]/ndim[i];

                if(xtal[i] + symprec  < 0.0) xtal[i]= xtal[i] + 1.0;
                if(xtal[i] + symprec >= 1.0) xtal[i]= xtal[i] - 1.0;
            }

            bool find_atom = false;
            bool cond_x = false;
            bool cond_y = false;
            bool cond_z = false;
            int ionb;
            for (ionb = 0; ionb < ct.num_ions; ionb++)
            {
                if(ityp[ion] == ityp[ionb])
                {
                    //                    r =  (xtal[0] - Atoms[ionb].xtal[0]) *(xtal[0] - Atoms[ionb].xtal[0])
                    //                        +(xtal[1] - Atoms[ionb].xtal[1]) *(xtal[1] - Atoms[ionb].xtal[1])
                    //                        +(xtal[2] - Atoms[ionb].xtal[2]) *(xtal[2] - Atoms[ionb].xtal[2]);
                    double mod_x = (xtal[0] - Atoms[ionb].xtal[0]) *(xtal[0] - Atoms[ionb].xtal[0]);
                    double mod_y = (xtal[1] - Atoms[ionb].xtal[1]) *(xtal[1] - Atoms[ionb].xtal[1]);
                    double mod_z = (xtal[2] - Atoms[ionb].xtal[2]) *(xtal[2] - Atoms[ionb].xtal[2]);

                    cond_x = fabs(mod_x - (int) mod_x) < symprec*10 || fabs(mod_x - (int)mod_x) > 1.0-symprec*10;
                    cond_y = fabs(mod_y - (int) mod_y) < symprec*10 || fabs(mod_y - (int)mod_y) > 1.0-symprec*10;
                    cond_z = fabs(mod_z - (int) mod_z) < symprec*10 || fabs(mod_z - (int)mod_z) > 1.0-symprec*10;
                    if(cond_x && cond_y && cond_z)
                    {
                        sym_atom[isym * ct.num_ions + ion] = ionb;
                        find_atom = true;
                    }

                }

                if(find_atom) break;

            }
            if(!find_atom)
            {
                printf("\n Equivalent atom not found %d %d %d %d %e %e %e \n", ion, ionb,isym, nsym, xtal[0],xtal[1], xtal[2]);
                rmg_error_handler(__FILE__, __LINE__, "Exiting.\n");
            }
        }
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

void Symmetry::symforce (void)
{


    double *force = new double[3* ct.num_ions];

    for (int ion = 0; ion < ct.num_ions; ion++)
    {
        for (int ir = 0; ir < 3; ir++)
        {
            force[ion *3 + ir] = Atoms[ion].force[ct.fpt[0]][0] * L.b0[ir] +
                Atoms[ion].force[ct.fpt[0]][1] * L.b1[ir] +
                Atoms[ion].force[ct.fpt[0]][2] * L.b2[ir];
        }                       /* end for ir */

        Atoms[ion].force[ct.fpt[0]][0] = 0.0;
        Atoms[ion].force[ct.fpt[0]][1] = 0.0;
        Atoms[ion].force[ct.fpt[0]][2] = 0.0;
    }
    for (int ion = 0; ion < ct.num_ions; ion++)
    {
        for(int isy = 0; isy < nsym; isy++)
        {
            int ion1 = sym_atom[isy * ct.num_ions + ion];
            for(int i = 0; i < 3; i++)
                for(int j = 0; j < 3; j++)
                    Atoms[ion1].force[ct.fpt[0]][i] += s[isy *9 + i* 3 + j] * force[ion*3 + j];
        }

    }
    for (int ion = 0; ion < ct.num_ions; ion++)
    {
        for (int ir = 0; ir < 3; ir++)
        {
            force[ion *3 + ir] = Atoms[ion].force[ct.fpt[0]][0] * L.a0[ir] +
                Atoms[ion].force[ct.fpt[0]][1] * L.a1[ir] +
                Atoms[ion].force[ct.fpt[0]][2] * L.a2[ir];
        }                       /* end for ir */
    }

    for (int ion = 0; ion < ct.num_ions; ion++)
        for(int i = 0; i < 3; i++)
            Atoms[ion].force[ct.fpt[0]][i] = force[3*ion + i] /nsym;

    delete [] force;
}                               /* end symforce */

void Symmetry::symmetrize_tensor(double *mat_tensor)
{
    // symmetrize the stress tensor matrix
    double work[9], b[9], a[9];
    for (int i = 0; i < 3; i++)
    {
        b[0 * 3 + i] = L.b0[i];
        b[1 * 3 + i] = L.b1[i];
        b[2 * 3 + i] = L.b2[i];
        a[0 * 3 + i] = L.a0[i];
        a[1 * 3 + i] = L.a1[i];
        a[2 * 3 + i] = L.a2[i];
    }

    for(int i = 0; i < 9; i++) work[i] = 0.0;

    // transfer to crystal coordinate
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for(int k = 0; k < 3; k++)
                for(int l = 0; l < 3; l++)
                {
                    work[i*3 + j] += mat_tensor[k * 3 +l] *b[i*3 + k] * b[j*3 + l];
                }

    for(int i = 0; i < 9; i++) mat_tensor[i] = 0.0;
    for(int isy = 0; isy < nsym; isy++)
    {

        for(int i = 0; i < 3; i++)
            for(int j = 0; j < 3; j++)
                for(int k = 0; k < 3; k++)
                    for(int l = 0; l < 3; l++)
                    {
                        mat_tensor[i*3+j] += s[isy * 9 + i * 3 + k] * work[k * 3 + l] * s[isy*9 + j*3 +l];
                    }
    }

    for(int i = 0; i < 9; i++) work[i] = 0.0;
    //transfer bact to cartesian 
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for(int k = 0; k < 3; k++)
                for(int l = 0; l < 3; l++)
                {
                    work[i*3 + j] += mat_tensor[k * 3 +l] *a[k*3 + i] * a[l*3 + j];
                }

    for(int i = 0; i < 9; i++) mat_tensor[i] = work[i] / nsym;
}


Symmetry::~Symmetry(void)
{
}
