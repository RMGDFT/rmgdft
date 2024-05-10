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
#include "blas.h"
#include "transition.h"
#include "qk_compat.h"

extern "C" int *s_fortran_ptr;
extern "C" double *ft_fortran_ptr;
extern "C" int *irt_fortran_ptr;


Symmetry *Rmg_Symm;


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


Symmetry::Symmetry ( Lattice &L_in, int NX, int NY, int NZ, int density) : L(L_in)
{
    symprec = 1.0e-4, angprec = 1.0;


    time_reversal = true;
    int nx_grid_c = NX;
    int ny_grid_c = NY;
    int nz_grid_c = NZ;

    nx_grid = NX * density;
    ny_grid = NY * density;
    nz_grid = NZ * density;

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

    double sr[3][3];

    double *tau = new double[4*3 * ct.num_ions];
    int *ityp = new int[4 * ct.num_ions];

    /* Set up atomic positions and species for external routines */
    for (int ion = 0; ion < ct.num_ions; ion++)
    {
        L.to_crystal (&tau[ion * 3], Atoms[ion].crds);
        ityp[ion] = Atoms[ion].species;
    }


    int nsym_atom=1;
    nsym_atom = spg_get_multiplicity(lattice, tau, ityp, ct.num_ions, symprec, angprec);

    int nrot = nsym_atom;

    int *sa = new int[9 * nrot]();
    std::vector<double> translation(3 * nrot);

    nsym_atom = spg_get_symmetry(sa, translation.data(),  nsym_atom, lattice, tau, ityp, ct.num_ions, symprec, angprec);

    nrot = nsym_atom;

    ftau.resize(3 * nrot);
    ftau_wave.resize(3 * nrot);
    inv_type.resize(nrot);
    time_rev.resize(nrot);

    sym_trans.resize(3 * nrot);
    sym_rotate.resize(9 * nrot);


    if(ct.verbose && pct.gridpe==0) rmg_printf("nsym_atom = %d\n",nrot);

    if(!ct.time_reversal) time_reversal = false;

    for(int isym = 0; isym < nrot; isym++)
    {

        for(int i = 0; i < 3; i++)
        {
            while(std::abs(translation[isym*3+i] - 1.0) < symprec ) translation[isym*3+i] -=1.0;
            while(std::abs(translation[isym*3+i] + 1.0) < symprec ) translation[isym*3+i] +=1.0;
        }



        for(int i = 0; i < 3; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                sym_rotate[isym * 9 + i *3 + j] = sa[isym * 9 + i *3 + j];
                sr[i][j] = sa[isym * 9 + i *3 + j];
            }
        }


        ftau[isym*3 + 0] = std::round(translation[isym*3 + 0] * nx_grid);
        ftau[isym*3 + 1] = std::round(translation[isym*3 + 1] * ny_grid);
        ftau[isym*3 + 2] = std::round(translation[isym*3 + 2] * nz_grid);
        ftau_wave[isym*3 + 0] = std::round(translation[isym*3 + 0] * nx_grid_c);
        ftau_wave[isym*3 + 1] = std::round(translation[isym*3 + 1] * ny_grid_c);
        ftau_wave[isym*3 + 2] = std::round(translation[isym*3 + 2] * nz_grid_c);

        sym_trans[isym*3+0] = translation[isym*3+0];
        sym_trans[isym*3+1] = translation[isym*3+1];
        sym_trans[isym*3+2] = translation[isym*3+2];
        inv_type[isym] = false;
        if(type_symm(sr) == 2 || type_symm(sr) == 5 || type_symm(sr) == 6)
            inv_type[isym] = true;

        time_rev[isym] = false;

    }   


    sym_atom.resize(ct.num_ions* nsym_atom);

    //  determine equivenlent ions after symmetry operation
    double xtal[3];
    double ndim[3];
    ndim[0] = (double)nx_grid;
    ndim[1] = (double)ny_grid;
    ndim[2] = (double)nz_grid;
    for(int isym = 0; isym < nsym_atom; isym++)
    {
        for (int ion = 0; ion < ct.num_ions; ion++)
        {

            // xtal is the coordinates of atom operated by symmetry operatio isym.
            for(int i = 0; i < 3; i++)
            {
                xtal[i] = sym_rotate[isym *9 + i *3 + 0] * Atoms[ion].xtal[0]
                    + sym_rotate[isym *9 + i *3 + 1] * Atoms[ion].xtal[1]
                    + sym_rotate[isym *9 + i *3 + 2] * Atoms[ion].xtal[2] + translation[isym *3 + i];

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
                rmg_printf("\n Equivalent atom not found %d %d %d %d %e %e %e \n", ion, ionb,isym, nsym, xtal[0],xtal[1], xtal[2]);
                rmg_error_handler(__FILE__, __LINE__, "Exiting.\n");
            }
        }
    }


    //remove the symmetry operaion which break the symmetry by spin polarization
    if(ct.nspin == 2)
    {

        bool sym_break, mag_same, mag_oppo;
        for(int isym = nsym_atom-1; isym > 0; isym--)
        {
            sym_break=false;
            mag_same = true;
            mag_oppo = true;

            for (int ion1 = 0; ion1 < ct.num_ions; ion1++)
            {
                int ion2 = sym_atom[isym * ct.num_ions + ion1];
                double diff1 = std::abs(Atoms[ion1].init_spin_rho - Atoms[ion2].init_spin_rho);
                double diff2 = std::abs(Atoms[ion1].init_spin_rho + Atoms[ion2].init_spin_rho);
                mag_same = mag_same && (diff1 < symprec);
                mag_oppo = mag_oppo && (diff2 < symprec);

                if(!mag_same && !mag_oppo)
                {
                    sym_break = true;
                    break;
                }
                if(!mag_same && !ct.AFM)
                {
                    sym_break = true;
                    break;
                }

            }

            if(mag_oppo && !mag_same && ct.AFM) time_rev[isym] = true;
            if(sym_break)
            {
                ftau.erase(ftau.begin() + isym *3, ftau.begin() + isym * 3 + 3); 
                ftau_wave.erase(ftau_wave.begin() + isym *3, ftau_wave.begin() + isym * 3 + 3); 
                sym_trans.erase(sym_trans.begin() + isym *3, sym_trans.begin() + isym * 3 + 3); 
                sym_rotate.erase(sym_rotate.begin() + isym *9, sym_rotate.begin() + isym * 9 + 9); 
                sym_atom.erase(sym_atom.begin() + isym * ct.num_ions, sym_atom.begin() + (isym+1) * ct.num_ions);
                inv_type.erase(inv_type.begin() + isym, inv_type.begin() + isym + 1);
                time_rev.erase(time_rev.begin() + isym, time_rev.begin() + isym + 1);
                translation.erase(translation.begin() + isym *3, translation.begin() + isym * 3 + 3); 
            }
        }

        if(ct.verbose && pct.imgpe == 0)
        {
            rmg_printf("\n sym operation after considering spin polarization# %d",(int) sym_rotate.size()/9);
            for(int isym = 0; isym <(int) sym_rotate.size()/9; isym++)
            {
                rmg_printf("\n symmetry operation # %d:", isym);
                for(int i = 0; i < 3; i++)
                {
                    rmg_printf("\n      %3d  %3d  %3d", sym_rotate[isym * 9 + i *3 + 0],sym_rotate[isym * 9 + i *3 + 1],sym_rotate[isym * 9 + i *3 + 2]);
                }
                rmg_printf("  with translation of (%d %d %d) grids %d", ftau[isym*3 + 0],ftau[isym*3 + 1],ftau[isym*3 + 2], int(time_rev[isym]));
            }
        }

    }

    //remove the symmetry operaion which break the symmetry by noncollinear spin
    if(ct.noncoll)
    {


        double vec1[3], vec2[3];
        bool sym_break, mag_same, mag_oppo;
        for(int isym = nsym_atom-1; isym > 0; isym--)
        {
            sym_break=false;
            mag_same = true;
            mag_oppo = true;

            for (int ion1 = 0; ion1 < ct.num_ions; ion1++)
            {
                int ion2 = sym_atom[isym * ct.num_ions + ion1];
                vec1[0] = Atoms[ion1].init_spin_x;
                vec1[1] = Atoms[ion1].init_spin_y;
                vec1[2] = Atoms[ion1].init_spin_z;
                vec2[0] = Atoms[ion2].init_spin_x;
                vec2[1] = Atoms[ion2].init_spin_y;
                vec2[2] = Atoms[ion2].init_spin_z;
                symm_vec(isym, vec1);
                double diff1 = std::abs(vec1[0] - vec2[0]) + std::abs(vec1[1] - vec2[1]) + std::abs(vec1[2] - vec2[2]);
                double diff2 = std::abs(vec1[0] + vec2[0]) + std::abs(vec1[1] + vec2[1]) + std::abs(vec1[2] + vec2[2]);
                mag_same = mag_same && (diff1 < symprec);
                mag_oppo = mag_oppo && (diff2 < symprec);

                if(!mag_same && !mag_oppo)
                {
                    sym_break = true;
                    break;
                }

            }

            if(mag_oppo && !mag_same) time_rev[isym] = true;
            if(sym_break)
            {
                ftau.erase(ftau.begin() + isym *3, ftau.begin() + isym * 3 + 3); 
                ftau_wave.erase(ftau_wave.begin() + isym *3, ftau_wave.begin() + isym * 3 + 3); 
                sym_trans.erase(sym_trans.begin() + isym *3, sym_trans.begin() + isym * 3 + 3); 
                sym_rotate.erase(sym_rotate.begin() + isym *9, sym_rotate.begin() + isym * 9 + 9); 
                sym_atom.erase(sym_atom.begin() + isym * ct.num_ions, sym_atom.begin() + (isym+1) * ct.num_ions);
                inv_type.erase(inv_type.begin() + isym, inv_type.begin() + isym + 1);
                time_rev.erase(time_rev.begin() + isym, time_rev.begin() + isym + 1);
                translation.erase(translation.begin() + isym *3, translation.begin() + isym * 3 + 3); 
            }
        }

        if(ct.verbose && pct.imgpe == 0)
        {
            rmg_printf("\n sym operation after considering noncollinear spin # %d",(int) sym_rotate.size()/9);
            for(int isym = 0; isym <(int) sym_rotate.size()/9; isym++)
            {
                rmg_printf("\n symmetry operation # %d:", isym);
                for(int i = 0; i < 3; i++)
                {
                    rmg_printf("\n      %3d  %3d  %3d", sym_rotate[isym * 9 + i *3 + 0],sym_rotate[isym * 9 + i *3 + 1],sym_rotate[isym * 9 + i *3 + 2]);
                }
                rmg_printf("  with translation of (%d %d %d) grids ", ftau[isym*3 + 0],ftau[isym*3 + 1],ftau[isym*3 + 2]);
            }
        }

    }

    nsym = sym_rotate.size()/9;
    std::vector<int> sym_to_be_removed;
    for (int isym = 0; isym < nsym; isym++)
    {
        bool already_in = false;
        for (int jsym = 0; jsym < isym; jsym++)
        {
            double delta = 0.0;
            for(int i = 0; i < 3; i++)
            {
                for(int j = 0; j < 3; j++)
                {
                    delta += std::abs(sym_rotate[isym * 9 + i *3 + j] - sym_rotate[jsym * 9 + i *3 + j]);
                }
            }

            if (delta < symprec && time_rev[isym] == time_rev[jsym])
            {
                already_in = true;
                break;
            }
        }

        if(already_in)
        {
            sym_to_be_removed.push_back(isym);
        }
    }

    if (sym_to_be_removed.size() > 0)
    {
        for (auto it = sym_to_be_removed.rbegin(); it!= sym_to_be_removed.rend(); ++it)
        {
            int isym = *it;
            ftau.erase(ftau.begin() + isym *3, ftau.begin() + isym * 3 + 3); 
            ftau_wave.erase(ftau_wave.begin() + isym *3, ftau_wave.begin() + isym * 3 + 3); 
            sym_trans.erase(sym_trans.begin() + isym *3, sym_trans.begin() + isym * 3 + 3); 
            sym_rotate.erase(sym_rotate.begin() + isym *9, sym_rotate.begin() + isym * 9 + 9); 
            sym_atom.erase(sym_atom.begin() + isym * ct.num_ions, sym_atom.begin() + (isym+1) * ct.num_ions);
            inv_type.erase(inv_type.begin() + isym, inv_type.begin() + isym + 1);
            time_rev.erase(time_rev.begin() + isym, time_rev.begin() + isym + 1);
            translation.erase(translation.begin() + isym *3, translation.begin() + isym * 3 + 3); 
        }
    }

    if(ct.verbose && pct.imgpe == 0)
    {
        rmg_printf("\n sym operation after removing duplicating rotation symmetries with different translation # %d",(int) sym_rotate.size()/9);
        for(int isym = 0; isym <(int) sym_rotate.size()/9; isym++)
        {
            rmg_printf("\n symmetry operation # %d:", isym);
            for(int i = 0; i < 3; i++)
            {
                rmg_printf("\n      %3d  %3d  %3d", sym_rotate[isym * 9 + i *3 + 0],sym_rotate[isym * 9 + i *3 + 1],sym_rotate[isym * 9 + i *3 + 2]);
            }
            rmg_printf("  with translation of (%d %d %d) grids, time_rev %d", ftau[isym*3 + 0],ftau[isym*3 + 1],ftau[isym*3 + 2], time_rev[isym]);
        }
    }
    nsym = (int)sym_rotate.size()/9;

    // if use_symmetry = false and AFM = true, just use one time_rev symmetry + unitary
    if(!ct.is_use_symmetry)
    {
        sym_to_be_removed.clear();
        if(ct.AFM)
        {
            bool find_AFM_sym = false;
            for(int isy = 1; isy < nsym; isy++)
            {
                if(time_rev[isy] && !find_AFM_sym) 
                {
                    find_AFM_sym = true;
                }
                else
                {
                    sym_to_be_removed.push_back(isy);
                }

            }
        }
        else 
        {
            for(int isy = 1; isy < nsym; isy++)
                sym_to_be_removed.push_back(isy);
        }

        if (sym_to_be_removed.size() > 0)
        {
            for (auto it = sym_to_be_removed.rbegin(); it!= sym_to_be_removed.rend(); ++it)
            {
                int isym = *it;
                ftau.erase(ftau.begin() + isym *3, ftau.begin() + isym * 3 + 3); 
                ftau_wave.erase(ftau_wave.begin() + isym *3, ftau_wave.begin() + isym * 3 + 3); 
                sym_trans.erase(sym_trans.begin() + isym *3, sym_trans.begin() + isym * 3 + 3); 
                sym_rotate.erase(sym_rotate.begin() + isym *9, sym_rotate.begin() + isym * 9 + 9); 
                sym_atom.erase(sym_atom.begin() + isym * ct.num_ions, sym_atom.begin() + (isym+1) * ct.num_ions);
                inv_type.erase(inv_type.begin() + isym, inv_type.begin() + isym + 1);
                time_rev.erase(time_rev.begin() + isym, time_rev.begin() + isym + 1);
                translation.erase(translation.begin() + isym *3, translation.begin() + isym * 3 + 3); 
            }
        }
    }

    nsym_full = (int)sym_rotate.size()/9;
    if(nsym > 48)
    {
        std::cout << " too many rotation symmetries nsym = "<< nsym << std::endl;
        rmg_error_handler(__FILE__, __LINE__, "too many  symmetry .\n");
    }

    full_sym_rotate.resize(9 * nsym_full);
    full_sym_atom.resize(ct.num_ions * nsym_full);
    full_time_rev.resize(nsym_full);

    full_sym_rotate = sym_rotate;
    full_sym_atom = sym_atom;
    full_time_rev = time_rev;

    sym_to_be_removed.clear();
    std::vector<int> sym_keep;
    for(int kpt = 0; kpt < nsym_full; kpt++)
    {

        double intpart;

        double frac_sum = std::abs(translation[kpt*3+0]) + std::abs(translation[kpt*3+1]) + std::abs(translation[kpt*3+2]);
        if(!ct.frac_symm  && frac_sum > symprec) 
        {
            sym_to_be_removed.push_back(kpt);
            continue;
        }
        double frac1 = modf(translation[kpt*3 + 0] * nx_grid, &intpart);
        double frac2 = modf(translation[kpt*3 + 1] * ny_grid, &intpart);
        double frac3 = modf(translation[kpt*3 + 2] * nz_grid, &intpart);
        if(frac1 > 0.5) frac1 = 1.0-frac1;
        if(frac2 > 0.5) frac2 = 1.0-frac2;
        if(frac3 > 0.5) frac3 = 1.0-frac3;
        if(frac1 > symprec || frac2 > symprec ||frac3 > symprec)
        {
            sym_to_be_removed.push_back(kpt);
            continue;
        }
        sym_keep.push_back(kpt);
    }



    for(size_t isym = 0; isym < sym_keep.size(); isym++)
    {
        int jsym = sym_keep[isym];
        for(int i = 0; i < 9; i++)
        {
            full_sym_rotate[isym*9+i] = sym_rotate[jsym*9+i];
        }
        full_time_rev[isym] = time_rev[jsym];

        for(int i = 0; i < ct.num_ions; i++)
        {
            full_sym_atom[isym*ct.num_ions+i] = sym_atom[jsym*ct.num_ions+i];
        }
    }
    for(size_t isym = 0; isym < sym_to_be_removed.size(); isym++)
    {
        int jsym = sym_to_be_removed[isym];
        int iisym = isym + sym_keep.size();
        for(int i = 0; i < 9; i++)
        {
            full_sym_rotate[iisym*9+i] = sym_rotate[jsym*9+i];
        }
        full_time_rev[iisym] = time_rev[jsym];

        for(int i = 0; i < ct.num_ions; i++)
        {
            full_sym_atom[iisym*ct.num_ions+i] = sym_atom[jsym*ct.num_ions+i];
        }
    }

    if (sym_to_be_removed.size() > 0)
    {
        for (auto it = sym_to_be_removed.rbegin(); it!= sym_to_be_removed.rend(); ++it)
        {
            int isym = *it;
            ftau.erase(ftau.begin() + isym *3, ftau.begin() + isym * 3 + 3); 
            ftau_wave.erase(ftau_wave.begin() + isym *3, ftau_wave.begin() + isym * 3 + 3); 
            sym_trans.erase(sym_trans.begin() + isym *3, sym_trans.begin() + isym * 3 + 3); 
            sym_rotate.erase(sym_rotate.begin() + isym *9, sym_rotate.begin() + isym * 9 + 9); 
            sym_atom.erase(sym_atom.begin() + isym * ct.num_ions, sym_atom.begin() + (isym+1) * ct.num_ions);
            inv_type.erase(inv_type.begin() + isym, inv_type.begin() + isym + 1);
            time_rev.erase(time_rev.begin() + isym, time_rev.begin() + isym + 1);
            translation.erase(translation.begin() + isym *3, translation.begin() + isym * 3 + 3); 
        }
    }

    nsym = (int)sym_rotate.size()/9;

    if(ct.verbose && pct.imgpe == 0)
    {
        rmg_printf("\n sym operation after real space grid  # %d",(int) sym_rotate.size()/9);
        for(int isym = 0; isym <(int) sym_rotate.size()/9; isym++)
        {
            rmg_printf("\n symmetry operation # %d:", isym);
            for(int i = 0; i < 3; i++)
            {
                rmg_printf("\n      %3d  %3d  %3d", sym_rotate[isym * 9 + i *3 + 0],sym_rotate[isym * 9 + i *3 + 1],sym_rotate[isym * 9 + i *3 + 2]);
            }
            rmg_printf("  with translation of (%d %d %d) grids, time_rev %d", ftau[isym*3 + 0],ftau[isym*3 + 1],ftau[isym*3 + 2], time_rev[isym]);
        }
    }


    n_time_rev = 0;
    for (auto t_r = time_rev.begin(); t_r != time_rev.end(); ++t_r)
    {
        if(*t_r) n_time_rev++;
    }

    if(n_time_rev == 0 && ct.AFM)
    {
        std::cout << " no AFM symmetry is found from initial spin set up " << std::endl;
        rmg_error_handler(__FILE__, __LINE__, "no AFM symmetry found Exiting.\n");
    }
    rotate_ylm();
    rotate_spin(); 


    delete [] sa;
    delete [] tau;
    delete [] ityp;
}

void Symmetry::init_symm_ijk()
{
    int incx = py_grid * pz_grid;
    int incy = pz_grid;

    int incx1 = ny_grid * nz_grid;
    int incy1 = nz_grid;
    int incz1 = 1;

    int ixx, iyy, izz;
    for(int isy = 0; isy < nsym; isy++)
    {
        for (int ix = 0; ix < px_grid; ix++) {
            for (int iy = 0; iy < py_grid; iy++) {
                for (int iz = 0; iz < pz_grid; iz++) {
                    int ix1 = ix + xoff;
                    int iy1 = iy + yoff;
                    int iz1 = iz + zoff;

                    symm_ijk(&sym_rotate[isy *9], &ftau[isy*3], ix1, iy1, iz1, ixx, iyy, izz, nx_grid, ny_grid, nz_grid);
                    int idx = izz *incz1 + iyy *incy1 + ixx *incx1;
                    sym_idx[isy * pbasis + ix * incx + iy * incy + iz] = idx;
                }
            }
        }
    }
}
void Symmetry::symmetrize_grid_vector(double *object)
{
    int incx = py_grid * pz_grid;
    int incy = pz_grid;

    int incx1 = ny_grid * nz_grid;
    int incy1 = nz_grid;
    int incz1 = 1;

    // Allocate a global array object and put this processors object into the correct location
    double *da = new double[nbasis*3]();

    for(int is = 0; is < 3; is++)
    {
        for (int ix = 0; ix < px_grid; ix++) {
            for (int iy = 0; iy < py_grid; iy++) {
                for (int iz = 0; iz < pz_grid; iz++) {
                    da[is * nbasis + (iz + zoff)*incz1 + (iy + yoff)*incy1 + (ix + xoff)*incx1] 
                        = object[is * pbasis + ix * incx + iy*incy + iz];
                }
            }
        }
    }

    /* Call global sums to give everyone the full array */
    size_t length = (size_t)nbasis * 3;
    BlockAllreduce(da, length, pct.grid_comm);

    for(int ix=0;ix < 3 * pbasis;ix++) object[ix] = 0.0;

    double vec[3];
    for(int isy = 0; isy < nsym; isy++)
    {
        for (int ix = 0; ix < px_grid; ix++) {
            for (int iy = 0; iy < py_grid; iy++) {
                for (int iz = 0; iz < pz_grid; iz++) {

                    int idx = sym_idx[isy * pbasis + ix * incx + iy * incy + iz] ;

                    vec[0] = da[idx + 0 * nbasis];
                    vec[1] = da[idx + 1 * nbasis];
                    vec[2] = da[idx + 2 * nbasis];
                    symm_vec(isy, vec);
                    if(time_rev[isy]) 
                    {
                        vec[0] = - vec[0];
                        vec[1] = - vec[1];
                        vec[2] = - vec[2];
                    }
                    object[0*pbasis + ix * incx + iy*incy + iz] += vec[0];
                    object[1*pbasis + ix * incx + iy*incy + iz] += vec[1];
                    object[2*pbasis + ix * incx + iy*incy + iz] += vec[2];
                }
            }
        }
    }

    double t1 = (double) nsym;
    t1 = 1.0 / t1;
    for(int ix = 0; ix < 3*pbasis; ix++) object[ix] = object[ix] * t1;

    delete [] da;

}

void Symmetry::symm_vec(int isy, double *vec)
{
    double vec_tem[3], vec_rot[3];
    for (int ir = 0; ir < 3; ir++)
    {
        vec_tem[ir] = vec[0] * L.b0[ir] +vec[1] * L.b1[ir] +vec[2] * L.b2[ir];
    }                       /* end for ir */

    vec_rot[0] = 0.0;
    vec_rot[1] = 0.0;
    vec_rot[2] = 0.0;
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
            vec_rot[i] += sym_rotate[isy *9 + i* 3 + j] * vec_tem[j];

    for (int ir = 0; ir < 3; ir++)
    {
        vec[ir] = vec_rot[0] * L.a0[ir] + vec_rot[1] * L.a1[ir] + vec_rot[2] * L.a2[ir];
        if(inv_type[isy]) vec[ir] = - vec[ir];
    }                       /* end for ir */

}

void Symmetry::symmetrize_grid_object(double *object)
{
    int incx = py_grid * pz_grid;
    int incy = pz_grid;

    int incx1 = ny_grid * nz_grid;
    int incy1 = nz_grid;
    int incz1 = 1;

    // Allocate a global array object and put this processors object into the correct location
    double *da = new double[nbasis]();

    for (int ix = 0; ix < px_grid; ix++) {
        for (int iy = 0; iy < py_grid; iy++) {
            for (int iz = 0; iz < pz_grid; iz++) {
                da[(iz + zoff)*incz1 + (iy + yoff)*incy1 + (ix + xoff)*incx1] = object[ix * incx + iy*incy + iz];
            }
        }
    }

    /* Call global sums to give everyone the full array */
    int length = nbasis;
    GlobalSums ((double *)da, length, pct.grid_comm);

    for(int ix=0;ix < pbasis;ix++) object[ix] = 0.0;

    int nsym_rho = 0;
    for(int isy = 0; isy < nsym; isy++)
    {
        if(time_rev[isy]) continue;
        nsym_rho++;
        for (int ix = 0; ix < px_grid; ix++) {
            for (int iy = 0; iy < py_grid; iy++) {
                for (int iz = 0; iz < pz_grid; iz++) {

                    int idx = sym_idx[isy * pbasis + ix * incx + iy * incy + iz] ;

                    object[ix * incx + iy*incy + iz] += da[idx];
                }
            }
        }
    }

    double t1 = (double) nsym_rho;
    t1 = 1.0 / t1;
    for(int ix = 0; ix < pbasis; ix++) object[ix] = object[ix] * t1;

    delete [] da;

}

void Symmetry::symforce (void)
{


    double *force = new double[3* ct.num_ions];

    for (int ion = 0; ion < ct.num_ions; ion++)
    {
        Rmg_L.to_crystal_vector( &force[ion *3 ], Atoms[ion].force[ct.fpt[0]]);

        Atoms[ion].force[ct.fpt[0]][0] = 0.0;
        Atoms[ion].force[ct.fpt[0]][1] = 0.0;
        Atoms[ion].force[ct.fpt[0]][2] = 0.0;
    }
    for (int ion = 0; ion < ct.num_ions; ion++)
    {
        for(int isy = 0; isy < nsym_full; isy++)
        {
            int ion1 = full_sym_atom[isy * ct.num_ions + ion];
            for(int i = 0; i < 3; i++)
                for(int j = 0; j < 3; j++)
                    Atoms[ion1].force[ct.fpt[0]][i] += full_sym_rotate[isy *9 + i* 3 + j] * force[ion*3 + j];
        }

    }
    for (int ion = 0; ion < ct.num_ions; ion++)
    {
        Rmg_L.to_cartesian( Atoms[ion].force[ct.fpt[0]], &force[ion *3 ]);
    }

    for (int ion = 0; ion < ct.num_ions; ion++)
        for(int i = 0; i < 3; i++)
            Atoms[ion].force[ct.fpt[0]][i] = force[3*ion + i] /nsym_full;

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
    for(int isy = 0; isy < nsym_full; isy++)
    {

        for(int i = 0; i < 3; i++)
            for(int j = 0; j < 3; j++)
                for(int k = 0; k < 3; k++)
                    for(int l = 0; l < 3; l++)
                    {
                        mat_tensor[i*3+j] += full_sym_rotate[isy * 9 + i * 3 + k] * work[k * 3 + l] * full_sym_rotate[isy*9 + j*3 +l];
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

    for(int i = 0; i < 9; i++) mat_tensor[i] = work[i] / nsym_full;
}

void Symmetry::rotate_ylm()
{
    //dimension: 1: num of symmetry operation
    //           2: l-value
    //           3,4: m-value
    rot_ylm.resize(boost::extents[nsym][4][7][7]);


    int lmax = 3; 
    int mmax = 2*lmax + 1;
    int info, ipvt[mmax];
    double r_crds[3], r_xtal[3], r_rand[mmax][3];
    double ylm_array[mmax * mmax];
    double ylm_invert[mmax * mmax];
    double rot_tem[mmax * mmax];

    // for s orbital, no rotateion is needed.
    for(int isy = 0; isy < nsym; isy++)
        rot_ylm[isy][0][0][0] = 1.0;

    std::srand(224);
    for(int m = 0; m < mmax; m++)
    {
        r_rand[m][0] = (2.0 *std::rand())/RAND_MAX -1.0;
        r_rand[m][1] = (2.0 *std::rand())/RAND_MAX -1.0;
        r_rand[m][2] = (2.0 *std::rand())/RAND_MAX -1.0;

    }

    double sr[3][3];
    double work1[9], work2[9], b[9], a[9], one= 1.0, zero = 0.0;
    int three = 3;
    for (int i = 0; i < 3; i++)
    {
        b[0 * 3 + i] = L.b0[i];
        b[1 * 3 + i] = L.b1[i];
        b[2 * 3 + i] = L.b2[i];
        a[0 * 3 + i] = L.a0[i];
        a[1 * 3 + i] = L.a1[i];
        a[2 * 3 + i] = L.a2[i];
    }

    // get the symmetry operation in cartesian
    for(int l = 1; l < 4; l++)
    {

        int lm = 2*l + 1;

        for (int m1 = 0; m1 < 2*l+1; m1++)
        {

            L.to_cartesian(r_rand[m1], r_crds);
            for (int m2 = 0; m2 < lm; m2++)
            {
                ylm_array[m1 * lm  + m2] = Ylm(l, m2, r_rand[m1]);
            }
        }

        for(int i = 0; i < lm * lm; i++) ylm_invert[i] = 0.0;
        for(int i = 0; i < lm; i++) ylm_invert[i * lm + i] = 1.0;

        dgesv (&lm, &lm, ylm_array, &lm, ipvt, ylm_invert, &lm, &info);

        for(int isym = 0; isym < nsym; isym++)
        {
            for(int i = 0; i < 9; i++) work2[i] = sym_rotate[isym*9+i];
            dgemm("N", "N", &three, &three, &three, &one, b, &three, work2, &three, &zero, work1, &three);
            dgemm("N", "T", &three, &three, &three, &one, work1, &three, a, &three, &zero, work2, &three);
            if(pct.imgpe == 0 && ct.verbose) rmg_printf("\n rotate ylm symm op %d l=%d", isym, l);
            for(int i = 0; i < 3; i++) 
            {
                for(int j = 0; j < 3; j++) 
                {
                    sr[i][j] = work2[i*3+j];
                }
            }
            for (int m1 = 0; m1 < 2*l+1; m1++)
            {
                r_xtal[0] = sr[0][0] * r_rand[m1][0] +sr[0][1] * r_rand[m1][1] +sr[0][2] * r_rand[m1][2];
                r_xtal[1] = sr[1][0] * r_rand[m1][0] +sr[1][1] * r_rand[m1][1] +sr[1][2] * r_rand[m1][2];
                r_xtal[2] = sr[2][0] * r_rand[m1][0] +sr[2][1] * r_rand[m1][1] +sr[2][2] * r_rand[m1][2];

                for (int m2 = 0; m2 < lm; m2++)
                {
                    ylm_array[m1 * lm  + m2] = Ylm(l, m2, r_xtal);
                }

            }

            double one = 1.0, zero = 0.0;
            dgemm("T", "T", &lm, &lm, &lm, &one, ylm_invert, &lm, ylm_array, &lm, &zero, rot_tem, &lm);
            for (int m1 = 0; m1 < 2*l+1; m1++)
            {
                if(pct.imgpe == 0 && ct.verbose) rmg_printf("\n ");
                for (int m2 = 0; m2 < lm; m2++)
                {
                    rot_ylm[isym][l][m1][m2] = rot_tem[m1 * lm + m2];
                    if(pct.imgpe == 0 && ct.verbose) rmg_printf(" %7.4f ", rot_ylm[isym][l][m1][m2]);
                }
            }

        }

    }
}

void Symmetry::rotate_spin()
{
    //dimension: 1: num of symmetry operation
    //           2: l-value
    //           3,4: m-value
    rot_spin.resize(boost::extents[nsym][2][2]);
    rot_spin_wave.resize(boost::extents[nsym][2][2]);
    if(!ct.noncoll)
    {
        for(int isym = 0; isym < nsym; isym++)
        {
            rot_spin[isym][0][0] = 1.0;
            rot_spin[isym][1][0] = 0.0;
            rot_spin[isym][0][1] = 0.0;
            rot_spin[isym][1][1] = 1.0;
            rot_spin_wave[isym][0][0] = 1.0;
            rot_spin_wave[isym][1][0] = 0.0;
            rot_spin_wave[isym][0][1] = 0.0;
            rot_spin_wave[isym][1][1] = 1.0;
        }
        return;
    }
    double sr[3][3];
    double work1[9], work2[9], b[9], a[9], one= 1.0, zero = 0.0;
    int three = 3;
    double angle(0.0), tem, sint, cost, axis[3];
    double eps = 1.0e-5;
    for (int i = 0; i < 3; i++)
    {
        b[0 * 3 + i] = L.b0[i];
        b[1 * 3 + i] = L.b1[i];
        b[2 * 3 + i] = L.b2[i];
        a[0 * 3 + i] = L.a0[i];
        a[1 * 3 + i] = L.a1[i];
        a[2 * 3 + i] = L.a2[i];
    }

    for(int isym = 0; isym < nsym; isym++)
    {
        // get the symmetry operation in cartesian
        for(int i = 0; i < 9; i++) work2[i] = sym_rotate[isym*9+i];
        dgemm("N", "N", &three, &three, &three, &one, b, &three, work2, &three, &zero, work1, &three);
        dgemm("N", "T", &three, &three, &three, &one, work1, &three, a, &three, &zero, work2, &three);
        for(int i = 0; i < 3; i++) 
        {
            for(int j = 0; j < 3; j++) 
            {
                sr[i][j] = work2[i*3+j];
            }
        }
        if(pct.imgpe == 0 && ct.verbose)
        {
            rmg_printf("\n\n rotate spin for symm op %d %d %d %d", isym, type_symm(sr), (int)inv_type[isym], (int)time_rev[isym] );
            for(int i = 0; i < 3; i++) rmg_printf("\n sym %5.2f  %5.2f  %5.2f", sr[i][0], sr[i][1],sr[i][2]);
            rmg_printf("  with translation of (%d %d %d) grids ", ftau[isym*3 + 0],ftau[isym*3 + 1],ftau[isym*3 + 2]);
        }

        // improper rotation to proper rotation
        if(type_symm(sr) == 5 || type_symm(sr) == 6) 
        {
            for(int i = 0; i < 3; i++) 
            {
                for(int j = 0; j < 3; j++) 
                {
                    sr[i][j] = -sr[i][j];
                }
            }
        }

        if(type_symm(sr) == 1 ||type_symm(sr) == 2)
        {
            rot_spin[isym][0][0] = 1.0;
            rot_spin[isym][0][1] = 0.0;
            rot_spin[isym][1][0] = 0.0;
            rot_spin[isym][1][1] = 1.0;
            rot_spin_wave[isym][0][0] = 1.0;
            rot_spin_wave[isym][0][1] = 0.0;
            rot_spin_wave[isym][1][0] = 0.0;
            rot_spin_wave[isym][1][1] = 1.0;
            angle = 0.0;
            axis[0] = 0.0;
            axis[1] = 0.0;
            axis[2] = 0.0;
        }
        else if(type_symm(sr) == 4)
        {
            angle = PI;
            for(int i = 0; i < 3; i++)
                axis[i] = std::sqrt( std::abs( sr[i][i] + 1.0)/2.0);

            for(int i = 0; i < 3; i++)
            {
                for(int j = i+1; j < 3; j++)
                {
                    if(std::abs(axis[i] * axis[j]) > eps )
                    {
                        axis[i] = 0.5 * sr[i][j]/axis[j];
                    }
                }
            }

        }
        else if(type_symm(sr) == 3)
        {
            axis[0] = -sr[1][2] + sr[2][1];
            axis[1] = -sr[2][0] + sr[0][2];
            axis[2] = -sr[0][1] + sr[1][0];

            tem =  std::sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
            sint = 0.5 * tem;
            if(tem < eps || sint > 1.0 +eps)
            {
                rmg_printf("\n rotation matrix error: norm = %f \n", tem);
                fflush(NULL);
                rmg_error_handler(__FILE__, __LINE__, " type of symmetry not defined Exiting.\n");
            }
            axis[0] /= tem;
            axis[1] /= tem;
            axis[2] /= tem;
            if(1.0 - axis[0] * axis[0] > eps)
            {
                cost = (sr[0][0]- axis[0] * axis[0])/(1.0 -axis[0] * axis[0]);
            }
            else if (1.0 - axis[1] * axis[1] > eps)
            {
                cost = (sr[1][1]- axis[1] * axis[1])/(1.0 -axis[1] * axis[1]);
            }
            else 
            {
                cost = (sr[2][2]- axis[2] * axis[2])/(1.0 -axis[2] * axis[2]);
            }

            angle = std::atan2(sint,cost);

        }
        else
        {
            rmg_printf("\n rotation matrix error \n");
            fflush(NULL);
            rmg_error_handler(__FILE__, __LINE__, " type of symmetry not defined Exiting.\n");
        }

        angle = 0.5 * angle;
        sint = std::sin(angle);
        cost = std::cos(angle);
        rot_spin[isym][0][0] = std::complex<double>(cost, -axis[2] * sint);
        rot_spin[isym][0][1] = std::complex<double>(-axis[1] * sint, -axis[0] * sint);
        rot_spin[isym][1][0] = -std::conj(rot_spin[isym][0][1]);
        rot_spin[isym][1][1] = std::conj(rot_spin[isym][0][0]);
        rot_spin_wave[isym][0][0] = std::complex<double>(cost, -axis[2] * sint);
        rot_spin_wave[isym][0][1] = std::complex<double>(-axis[1] * sint, -axis[0] * sint);
        rot_spin_wave[isym][1][0] = -std::conj(rot_spin[isym][0][1]);
        rot_spin_wave[isym][1][1] = std::conj(rot_spin[isym][0][0]);


        if(time_rev[isym])
        {
            rot_spin[isym][0][0] = -std::complex<double>(-axis[1] * sint, -axis[0] * sint);
            rot_spin[isym][0][1] = std::complex<double>(cost, -axis[2] * sint);
            rot_spin[isym][1][0] = std::conj(rot_spin[isym][0][1]);
            rot_spin[isym][1][1] = std::conj(rot_spin[isym][0][0]);
        }


        if(pct.imgpe == 0 && ct.verbose)
        {
            rmg_printf("\n angle: %7.4f  axis: %7.4f %7.4f  %7.4f", angle, axis[0], axis[1], axis[2]);
            rmg_printf("\n (%7.4f  %7.4f)  ", std::real(rot_spin[isym][0][0]), std::imag(rot_spin[isym][0][0])); 
            rmg_printf(  " (%7.4f  %7.4f)  ", std::real(rot_spin[isym][0][1]), std::imag(rot_spin[isym][0][1])); 
            rmg_printf("\n (%7.4f  %7.4f)  ", std::real(rot_spin[isym][1][0]), std::imag(rot_spin[isym][1][0])); 
            rmg_printf(  " (%7.4f  %7.4f)  ", std::real(rot_spin[isym][1][1]), std::imag(rot_spin[isym][1][1])); 
        }

    }
}

int Symmetry::type_symm(double sr[3][3])
{
    // follow the QE notation for symmetry type
    // 1: identity, 2: inversion, 3, proper rotation,, 4 180 degree rotation,, 5: mirror, 6: improper rotation

    double eps = 1.0e-5;
    double det, det1, sum;
    sum = 0.0;
    for(int i = 0; i < 3; i++) 
    {
        for(int j = 0; j < 3; j++) 
        {
            if(i == j) sum += std::abs(sr[i][j]-1.0);
            else    sum += std::abs(sr[i][j]);
        }
    }
    if( std::abs(sum) < eps ) return 1;
    sum = 0.0;
    for(int i = 0; i < 3; i++) 
        for(int j = 0; j < 3; j++) 
        {
            if(i == j) sum += std::abs(sr[i][j]+1.0);
            else    sum += std::abs(sr[i][j]);
        }
    if( std::abs(sum) < eps ) return 2;

    det = sr[0][0] * ( sr[1][1] * sr[2][2] - sr[2][1] * sr[1][2] ) -  
        sr[0][1] * ( sr[1][0] * sr[2][2] - sr[2][0] * sr[1][2] ) +  
        sr[0][2] * ( sr[1][0] * sr[2][1] - sr[2][0] * sr[1][1] );

    // determinant of sr  = 1.0: proper rotation =-1 improper rotation (mirror symmetry)
    if(std::abs(det - 1.0) < eps)
    {
        // eigenvalue = -1: 180 degree rotation 
        det1 = (sr[0][0]+1.0) * ( (sr[1][1]+1.0) * (sr[2][2]+1.0) - sr[2][1] * sr[1][2] ) -  
            sr[0][1] * ( sr[1][0] * (sr[2][2]+1.0) - sr[2][0] * sr[1][2] ) +  
            sr[0][2] * ( sr[1][0] * sr[2][1] - sr[2][0] * (sr[1][1]+1.0) );
        if(std::abs(det1) < eps) return 4;
        else return 3;
    }
    else if(std::abs(det + 1.0) < eps)
    {
        // eigenvalue = +1: mirror symmetry
        det1 = (sr[0][0]-1.0) * ( (sr[1][1]-1.0) * (sr[2][2]-1.0) - sr[2][1] * sr[1][2] ) -  
            sr[0][1] * ( sr[1][0] * (sr[2][2]+1.0) - sr[2][0] * sr[1][2] ) +  
            sr[0][2] * ( sr[1][0] * sr[2][1] - sr[2][0] * (sr[1][1]+1.0) );
        if(std::abs(det1) < eps) return 5;
        else return 6;

    }
    else
    {
        rmg_printf("\n determinant of the rotation matrix %f\n", det);
        fflush(NULL);
        rmg_error_handler(__FILE__, __LINE__, " type of symmetry not defined Exiting.\n");
    }
    return 0; 

}
void Symmetry::nsocc_AFM(boost::multi_array<std::complex<double>, 4> &ns_occ, int mmax, std::vector<int> map_to_ldaU_ion, std::vector<int> ldaU_ion_index)
{
    int num_ldaU_ions = ldaU_ion_index.size();


    for(int ion = 0; ion < num_ldaU_ions; ion++)
    {
        for(int i = 0; i< mmax; i++)
        {
            for(int j = 0; j < mmax; j++)
            {
                ns_occ[1][ion][i][j] = 0.0;
            }
        }
    }


    //  the loops below can be optimized if it is slow    
    for (int ion = 0; ion < num_ldaU_ions; ion++)
    {

        int ion_idx = ldaU_ion_index[ion];
        int num_orb = Species[Atoms[ion_idx].species].num_ldaU_orbitals;
        int l_val = Species[Atoms[ion_idx].species].ldaU_l;
        for(int i1 = 0; i1 < num_orb; i1++)
        {
            for(int i2 = 0; i2 < num_orb; i2++)
            { 
                for(int isy = 0; isy < nsym; isy++)
                {
                    if(!time_rev[isy] ) continue;
                    int ion1_idx = sym_atom[isy * ct.num_ions + ion_idx];
                    int ion1 = map_to_ldaU_ion[ion1_idx];

                    for(int i3 = 0; i3 < num_orb; i3++)
                    {
                        for(int i4 = 0; i4 < num_orb; i4++)
                        { 
                            ns_occ[1][ion][i1][i2] += rot_ylm[isy][l_val][i1][i3] *
                                ns_occ[0][ion1][i4][i3] * rot_ylm[isy][l_val][i2][i4];
                        }
                    }

                }
            }
        }
    }

    for(int ion = 0; ion < num_ldaU_ions; ion++)
    {
        for(int i = 0; i< mmax; i++)
        {
            for(int j = 0; j < mmax; j++)
            {
                ns_occ[1][ion][i][j] /= (double)n_time_rev;
            }
        }
    }
}

void Symmetry::symm_nsocc(std::complex<double> *ns_occ_g, int mmax, std::vector<int> map_to_ldaU_ion, std::vector<int> ldaU_ion_index)
{
    boost::multi_array<std::complex<double>, 5> ns_occ_sum, ns_occ;
    int num_spin = 1;
    if (ct.nspin != 1) num_spin = 2;
    int num_ldaU_ions = ldaU_ion_index.size();
    ns_occ_sum.resize(boost::extents[num_spin][num_spin][num_ldaU_ions][mmax][mmax]);
    ns_occ.resize(boost::extents[num_spin][num_spin][num_ldaU_ions][mmax][mmax]);

    int size_each_spin = num_ldaU_ions * mmax * mmax;
    if(ct.nspin == 1)
    {
        for(int idx = 0; idx < size_each_spin; idx++)
        {
            ns_occ.data()[idx] = ns_occ_g[idx];
        }

    }
    else if(ct.nspin == 2)
    {
        for(int idx = 0; idx < size_each_spin; idx++)
        {
            ns_occ.data()[idx] = ns_occ_g[idx];
            ns_occ.data()[idx + size_each_spin] = 0.0;
            ns_occ.data()[idx + 2*size_each_spin] = 0.0;
            ns_occ.data()[idx + 3*size_each_spin] = ns_occ_g[idx + size_each_spin];
        }

    }
    else 
    {
        for(int idx = 0; idx < size_each_spin * ct.nspin; idx++)
        {
            ns_occ.data()[idx] = ns_occ_g[idx];
        }
    }

    size_t occ_size = num_spin * num_spin * num_ldaU_ions * mmax * mmax;
    for(size_t idx = 0; idx < occ_size; idx++)ns_occ_sum.data()[idx] = 0.0;
    //  the loops below can be optimized if it is slow    
    for (int ion = 0; ion < num_ldaU_ions; ion++)
    {

        int ion_idx = ldaU_ion_index[ion];
        int num_orb = Species[Atoms[ion_idx].species].num_ldaU_orbitals;
        int l_val = Species[Atoms[ion_idx].species].ldaU_l;
        for(int is1 = 0; is1 < num_spin; is1++)
        {
            for(int is2 = 0; is2 < num_spin; is2++)
            {
                for(int i1 = 0; i1 < num_orb; i1++)
                {
                    for(int i2 = 0; i2 < num_orb; i2++)
                    { 
                        for(int isy = 0; isy < nsym; isy++)
                        {
                            int ion1_idx = sym_atom[isy * ct.num_ions + ion_idx];
                            int ion1 = map_to_ldaU_ion[ion1_idx];

                            for(int is3 = 0; is3 < num_spin; is3++)
                            {     
                                for(int is4 = 0; is4 < num_spin; is4++)
                                { 
                                    for(int i3 = 0; i3 < num_orb; i3++)
                                    {
                                        for(int i4 = 0; i4 < num_orb; i4++)
                                        { 
                                            if(time_rev[isy] )
                                            {
                                                if(ct.nspin == 4)
                                                {
                                                    ns_occ_sum[is1][is2][ion][i1][i2] += 
                                                        std::conj(rot_spin[isy][is1][is3]) * rot_ylm[isy][l_val][i1][i3] *
                                                        ns_occ[is4][is3][ion1][i4][i3] *
                                                        rot_spin[isy][is2][is4] * rot_ylm[isy][l_val][i2][i4];
                                                }
                                            }
                                            else
                                            {
                                                ns_occ_sum[is1][is2][ion][i1][i2] += 
                                                    std::conj(rot_spin[isy][is1][is3]) * rot_ylm[isy][l_val][i1][i3] *
                                                    ns_occ[is3][is4][ion1][i3][i4] *
                                                    rot_spin[isy][is2][is4] * rot_ylm[isy][l_val][i2][i4];
                                            }

                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if(ct.nspin == 1)
    {
        for(int idx = 0; idx < size_each_spin; idx++)
        {
            ns_occ_g[idx] = ns_occ_sum.data()[idx]/(double)nsym;
        }

    }
    else if(ct.nspin == 2)
    {
        for(int idx = 0; idx < size_each_spin; idx++)
        {
            ns_occ_g[idx] = ns_occ_sum.data()[idx]/(double)(nsym-n_time_rev);
            ns_occ_g[idx + size_each_spin] = ns_occ_sum.data()[idx + 3*size_each_spin]/(double)(nsym-n_time_rev);
        }

    }
    else 
    {
        for(int idx = 0; idx < size_each_spin * ct.nspin; idx++)
        {
            ns_occ_g[idx] = ns_occ_sum.data()[idx]/(double)nsym;
        }
    }

}

void Symmetry::setgrid(BaseGrid &G, int density)
{
    px_grid = G.get_PX0_GRID(density);
    py_grid = G.get_PY0_GRID(density);
    pz_grid = G.get_PZ0_GRID(density);

    nx_grid = G.get_NX_GRID(density);
    ny_grid = G.get_NY_GRID(density);
    nz_grid = G.get_NZ_GRID(density);

    max_pdim = std::max(nx_grid, ny_grid);
    max_pdim = std::max(max_pdim, nz_grid);
    //int nx_grid_c = G.get_NX_GRID(1);
    //int ny_grid_c = G.get_NY_GRID(1);
    //int nz_grid_c = G.get_NZ_GRID(1);

    xoff = G.get_PX_OFFSET(density);
    yoff = G.get_PY_OFFSET(density);
    zoff = G.get_PZ_OFFSET(density);

    pbasis = px_grid * py_grid * pz_grid;
    nbasis = nx_grid * ny_grid * nz_grid;
    // sym index arrays dimensioned to size of smallest possible integer type
    sym_idx.resize(nsym * pbasis);
    init_symm_ijk();

    ct.nsym = nsym;
}
Symmetry::~Symmetry(void)
{
}


void Symmetry::symmetrize_rho_AFM(double *rho,double *rho_oppo)
{
    int incx = py_grid * pz_grid;
    int incy = pz_grid;

    int incx1 = ny_grid * nz_grid;
    int incy1 = nz_grid;
    int incz1 = 1;

    // Allocate a global array object and put this processors object into the correct location
    double *da1 = new double[nbasis]();

    for (int ix = 0; ix < px_grid; ix++) {
        for (int iy = 0; iy < py_grid; iy++) {
            for (int iz = 0; iz < pz_grid; iz++) {
                da1[ (iz + zoff)*incz1 + (iy + yoff)*incy1 + (ix + xoff)*incx1] 
                    = rho[ ix * incx + iy*incy + iz];
            }
        }
    }

    /* Call global sums to give everyone the full array */
    size_t length = (size_t)nbasis;
    BlockAllreduce(da1, length, pct.grid_comm);

    for(int idx = 0; idx < px_grid * py_grid * pz_grid; idx++) rho_oppo[idx] = 0.0;
    for(int isy = 0; isy < nsym; isy++)
    {
        if(time_rev[isy]) 
        {
            for (int ix = 0; ix < px_grid; ix++) {
                for (int iy = 0; iy < py_grid; iy++) {
                    for (int iz = 0; iz < pz_grid; iz++) {

                        int idx = sym_idx[isy * pbasis + ix * incx + iy * incy + iz] ;
                        rho_oppo[ix * incx + iy*incy + iz] += da1[idx];
                    }
                }
            }
        }
    }

    for(int idx = 0; idx < px_grid * py_grid * pz_grid; idx++) rho_oppo[idx] /= (double)n_time_rev;
    delete [] da1;

}

