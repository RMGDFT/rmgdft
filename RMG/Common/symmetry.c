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
#include <assert.h>

#include "main.h"


static int *s;
//static int irg[MAX_SYMMETRY], irt[MAX_IONS][MAX_SYMMETRY];
static int *sym_atom;  //  atom B = sym(atom A)
static int *ftau;
static int nsym;


int spg_get_multiplicity(const double *lattice,
                         const double *position,
                         const int *types,
                         const int num_atom,
                         const double symprec);

int spg_get_symmetry(int *rotation,
                     double *translation,
                     const int max_size,
                     const double *lattice,
                     const double *position,
                     const int types[],
                     const int num_atom,
                     const double symprec);

void  symm_ijk(int *srotate, int *strans, int ix, int iy, int iz, int *ixx, int *iyy, int *izz,
              int nx, int ny, int nz);

/* This routine is used to initialize the symmetry structures */
void init_sym (void)
{
    int nr1, nr2, nr3;
    int ion, kpt;
    double frac1, frac2, frac3, intpart, symprec = 1.0e-5;
    double *tau, *translation;
    double xtal[3];
    double ndim[3];
    int *ityp, *sa;
    //double omega = get_omega();
    //int ibrav = get_ibrav_type();
    /* Only have PE zero output symmetry information */
    //int wflag = pct.gridpe;



    double lattice[9];
    int nsym_atom;
    int i, j;

    lattice[0*3+0] = get_a0(0);
    lattice[0*3+1] = get_a0(1);
    lattice[0*3+2] = get_a0(2);
    lattice[1*3+0] = get_a1(0);
    lattice[1*3+1] = get_a1(1);
    lattice[1*3+2] = get_a1(2);
    lattice[2*3+0] = get_a2(0);
    lattice[2*3+1] = get_a2(1);
    lattice[2*3+2] = get_a2(2);


    /* This function uses MAX_IONS as a limit for array sizes.
     * It is, of course, possible to allocate these arrays dynamically,
     * but I am not sure if this would not break FORTRAN routines to which
     * those arrays are passed*/

    if (ct.num_ions >= MAX_IONS)
        error_handler ("Too many ions, increase MAX_IONS in params.h");

    nr1 = get_FNX_GRID();
    nr2 = get_FNY_GRID();
    nr3 = get_FNZ_GRID();
    ndim[0] = nr1;
    ndim[1] = nr2;
    ndim[2] = nr3;


    my_malloc (tau, 3 * ct.num_ions, double);
    my_malloc (ityp, ct.num_ions, int);

    /* Set up atomic positions and species for fortran routines */
    for (ion = 0; ion < ct.num_ions; ion++)
    {
        to_crystal (&tau[ion * 3], ct.ions[ion].crds);
        ityp[ion] = ct.ions[ion].species;
    }


    nsym_atom = spg_get_multiplicity(lattice, tau, ityp, ct.num_ions, symprec);


    my_malloc (sa, 9 * nsym_atom, int);
    my_malloc (translation, 3 * nsym_atom, double);
    my_malloc (s, 9 * nsym_atom, int);
    my_malloc (ftau, 3 * nsym_atom, int);

    my_malloc (ct.sym_trans, 3 * nsym_atom, double);
    my_malloc (ct.sym_rotate, 9 * nsym_atom, int);

    nsym_atom = spg_get_symmetry(sa, translation,  nsym_atom, lattice, tau, ityp, ct.num_ions, symprec);

//    for(kpt = 0; kpt < nsym_atom; kpt++)
//    {
 //       printf("\n sym operation considering atom only # %d", kpt);
  //      for(i = 0; i < 3; i++)
   //     {
    //        printf("\n      %3d  %3d  %3d", sa[kpt][i][0],sa[kpt][i][1],sa[kpt][i][2]);
     //   }
      //  printf("\n  translation in crystal unit:    %f %f %f ", translation[kpt][0], translation[kpt][1], translation[kpt][2]);
    //}
    //   check if the real space grid fit the symmetry operations, if not, kick it out

    for(kpt = 0; kpt < nsym_atom; kpt++)
    {
        frac1 = modf(translation[kpt*3 + 0] * nr1, &intpart);
        frac2 = modf(translation[kpt*3 + 1] * nr2, &intpart);
        frac3 = modf(translation[kpt*3 + 2] * nr3, &intpart);
        if(frac1 > 0.5) frac1 = 1.0-frac1;
        if(frac2 > 0.5) frac2 = 1.0-frac2;
        if(frac3 > 0.5) frac3 = 1.0-frac3;
        if(frac1 < symprec && frac2 < symprec &&frac3 < symprec)
        {

            printf("\n sym operation after considering real space grid # %d",nsym);
            for(i = 0; i < 3; i++)
                for(j = 0; j < 3; j++)
                {
                    s[nsym * 9 + i *3 + j] = sa[kpt * 9 + i *3 + j];
                    ct.sym_rotate[nsym * 9 + i *3 + j] = sa[kpt * 9 + i *3 + j];
                }

            ftau[nsym*3 + 0] = translation[kpt*3 + 0] * nr1 + symprec;
            ftau[nsym*3 + 1] = translation[kpt*3 + 1] * nr2 + symprec;
            ftau[nsym*3 + 2] = translation[kpt*3 + 2] * nr3 + symprec;

            ct.sym_trans[nsym*3+0] = translation[kpt*3+0];
            ct.sym_trans[nsym*3+1] = translation[kpt*3+1];
            ct.sym_trans[nsym*3+2] = translation[kpt*3+2];
            


            for(i = 0; i < 3; i++)
            {
                printf("\n      %3d  %3d  %3d", s[nsym * 9 + i *3 + 0],s[nsym * 9 + i *3 + 1],s[nsym * 9 + i *3 + 2]);
            }
            printf("  with translation of (%d %d %d) grids ", ftau[nsym*3 + 0],ftau[nsym*3 + 1],ftau[nsym*3 + 2]);
            nsym++;
        }
        else
        {
//            printf("\n translation break a symmetry") ;
 //           for(i = 0; i < 3; i++)
  //          {
   //             printf("\n      %3d  %3d  %3d", sa[kpt * 9 + i *3 + 0],sa[kpt * 9 + i *3 + 1],sa[kpt * 9 + i *3 + 2]);
    //        }
     //       printf("  with translation of (%f %f %f) grids ", frac1, frac2, frac3);
        }
    }

    ct.nsym = nsym;
    printf("\n number of sym operation before considering real space grid: %d",nsym_atom);
    printf("\n number of sym operation  after considering real space grid: %d",nsym);
    assert(nsym >0);
    if(nsym == 1) ct.is_use_symmetry = 0;

    my_malloc (sym_atom,  ct.num_ions* nsym, int);

//  determine equivenlent ions after symmetry operation
    bool find_atom, cond_x, cond_y, cond_z;
    double mod_x, mod_y, mod_z;
    int ionb, isym;
    for(isym = 0; isym < nsym; isym++)
    {
        for (ion = 0; ion < ct.num_ions; ion++)
        {

            // xtal is the coordinates of atom operated by symmetry operatio isym.
            for(i = 0; i < 3; i++)
            {
                xtal[i] = s[isym *9 + i *3 + 0] * ct.ions[ion].xtal[0]
                    + s[isym *9 + i *3 + 1] * ct.ions[ion].xtal[1]
                    + s[isym *9 + i *3 + 2] * ct.ions[ion].xtal[2] +ftau[isym *3 + i]/ndim[i];

                if(xtal[i] + symprec  < 0.0) xtal[i]= xtal[i] + 1.0;
                if(xtal[i] + symprec >= 1.0) xtal[i]= xtal[i] - 1.0;
            }

            find_atom = false;
            cond_x = false;
            cond_y = false;
            cond_z = false;
            for (ionb = 0; ionb < ct.num_ions; ionb++)
            {
                if(ityp[ion] == ityp[ionb])
                {
//                    r =  (xtal[0] - ct.ions[ionb].xtal[0]) *(xtal[0] - ct.ions[ionb].xtal[0])
//                        +(xtal[1] - ct.ions[ionb].xtal[1]) *(xtal[1] - ct.ions[ionb].xtal[1])
//                        +(xtal[2] - ct.ions[ionb].xtal[2]) *(xtal[2] - ct.ions[ionb].xtal[2]);
                    mod_x = (xtal[0] - ct.ions[ionb].xtal[0]) *(xtal[0] - ct.ions[ionb].xtal[0]);
                    mod_y = (xtal[1] - ct.ions[ionb].xtal[1]) *(xtal[1] - ct.ions[ionb].xtal[1]);
                    mod_z = (xtal[2] - ct.ions[ionb].xtal[2]) *(xtal[2] - ct.ions[ionb].xtal[2]);

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
            {   printf("\n Equivalent atom not found %d %d %d %d %e %e %e \n", ion, ionb,isym, nsym, xtal[0],xtal[1], xtal[2]);
                fflush(NULL);
                exit(0);
            }
        }
    }


    my_free (tau);
    my_free (ityp);
    my_free (sa);
    my_free (translation);

    //    if (ct.boundaryflag != CLUSTER)
    //    {
    //       /* Call the symmetry routine */
    //      symmetry (&ibrav, &s[0][0][0], &nsym, irg, &irt[0][0],
    //             &ftau[0][0], &ct.num_ions, &tau[0][0], ityp,
    //            &ct.num_kpts, xk, wk, celldm, &nr1, &nr2, &nr3, a0, a1, a2, &omega, &wflag);
    //}

}                               /* end init_sym */




/* Symmetrizes the density */
void symmetrize_rho (double * rho)
{

    int idx, ix, iy, iz, xoff, yoff, zoff;
    int ix1, iy1, iz1;
    int isy, ixx, iyy, izz;
    double *da;
    double t1, tem;

    int FPX0_GRID = get_FPX0_GRID();
    int FPY0_GRID = get_FPY0_GRID();
    int FPZ0_GRID = get_FPZ0_GRID();
    int FNX_GRID = get_FNX_GRID();
    int FNY_GRID = get_FNY_GRID();
    int FNZ_GRID = get_FNZ_GRID();
    int FN_BASIS = FNX_GRID * FNY_GRID * FNZ_GRID;

    /* Wait until all processors arrive at this point */
    /*my_barrier(); */

    if(!ct.is_use_symmetry) return;


    /* Get some memory */
    my_malloc (da, FN_BASIS, double);


    for(idx = 0;idx < FN_BASIS;idx++) {
        da[idx] = 0.0;
    }



    /* Put this processors charge in the correct place */
    xoff = get_FPX_OFFSET();
    yoff = get_FPY_OFFSET();
    zoff = get_FPZ_OFFSET();

    int incx = FPY0_GRID * FPZ0_GRID;
    int incy = FPZ0_GRID;

    int incx1 = 1;
    int incy1 = FNX_GRID;
    int incz1 = FNX_GRID * FNY_GRID;

    for (ix = 0; ix < FPX0_GRID; ix++) {
        for (iy = 0; iy < FPY0_GRID; iy++) {
            for (iz = 0; iz < FPZ0_GRID; iz++) {
                da[(iz + zoff)*incz1 + (iy + yoff)*incy1 + (ix + xoff)*incx1] = rho[ix * incx + iy*incy + iz];
            }
        }
    }

    /* Call global sums to give everyone the full array */
    global_sums (da, &FN_BASIS, pct.grid_comm);


    /* Do the symmetrization on this processor */

    //    symrho (da, &FNX_GRID, &FNY_GRID, &FNZ_GRID, &nsym, &s[0][0][0], irg, &ftau[0][0]);


    /* Pack density back into correct place */
    t1 = (double) nsym;
    t1 = 1.0 / t1;

    for (ix = 0; ix < FPX0_GRID; ix++) {
        for (iy = 0; iy < FPY0_GRID; iy++) {
            for (iz = 0; iz < FPZ0_GRID; iz++) {
                tem = 0.0;
                for(isy = 0; isy < nsym; isy++)
                {
                    ix1 = ix + xoff;
                    iy1 = iy + yoff;
                    iz1 = iz + zoff;

                    symm_ijk(&s[isy *9], &ftau[isy*3], ix1, iy1, iz1, &ixx, &iyy, &izz, 
                            FNX_GRID, FNY_GRID, FNZ_GRID);
                    tem += da[izz *incz1 + iyy *incy1 + ixx *incx1];
                }
                rho[ix * incx + iy*incy + iz] = tem *t1;

            }
        }
    }


    /* Release our memory */
    my_free (da);


}                               /* end symmetrize_rho */




void symforce ()
{

    int ion, isy, i, j, ion1;

    double *force;

    my_malloc (force, 3* ct.num_ions, double);

    for(i = 0; i < 3 * ct.num_ions; i++ ) force[i] = 0.0;
    for (ion = 0; ion < ct.num_ions; ion++)
    {
        for(isy = 0; isy < nsym; isy++)
        {
            ion1 = sym_atom[isy * ct.num_ions + ion];
            for(i = 0; i < 3; i++)
                for(j = 0; j < 3; j++)
                    force[3* ion1 + i] += s[isy *9 + i* 3 + j] * ct.ions[ion].force[ct.fpt[0]][j];
        }



    }

    for (ion = 0; ion < ct.num_ions; ion++)
        for(i = 0; i < 3; i++) 
            ct.ions[ion].force[ct.fpt[0]][i] = force[3*ion + i] /nsym;
    my_free(force);
}                               /* end symforce */

void  symm_ijk(int *srotate, int *strans, int ix, int iy, int iz, int *ixx, int *iyy, int *izz,
        int nx, int ny, int nz)
{


    // rotate ijk by symmetry, and then translation
    int i, j, ipoint[3], opoint[3] ; 

    ipoint[0] = ix;
    ipoint[1] = iy;
    ipoint[2] = iz;

    for(i = 0; i < 3; i++)
    {
        opoint[i] = 0;
        for(j = 0; j < 3; j++)
            opoint[i] += srotate[i *3 + j] * ipoint[j] ;
        opoint[i] += strans[i];
    }



    *ixx = (opoint[0] + nx)%nx;
    *iyy = (opoint[1] + ny)%ny;
    *izz = (opoint[2] + nz)%nz;


}
/******/
