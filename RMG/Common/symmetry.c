/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/symmetry.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5
 * COPYRIGHT
 *   Copyright (C) 1995  Emil Briggs
 *   Copyright (C) 1998  Emil Briggs, Charles Brabec, Mark Wensell, 
 *                       Dan Sullivan, Chris Rapcewicz, Jerzy Bernholc
 *   Copyright (C) 2001  Emil Briggs, Wenchang Lu,
 *                       Marco Buongiorno Nardelli,Charles Brabec, 
 *                       Mark Wensell,Dan Sullivan, Chris Rapcewicz,
 *                       Jerzy Bernholc
 * FUNCTION
 *   C wrapper routines to call fortran functions which actually do the
 *   symmetrization.
 * 
 *   void init_sym(void)
 *     This routine is used to initialize the symmetry structures 
 *   void symmetrize_rho(P0_GRID *rho)
 *     Symmetrizes the density
 *   void symforce(void)
 *     This routine is used to symmetrize the forces (MBN)
 * INPUTS
 *   nothing
 * OUTPUT
 *   charge density and forces are updated 
 * PARENTS
 *   force.c get_rho.c init.c
 * CHILDREN
 *   to_crystal.c symmetry.f
 * SOURCE
 */

#include <float.h>
#include <math.h>
#include <assert.h>

#include "main.h"


static int *s;
static int irg[MAX_SYMMETRY], irt[MAX_IONS][MAX_SYMMETRY];
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
    int ion, kpt, wflag;
    int ibrav = get_ibrav_type();
    double a0[3], a1[3], a2[3], omega;
    double frac1, frac2, frac3, intpart, symprec = 1.0e-5;
    double *tau, *translation;
    int *ityp, *sa;


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

    a0[0] = get_a0(0);
    a0[1] = get_a0(1);
    a0[2] = get_a0(2);
    a1[0] = get_a1(0);
    a1[1] = get_a1(1);
    a1[2] = get_a1(2);
    a2[0] = get_a2(0);
    a2[1] = get_a2(1);
    a2[2] = get_a2(2);
    omega = get_omega();

    /* This function uses MAX_IONS as a limit for array sizes.
     * It is, of course, possible to allocate these arrays dynamically,
     * but I am not sure if this would not break FORTRAN routines to which
     * those arrays are passed*/

    if (ct.num_ions >= MAX_IONS)
        error_handler ("Too many ions, increase MAX_IONS in params.h");

    nr1 = get_FNX_GRID();
    nr2 = get_FNY_GRID();
    nr3 = get_FNZ_GRID();

    /* Only have PE zero output symmetry information */
    wflag = pct.gridpe;


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

    nsym = 0;
    for(kpt = 0; kpt < nsym_atom; kpt++)
    {
        frac1 = modf(translation[kpt*3 + 0] * nr1, &intpart);
        frac2 = modf(translation[kpt*3 + 1] * nr2, &intpart);
        frac3 = modf(translation[kpt*3 + 2] * nr3, &intpart);
        if(frac1 < symprec && frac2 < symprec &&frac3 < symprec)
        {

            printf("\n sym operation after considering real space grid # %d",nsym);
            for(i = 0; i < 3; i++)
                for(j = 0; j < 3; j++)
                    s[nsym * 9 + i *3 + j] = sa[kpt * 9 + i *3 + j];

            ftau[nsym*3 + 0] = translation[kpt*3 + 0] * nr1;
            ftau[nsym*3 + 1] = translation[kpt*3 + 1] * nr2;
            ftau[nsym*3 + 2] = translation[kpt*3 + 2] * nr3;

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

    assert(nsym >0);
    if(nsym == 1) ct.is_use_symmetry = 0;


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
    int FP0_BASIS = get_FP0_BASIS();
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
    int ion, isy, i, j;

    double force[3];


    for (ion = 0; ion < ct.num_ions; ion++)
    {
        for(i = 0; i < 3; i++) force[i] = 0.0;
        for(isy = 0; isy < nsym; isy++)
        for(i = 0; i < 3; i++)
        for(j = 0; j < 3; j++)
            force[i] += s[isy *9 + i* 3 + j] * ct.ions[ion].force[ct.fpt[0]][j];

       
        for(i = 0; i < 3; i++) 
            ct.ions[ion].force[ct.fpt[0]][i] = force[i] /nsym;

    }

}                               /* end symforce */

void  symm_ijk(int *srotate, int *strans, int ix, int iy, int iz, int *ixx, int *iyy, int *izz,
              int nx, int ny, int nz)
{


    // rotate ijk by symmetry, and then translation
    int i, j, k, ipoint[3], opoint[3] ; 

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
