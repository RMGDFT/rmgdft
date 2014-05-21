/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/md.c *****
 * NAME
 *   Ab initio O(n) real space code with 
 *   localized orbitals and multigrid acceleration
 *   Version: 3.0.0
 * COPYRIGHT
 *   Copyright (C) 2001  Wenchang Lu,
 *                       Jerzy Bernholc
 * FUNCTION
 *   int main(int argc, char **argv)
 *   Main program
 *   Read-in all informations, structures, pseudopentials, etc. 
 *   Then enters the main driver loop. 
 * INPUTS
 *   when we run it, we need to give the input control 
 *   file name in the first argument
 *   for example, md in.diamond8
 * OUTPUT
 *   nothing
 * PARENTS
 *   This is grand-grand-....
 * CHILDREN
 *   run.c
 * SEE ALSO
 *   main.h for structure defination
 * SOURCE
 */

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include "main.h"
//include "svnrev.h"



/* Main control structure which is declared extern in main.h so any module */
/* may access it.					                 */
//CONTROL ct;

/* PE control structure which is also declared extern in main.h */
//PE_CONTROL pct;

void eldyn_(int *num_states, double *, double *, double *, double *, int *, int*);

int main(int argc, char **argv)
{



    time_t tt;
    char *timeptr;

    double *Hmatrix, *Smatrix, *Xij, *Yij, *Zij;
    int Ieldyn=1, iprint = 6;
    int n2, num_states,i;
    double *Pn0, *Pn1;

    time(&tt);
    timeptr = ctime(&tt);
    ct.time0 = my_crtc();

    ct.images_per_node = 1;
    init_IO(argc, argv);


    my_barrier();

    

    /*  Begin to do the real calculations */
    init_TDDFT();
    get_cholesky_real(matB);
    get_dm_diag_p(states, l_s, mat_X, Hij);

    write_eigs(states);
    num_states = ct.num_states;
    n2 = num_states * num_states;
    my_malloc( Pn0, 2*n2, double );
    my_malloc( Pn1, 2*n2, double );

    for(i = 0; i < n2; i++) Pn0[i] = mat_X[i];
    for(i = 0; i < n2; i++) Pn0[i+n2] = 0.0;

/*  matB: overlap matrix, Hij:  Hamiltonian matrix */
    int j, k, idx;

    double x, y, z, dipole_ion_x = 0.0;
    double dipole_ion_y = 0.0;
    double dipole_ion_z = 0.0;

    double dipole_ele_x = 0.0;
    double dipole_ele_y = 0.0;
    double dipole_ele_z = 0.0;

    dipole_ele_x = 0.0;
    dipole_ele_y = 0.0;
    dipole_ele_z = 0.0;
    for(i = 0; i < pct.FPX0_GRID; i++)
    {
        x = (pct.FPX_OFFSET + i)*ct.hxxgrid * ct.xside;
        for(j = 0; j < pct.FPY0_GRID; j++)
        {
            y = (pct.FPY_OFFSET + j)*ct.hyygrid * ct.yside;

            for(k = 0; k < pct.FPZ0_GRID; k++)
            {
                z = (pct.FPZ_OFFSET + k)*ct.hzzgrid * ct.zside;
        
                idx = i * pct.FPY0_GRID * pct.FPZ0_GRID + j*pct.FPZ0_GRID + k;
                dipole_ion_x += x * rhoc[idx];
                dipole_ion_y += y * rhoc[idx];
                dipole_ion_z += z * rhoc[idx];
                dipole_ele_x += x * rho[idx];
                dipole_ele_y += y * rho[idx];
                dipole_ele_z += z * rho[idx];
            }
        }
    }

    dipole_ion_x *= ct.vel_f;
    dipole_ion_y *= ct.vel_f;
    dipole_ion_z *= ct.vel_f;
    dipole_ele_x *= ct.vel_f;
    dipole_ele_y *= ct.vel_f;
    dipole_ele_z *= ct.vel_f;

printf("\n  x dipolll  %f %f", dipole_ion_x, dipole_ele_x);
printf("\n  y dipolll  %f %f", dipole_ion_y, dipole_ele_y);
printf("\n  z dipolll  %f %f", dipole_ion_z, dipole_ele_z);

    my_malloc( Xij, n2, double );
    my_malloc( Yij, n2, double );
    my_malloc( Zij, n2, double );
    get_phi_xyz_phi(states, Xij, Yij, Zij);





    //for(i = 0; i < n2; i++) mat_X[i]= Pn1[i];
    //update_TDDFT(mat_X);
    //  get_dm_diag_p(states, l_s, mat_X, Hij);

    //  write_eigs(states);
    /*  Xij = <phi|x|phi>, Yij = <phi|y|phi>, Zij = <phi|z|phi>  */ 

    double time_step = 0.20;
    double dipole_m,  efield = 0.001;
    int ione = 1;
    double fs= 0.02418884;  // 1fs = 0.02418884 *10^-15 second 
FILE *dfi;
 dfi = fopen("dipole.dat.ykick", "w+");
    for(ct.scf_steps = 0; ct.scf_steps < ct.max_scf_steps; ct.scf_steps++)
    {

        if(ct.scf_steps > 0) 
        {
            efield = 0.0;
        }

        for(i = 0; i < n2; i++) Hij_00[i] = time_step*Hij_00[i] + efield * Yij[i];
        eldyn_(&num_states, Bij_00, Hij_00, Pn0, Pn1, &Ieldyn, &iprint);

        for(i = 0; i < n2; i++) mat_X[i]= Pn1[i];
        for(i = 0; i < 2*n2; i++) Pn0[i]= Pn1[i];
        update_TDDFT(mat_X);

    dipole_ele_x = 0.0;
    dipole_ele_y = 0.0;
    dipole_ele_z = 0.0;
    for(i = 0; i < pct.FPX0_GRID; i++)
    {
        x = (pct.FPX_OFFSET + i)*ct.hxxgrid * ct.xside;
        for(j = 0; j < pct.FPY0_GRID; j++)
        {
            y = (pct.FPY_OFFSET + j)*ct.hyygrid * ct.yside;

            for(k = 0; k < pct.FPZ0_GRID; k++)
            {
                z = (pct.FPZ_OFFSET + k)*ct.hzzgrid * ct.zside;
        
                idx = i * pct.FPY0_GRID * pct.FPZ0_GRID + j*pct.FPZ0_GRID + k;
                dipole_ele_x += x * rho[idx];
                dipole_ele_y += y * rho[idx];
                dipole_ele_z += z * rho[idx];
            }
        }
    }

    dipole_ele_x *= ct.vel_f;
    dipole_ele_y *= ct.vel_f;
    dipole_ele_z *= ct.vel_f;

//        dipole_m = ddot(&n2, Pn1, &ione, Xij, &ione);
    dipole_ele_x -= dipole_ion_x;
    dipole_ele_y -= dipole_ion_y;
    dipole_ele_z -= dipole_ion_z;
 

        fprintf(dfi, "\n  %f  %18.10f  %18.10f  %18.10f ",
ct.scf_steps*time_step, dipole_ele_x, dipole_ele_y, dipole_ele_z);
        // get_dm_diag_p(states, l_s, mat_X, Hij);
        // write_eigs(states);


    }

    //    my_free(Xij);
    //    my_free(Yij);
    //    my_free(Zij);

    MPI_Finalize();

    return 0;
}


