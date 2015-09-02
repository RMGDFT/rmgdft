/************************** SVN Revision Information **************************
 **    $Id: main.c 2671 2014-10-13 19:50:57Z luw $    **
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




#include <sys/stat.h>
#include <sys/types.h>
#include <float.h>
#include <fcntl.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include "svnrev.h"
#include <unordered_map>
#include "const.h"
#include "RmgTimer.h"
#include "RmgException.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "Kpoint.h"
#include "InputKey.h"
#include "blas.h"
#include "prototypes_on.h"
#include "prototypes_tddft.h"
#include "init_var.h"
#include "transition.h"

#include "Kbpsi.h"

#include "../Headers/common_prototypes.h"
#include "../Headers/common_prototypes1.h"



/* Main control structure which is declared extern in main.h so any module */
/* may access it.					                 */
CONTROL ct;

/* PE control structure which is also declared extern in main.h */
PE_CONTROL pct;

KBPSI Kbpsi_str;
unsigned int *perm_ion_index, *perm_state_index, *rev_perm_state_index;

double *projectors, *projectors_x, *projectors_y, *projectors_z;
int *num_nonlocal_ion;
double *kbpsi, *kbpsi_comm, *partial_kbpsi_x, *partial_kbpsi_y, *partial_kbpsi_z;
double *kbpsi_res; 
int kbpsi_num_loop, *kbpsi_comm_send, *kbpsi_comm_recv;
char *state_overlap_or_not;
int *send_to, *recv_from, num_sendrecv_loop;
int *send_to1, *recv_from1, num_sendrecv_loop1;
int *ionidx_allproc;
int max_ion_nonlocal;
int NPES;
STATE *states_tem;
int *state_to_proc;
STATE *states;
STATE *states1;
STATE *states_distribute;
ION_ORBIT_OVERLAP *ion_orbit_overlap_region_nl;
double *rho, *rho_old, *rhoc, *vh, *vnuc, *vxc, *rhocore, *eig_rho, *vtot, *vtot_c;
double *rho_oppo, *rho_tot;
int MXLLDA, MXLCOL;
double *sg_twovpsi, *sg_res;
double *statearray, *l_s, *matB, *mat_hb, *mat_X, *Hij, *theta, *work_dis;
double *Hij_00, *Bij_00;
double *work_matrix_row, *coefficient_matrix_row, *nlarray1;
double *work_dis2, *zz_dis, *cc_dis, *gamma_dis, *uu_dis, *mat_Omega;
ORBIT_ORBIT_OVERLAP *orbit_overlap_region;
char *vloc_state_overlap_or_not;
double *orbit_tem;
double *sg_orbit;
double *sg_orbit_res;
int *state_begin;
int *state_end;
double *vtot_global;
double *work_memory;
double *wave_global;
double *rho_global;
double *vxc_old, *vh_old;


int mpi_nprocs;
int mpi_myrank;

//STATE *states, *states1;

/*Variables from recips.h*/
double b0[3], b1[3], b2[3];
double alat;

std::unordered_map<std::string, InputKey *> ControlMap;


int main(int argc, char **argv)
{



    double *Hmatrix, *Smatrix, *Xij_00, *Yij_00, *Zij_00;
    double *Xij_dist, *Yij_dist, *Zij_dist;
    double *rho_matrix;
    
    int Ieldyn=1, iprint = 0;
    int n2, numst,i;
    double *Pn0, *Pn1;
    double dipole_ion[3], dipole_ele[3];
    int IA=1, JA=1, IB=1, JB=1;
    int ione=1;
    FILE *dfi;
    int tot_steps, pre_steps;
    int amode, fhand;
    int FP0_BASIS;

    char filename[MAX_PATH+200];
    char char_tem[MAX_PATH+200];
    /* Define a default output stream, gets redefined to log file later */
    ct.logfile = stdout;
    amode = S_IREAD | S_IWRITE;


    FP0_BASIS = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);
    ct.images_per_node = 1;
    InitIo(argc, argv, ControlMap);


    ReadBranchON(ct.cfile, ct, ControlMap);
    allocate_states();
    get_state_to_proc(states);
    perm_ion_index = (unsigned int *) malloc(ct.num_ions * sizeof(int));
    for(int i = 0; i < ct.num_ions; i++) perm_ion_index[i] = i;

    perm_state_index = (unsigned int *) malloc(ct.num_states * sizeof(int));
    rev_perm_state_index = (unsigned int *) malloc(ct.num_states * sizeof(int));

    if(pct.gridpe == 0) 
    {
        ReadPermInfo(ct.infile, perm_ion_index);
    }
    MPI_Bcast(perm_ion_index, ct.num_ions, MPI_INT, 0, pct.grid_comm);
    ReadOrbitals (ct.cfile, states, ct.ions, pct.img_comm, perm_ion_index);

    init_states();
    my_barrier();




    RmgTimer *RT = new RmgTimer("Main");
    RmgTimer *RT1 = new RmgTimer("Main: init");



    init_dimension(&MXLLDA, &MXLCOL);
    init_pe_on();


    /* allocate memory for matrixs  */
    allocate_matrix();

    /* Perform some necessary initializations no matter localized or not  
     */
    vxc_old = new double[Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO)];
    vh_old = new double[Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO)];


    InitON(vh, rho, rho_oppo, rhocore, rhoc, states, states1, vnuc, vxc, vh_old, vxc_old);

    my_barrier();


    numst = ct.num_states;

    get_HS(states, states1, vtot_c, Hij_00, Bij_00);

    Cpdgemr2d(numst, numst, Hij_00, IA, JA, pct.descb, Hij, IB, JB,
            pct.desca, pct.desca[1]);
    Cpdgemr2d(numst, numst, Bij_00, IA, JA, pct.descb, matB, IB, JB,
            pct.desca, pct.desca[1]);


    get_dm_diag_p(states, matB, mat_X, Hij);

    write_eigs(states);

    n2 = (ct.state_end - ct.state_begin) * ct.num_states;
    Xij_00 = new double[n2];
    Yij_00 = new double[n2];
    Zij_00 = new double[n2];
    rho_matrix = new double[n2];


    n2 = MXLLDA * MXLCOL;
    Xij_dist = new double[n2];
    Yij_dist = new double[n2];
    Zij_dist = new double[n2];

    /*  Xij = <phi|x|phi>, Yij = <phi|y|phi>, Zij = <phi|z|phi>  */ 

    get_phi_xyz_phi(states, Xij_00, Yij_00, Zij_00);

    Cpdgemr2d(numst, numst, Xij_00, IA, JA, pct.descb, Xij_dist, IB, JB,
            pct.desca, pct.desca[1]);
    Cpdgemr2d(numst, numst, Yij_00, IA, JA, pct.descb, Yij_dist, IB, JB,
            pct.desca, pct.desca[1]);
    Cpdgemr2d(numst, numst, Zij_00, IA, JA, pct.descb, Zij_dist, IB, JB,
            pct.desca, pct.desca[1]);

    n2 = numst * numst;
    Pn0 = new double[2*n2];
    Pn1 = new double[2*n2];
    Hmatrix = new double[n2];
    Smatrix = new double[n2];

    /*  matB: overlap matrix, Hij:  Hamiltonian matrix  distributed in
     *  scalapack way*/

    dipole_calculation(rhoc, dipole_ion);
    dipole_calculation(rho, dipole_ele);


    rmg_printf("\n  x dipolll  %f %f", dipole_ion[0], dipole_ele[0]);
    rmg_printf("\n  y dipolll  %f %f", dipole_ion[1], dipole_ele[1]);
    rmg_printf("\n  z dipolll  %f %f", dipole_ion[2], dipole_ele[2]);



    mat_dist_to_global(mat_X, Pn0, pct.desca);
    mat_dist_to_global(matB, Smatrix, pct.desca);

    for(i = 0; i < n2; i++) Pn0[i+n2] = 0.0;

    delete(RT1);


    double time_step = 0.20;
    double efield[3];
    efield[0] = ct.x_field_0 * ct.e_field;
    efield[1] = ct.y_field_0 * ct.e_field;
    efield[2] = ct.z_field_0 * ct.e_field;


    double fs= 0.02418884;  // 1fs = 0.02418884 *10^-15 second 
    if(pct.gridpe == 0)
    {
        sprintf(char_tem, "%s%s", pct.image_path[pct.thisimg], "dipole.dat");
        int name_incr = FilenameIncrement(char_tem);
        sprintf(filename, "%s.%02d", char_tem, name_incr);

        dfi = fopen(filename, "w");
    }

    if(pct.gridpe == 0)fprintf(dfi, "\n  electric field:  %f  %f  %f ",efield[0], efield[1], efield[2]);

    pre_steps = 0;

    if(ct.runflag == Restart_TDDFT) 
    {
        sprintf(filename, "%s%s", ct.infile, ".TDDFT_restart");
        fhand = open(filename, O_RDWR);
        if (fhand < 0)
            rmg_error_handler(__FILE__, __LINE__, " Unable to write file ");
        read(fhand, &pre_steps,  sizeof(int));
        read(fhand, Pn0, 2* numst * numst*sizeof(double));
        close(fhand);
    }


    for(ct.scf_steps = 0; ct.scf_steps < ct.max_scf_steps; ct.scf_steps++)
    {
        tot_steps = ct.scf_steps + pre_steps;

        if(tot_steps == 0 ) 
        {
            for(i = 0; i < MXLLDA * MXLCOL; i++) Hij[i] = time_step*Hij[i] 
                + efield[0] * Xij_dist[i] + efield[1] * Yij_dist[i] + efield[2] * Zij_dist[i];
        }
        else
        {

            for(i = 0; i < MXLLDA * MXLCOL; i++) Hij[i] = time_step*Hij[i];

        }


        mat_dist_to_global(Hij, Hmatrix, pct.desca);

        RmgTimer *RT2 = new RmgTimer("Main: ELDYN");
        eldyn_(&numst, Smatrix, Hmatrix, Pn0, Pn1, &Ieldyn, &iprint);
        delete(RT2);

        //for(i = 0; i < n2; i++) mat_X[i]= Pn1[i];
        mat_global_to_dist(Pn1, mat_X, pct.desca);
        for(i = 0; i < 2*n2; i++) Pn0[i]= Pn1[i];

        Cpdgemr2d(numst, numst, mat_X, IA, JA, pct.desca, rho_matrix, IB, JB,
                pct.descb, pct.desca[1]);

        dcopy(&FP0_BASIS, rho, &ione, rho_old, &ione);

        GetNewRho_on(states, rho, rho_matrix);

        UpdatePot(vxc, vh, vxc_old, vh_old, vnuc, rho, rho_oppo, rhoc, rhocore);

        RmgTimer *RT4 = new RmgTimer("Main: get_HS");
        get_HS(states, states1, vtot_c, Hij_00, Bij_00);
        delete(RT4);

        Cpdgemr2d(numst, numst, Hij_00, IA, JA, pct.descb, Hij, IB, JB,
                pct.desca, pct.desca[1]);
        dipole_calculation(rho, dipole_ele);

        dipole_ele[0] -= dipole_ion[0];
        dipole_ele[1] -= dipole_ion[1];
        dipole_ele[2] -= dipole_ion[2];


        if(pct.gridpe == 0)fprintf(dfi, "\n  %f  %18.10f  %18.10f  %18.10f ",
                tot_steps*time_step, dipole_ele[0], dipole_ele[1], dipole_ele[2]);
        if((ct.scf_steps +1)%ct.checkpoint == 0)
        {
            write_data(ct.outfile, vh, vxc, vh_old, vxc_old, rho, &states[0]); 
            if(pct.gridpe == 0)
            {
                sprintf(filename, "%s%s", ct.outfile, ".TDDFT_restart");
                fhand = open(filename, O_CREAT |O_TRUNC |O_RDWR, amode);
                if (fhand < 0)
                    rmg_error_handler(__FILE__, __LINE__, " Unable to write file ");

                write(fhand, &tot_steps,  sizeof(int));
                write(fhand, Pn0, 2* numst * numst*sizeof(double));
                close(fhand);
            }

        }


    }


    delete(RT);

    tot_steps = ct.scf_steps + pre_steps;
    write_data(ct.outfile, vh, vxc, vh_old, vxc_old, rho, &states[0]); 
    if(pct.gridpe == 0)
    {
        sprintf(filename, "%s%s", ct.outfile, ".TDDFT_restart");
        fhand = open(filename, O_CREAT |O_TRUNC |O_RDWR, amode);
        if (fhand < 0)
            rmg_error_handler(__FILE__, __LINE__, " Unable to write file ");

        write(fhand, &tot_steps,  sizeof(int));
        write(fhand, Pn0, 2* numst * numst*sizeof(double));
        close(fhand);
    }

    get_dm_diag_p(states, matB, mat_X, Hij);

    write_eigs(states);

    if(pct.imgpe == 0) fclose(ct.logfile);
    RmgPrintTimings(Rmg_G, ct.logname, ct.scf_steps);

    MPI_Finalize();

    return 0;
}


