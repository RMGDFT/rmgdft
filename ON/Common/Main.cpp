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
#include "svnrev.h"
#include <float.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
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
//#include "main.h"
#include "prototypes_on.h"
#include "init_var.h"
#include "transition.h"


#include "../Headers/common_prototypes.h"
#include "../Headers/common_prototypes1.h"




/* Main control structure which is declared extern in main.h so any module */
/* may access it.					                 */
CONTROL ct;

/* PE control structure which is also declared extern in main.h */
PE_CONTROL pct;

double *projectors, *projectors_x, *projectors_y, *projectors_z;
int *num_nonlocal_ion;
double *kbpsi, *kbpsi_comm, *kbpsi_res, *partial_kbpsi_x, *partial_kbpsi_y, *partial_kbpsi_z;
int kbpsi_num_loop, *kbpsi_comm_send, *kbpsi_comm_recv;
char *state_overlap_or_not;
int *send_to, *recv_from, num_sendrecv_loop;
int *send_to1, *recv_from1, num_sendrecv_loop1;
int *ionidx_allproc;
int max_ion_nonlocal;
int NPES;
STATE *states_tem;
int *state_to_ion;
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

//#include "El.hpp"
//using namespace std;
//using namespace El;

int main(int argc, char **argv)
{




   // double a[4], b[4];
   // MPI_Init_thread(&argc, &argv, ct.mpi_threadlevel, &provided);
    //MPI_Init(&argc, &argv);

    //DiagElemental(4, a, b);    


    ct.mpi_threadlevel = MPI_THREAD_SERIALIZED;
    ct.mpi_threadlevel = 0;
   // mpi::Initialize(argc, argv);
    RmgTimer *RT = new RmgTimer("1-TOTAL");
    
    /* Define a default output stream, gets redefined to log file later */
    ct.logfile = stdout;
   // tem();

    ct.images_per_node = 1;
    InitIo(argc, argv, ControlMap);

    //  initialize for ELEMENTAl lib
    //Initialize( argc, argv );

    ReadBranchON(ct.cfile, ct, ControlMap);
    allocate_states();
    get_state_to_proc(states);
    ReadOrbitals (ct.cfile, states, state_to_ion, pct.img_comm);

    init_states();
    my_barrier();

    /*  Begin to do the real calculations */
    run(states, states1);


    delete(RT);


    if(pct.imgpe == 0) fclose(ct.logfile);
    RmgPrintTimings(Rmg_G, ct.logname, ct.scf_steps);


    MPI_Finalize();

    return 0;
}


