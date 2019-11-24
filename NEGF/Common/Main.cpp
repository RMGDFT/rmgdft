#include "negf_prototypes.h"
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


#include "LCR.h"
#include "prototypes_on.h"
#include "prototypes_negf.h"
#include "transition.h"


#include "../Headers/common_prototypes.h"
#include "../Headers/common_prototypes1.h"



#include "twoParts.h"
#include "pmo.h"
#include "cei.h"


#include "init_var.h"
#include "Kbpsi.h"

std::vector<ION> Atoms;
std::vector<SPECIES> Species;

/* Main control structure which is declared extern in main.h so any module */
/* may access it.					                 */
CONTROL ct;

parallel_matrix_operation pmo;

/* PE control structure which is also declared extern in main.h */
PE_CONTROL pct;

int mpi_nprocs;
int mpi_myrank;

int total_mem = 0;

/*Variables from recips.h*/
double b0[3], b1[3], b2[3];
double alat;

std::unordered_map<std::string, InputKey *> ControlMap;

COMPASS potentialCompass, chargeDensityCompass;


unsigned int *perm_ion_index;
double *projectors, *projectors_x, *projectors_y, *projectors_z;
int *num_nonlocal_ion;
double *kbpsi, *kbpsi_comm, *partial_kbpsi_x, *partial_kbpsi_y, *partial_kbpsi_z;
double *kbpsi_res;
int *kbpsi_comm_send, *kbpsi_comm_recv,  kbpsi_num_loop;
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
ION_ORBIT_OVERLAP *ion_orbit_overlap_region_nl;
double *rho, *rho_old, *rhoc, *vh, *vnuc, *vxc, *rhocore, *eig_rho, *vtot, *vtot_c, *rho_tf;
double *rho_oppo, *rho_tot;
int MXLLDA, MXLCOL;
double *sg_twovpsi, *sg_res;
double *statearray, *l_s, *matB, *mat_hb, *mat_X, *Hij, *theta, *work_dis;
double *Hij_00, *Bij_00;
double *work_matrix_row, *coefficient_matrix_row, *nlarray1;
double *work_dis2, *zz_dis, *cc_dis, *gamma_dis, *uu_dis, *mat_Omega;
ORBIT_ORBIT_OVERLAP *orbit_overlap_region;
std::vector<ORBITAL_PAIR> OrbitalPairs;
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
double *vh_old, *vxc_old;

double *vcomp, *peaks, *vext ;
ION_ORBIT_OVERLAP    *ion_orbit_overlap_region_loc;

double *work_matrix;
double *vnuc_x, *vnuc_y, *vnuc_z;
int peakNum;

DoubleC *sigma_all;

NON_LINEAR_THREE_PART lcr[NUM_SUBSYSTEM_MAX];

complex_energy_integral cei;

KBPSI Kbpsi_str;


MPI_Comm COMM_PEX, COMM_PEY, COMM_PEZ, COMM_3D;
MPI_Comm COMM_EN, COMM_EN1, COMM_EN2;



void ReadBranchNEGF(char *cfile, CONTROL& lc, complex_energy_integral& cei, COMPASS& potcompass, COMPASS& rhocompass);

int main (int argc, char **argv)
{


    RmgTimer *RT = new RmgTimer("1-TOTAL");

    // Set branch type
    ct.rmg_branch = RMG_NEGF;
    ct.save_args(argc, argv);


    /* Define a default output stream, gets redefined to log file later */
    ct.logfile = stdout;


    ct.images_per_node = 1;
    ct.proj_nophase = 1;
    InitIo(argc, argv, ControlMap);

    ReadBranchNEGF(ct.cfile, ct, cei, potentialCompass, chargeDensityCompass);
    states = new STATE[ct.num_states];
    allocate_states();

    perm_ion_index = (unsigned int *) malloc(ct.num_ions * sizeof(int));
    for(int i = 0; i < ct.num_ions; i++) perm_ion_index[i] = i;
    ReadOrbitals (ct.cfile, states, Atoms, pct.img_comm, perm_ion_index);
    get_state_to_proc(states);

    MPI_Barrier(pct.img_comm);

    /*  Begin to do the real calculations */
    Run (states, states1, ControlMap);


    delete(RT);

    if(pct.imgpe == 0) fclose(ct.logfile);
    RmgPrintTimings(pct.img_comm, ct.logname, ct.scf_steps, pct.num_owned_ions * ct.num_kpts_pe);


    MPI_Finalize ();


    return 0;
}                               /*   end main */

