/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/md.h *****
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
 *   
 * INPUTS
 *
 * OUTPUT
 *  
 * PARENTS
 *
 * CHILDREN
 * 
 * SEE ALSO
 *  
 * SOURCE
 */


#define ORDER_N 1


/* Version information */
#include "main.h"


/* Compile time parameters */
#include    "params_on.h"


/* Constants and symbolic definitions */
#include    "const_on.h"


/** Size of floating point variables used in QMD */
#include    "common_prototypes.h"


int MXLLDA, MXLCOL;
rmg_double_t *rho, *rho_old, *rhoc, *vh, *vnuc, *vxc, *rhocore, *eig_rho, *vtot,
    *vtot_c;
rmg_double_t *vh_old, *vxc_old;
rmg_double_t *statearray, *l_s, *matB, *mat_hb, *mat_X, *Hij, *theta, *work_dis;
double *Hij_00, *Bij_00;
rmg_double_t *work_dis2, *zz_dis, *cc_dis, *gamma_dis, *uu_dis, *mat_Omega;
rmg_double_t *work_matrix_row, *coefficient_matrix_row, *nlarray1;
rmg_double_t *projectors, *projectors_x, *projectors_y, *projectors_z;
rmg_double_t *sg_twovpsi, *sg_res;
int *nlindex;
rmg_double_t *work_memory;
rmg_double_t *sg_orbit;
rmg_double_t *sg_orbit_res;
rmg_double_t *orbit_tem;
rmg_double_t *vtot_global;


int NPES;
int NX_GRID, NY_GRID, NZ_GRID;
int P0_BASIS;
int S0_BASIS;
int PX0_GRID;
int PY0_GRID;
int PZ0_GRID;
int state_to_ion[MAX_STATES];
int state_to_proc[MAX_STATES];
char num_loc_st[MAX_IONS];
short int overlap_or_not[MAX_IONS * MAX_IONS];
short int state_overlap_or_not[MAX_STATES * MAX_STATES];


/* Array telling for each pair of ions if the localization regions
   around touch both a same region (anywhere on any PE) */
int array_both_not_zero[MAX_IONS][MAX_IONS];


STATE states[MAX_STATES];
STATE states1[MAX_STATES];
STATE states_res[MAX_STATES];
STATE states_res1[MAX_STATES];
STATE states_tem[MAX_STATES];


#include "overlap.h"

#include "macros.h"
#include    "prototypes_on.h"
