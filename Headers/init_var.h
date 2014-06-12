/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#define ORDER_N 1



int MXLLDA, MXLCOL;
rmg_double_t *rho, *rho_old, *rhoc, *vh, *vnuc, *vxc, *rhocore, *eig_rho, *vtot,
    *vtot_c;
double *rho_oppo, *rho_tot;
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
char *state_overlap_or_not;



int *state_to_ion;
int *state_to_proc;
STATE *states;
STATE *states1;
STATE *states_tem;
STATE *states_distribute;
char *vloc_state_overlap_or_not;



ORBIT_ORBIT_OVERLAP *orbit_overlap_region;
ION_ORBIT_OVERLAP *ion_orbit_overlap_region_nl;

rmg_double_t *rho_global;
rmg_double_t *wave_global;

rmg_double_t *kbpsi, *kbpsi_comm, *partial_kbpsi_x, *partial_kbpsi_y, *partial_kbpsi_z;
int *num_nonlocal_ion;
int *ionidx_allproc;
int max_ion_nonlocal;

int *state_begin;
int *state_end;


int *send_to, *recv_from, num_sendrecv_loop;
int *send_to1, *recv_from1, num_sendrecv_loop1;

rmg_double_t *vcomp, *peaks, *vext ;

double *mat_local;
double complex *sigma_all;
int peakNum;
rmg_double_t *work_matrix;
ION_ORBIT_OVERLAP 	*ion_orbit_overlap_region_loc; 
rmg_double_t *vnuc_x, *vnuc_y, *vnuc_z;



