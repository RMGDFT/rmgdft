/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#define ORDER_N 1

extern double *projectors, *projectors_x, *projectors_y, *projectors_z;
extern int *num_nonlocal_ion;
extern double *kbpsi, *kbpsi_comm, *partial_kbpsi_x, *partial_kbpsi_y, *partial_kbpsi_z;
extern char *state_overlap_or_not;
extern int *send_to, *recv_from, num_sendrecv_loop;
extern int *send_to1, *recv_from1, num_sendrecv_loop1;
extern int *ionidx_allproc;
extern int max_ion_nonlocal;
extern int NPES;
extern STATE *states_tem;
extern int *state_to_ion;
extern int *state_to_proc;
extern STATE *states;
extern STATE *states1;
extern STATE *states_distribute;
extern ION_ORBIT_OVERLAP *ion_orbit_overlap_region_nl;
extern double *rho, *rho_old, *rhoc, *vh, *vnuc, *vxc, *rhocore, *eig_rho, *vtot, *vtot_c;
extern double *rho_oppo, *rho_tot;
extern int MXLLDA, MXLCOL;
extern double *sg_twovpsi, *sg_res;
extern double *statearray, *l_s, *matB, *mat_hb, *mat_X, *Hij, *theta, *work_dis;
extern double *Hij_00, *Bij_00;
extern double *work_matrix_row, *coefficient_matrix_row, *nlarray1;
extern double *work_dis2, *zz_dis, *cc_dis, *gamma_dis, *uu_dis, *mat_Omega;
extern ORBIT_ORBIT_OVERLAP *orbit_overlap_region;
extern char *vloc_state_overlap_or_not;
extern double *orbit_tem;
extern double *sg_orbit;
extern double *sg_orbit_res;
extern int *state_begin;
extern int *state_end;
extern double *vtot_global;
extern double *work_memory;
extern double *wave_global;
extern double *rho_global;
extern double *vh_old, *vxc_old;

extern double *mat_local;

extern double *vcomp, *peaks, *vext ;
extern ION_ORBIT_OVERLAP 	*ion_orbit_overlap_region_loc; 
extern double *work_matrix;
extern double *vnuc_x, *vnuc_y, *vnuc_z;
extern int peakNum;
#if 0


int *nlindex;


int NX_GRID, NY_GRID, NZ_GRID;
int P0_BASIS;
int S0_BASIS;
int PX0_GRID;
int PY0_GRID;
int PZ0_GRID;













//double complex *sigma_all;
#endif


