
int MXLLDA, MXLCOL;
rmg_double_t *rho, *rho_old, *rhoc, *vh, *vnuc, *vxc, *vext,  *vcomp, *rhocore, *vtot, *vtot_old, *peaks;
rmg_double_t *vh_old, *vxc_old; 
rmg_double_t *statearray, *l_s, *matB, *mat_hb, *mat_X, *Hij, *theta, *work_dis;
double *mat_local;
rmg_double_t *work_dis2, *zz_dis, *gamma_dis, *uu_dis; 
rmg_double_t *work_matrix, *nlarray1;
rmg_double_t *projectors, *projectors_x, *projectors_y, *projectors_z;
rmg_double_t *sg_twovpsi, *sg_res;
int *nlindex;
rmg_double_t *work_memory;
rmg_double_t *sg_orbit;
rmg_double_t *sg_orbit_res;
rmg_double_t *orbit_tem;
rmg_double_t *vtot_global, *vtot_c;
double complex *sigma_all;


int NPES;
int peakNum;
int state_to_ion[MAX_STATES];
int state_to_proc[MAX_STATES];
char num_loc_st[MAX_IONS];
short int state_overlap_or_not[MAX_STATES * MAX_STATES];


STATE states[MAX_STATES];
STATE states1[MAX_STATES];
STATE states_distribute[MAX_STATES];


ORBIT_ORBIT_OVERLAP 	*orbit_overlap_region; 
ION_ORBIT_OVERLAP 	*ion_orbit_overlap_region_nl; 
ION_ORBIT_OVERLAP 	*ion_orbit_overlap_region_loc; 

rmg_double_t *rho_global;
rmg_double_t *kbpsi, *kbpsi_comm;
rmg_double_t *partial_kbpsi_x, *partial_kbpsi_y, *partial_kbpsi_z; 
int *num_nonlocal_ion;
int *ionidx_allproc;
int max_ion_nonlocal;

int *state_begin;
int *state_end;

int *send_to, *recv_from, num_sendrecv_loop;
int *send_to1, *recv_from1, num_sendrecv_loop1;


short int vloc_state_overlap_or_not[MAX_IONS * MAX_STATES];

rmg_double_t *vnuc_x, *vnuc_y, *vnuc_z;
