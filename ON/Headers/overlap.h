/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*  Overlap related struct and functions  

			-----Qingzhong 

*/


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


