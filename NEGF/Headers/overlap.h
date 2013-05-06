/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*  Overlap related struct and functions  

			-----Qingzhong 

*/


struct ION_ORBIT_OVERLAP
{
 short int xlow1;
 short int xhigh1;
 short int xlow2;
 short int xhigh2;
 short int xshift;
  
 short int ylow1;
 short int yhigh1;
 short int ylow2;
 short int yhigh2;
 short int yshift;

 short int zlow1;
 short int zhigh1;
 short int zlow2;
 short int zhigh2;
 short int zshift;
 
 short int flag;
};
typedef struct  ION_ORBIT_OVERLAP ION_ORBIT_OVERLAP;


struct ORBIT_ORBIT_OVERLAP
{
 short int xlow1;
 short int xhigh1;
 short int xlow2;
 short int xhigh2;
 short int xshift;
  
 short int ylow1;
 short int yhigh1;
 short int ylow2;
 short int yhigh2;
 short int yshift;

 short int zlow1;
 short int zhigh1;
 short int zlow2;
 short int zhigh2;
 short int zshift;
 
 short int flag;
};
typedef struct  ORBIT_ORBIT_OVERLAP ORBIT_ORBIT_OVERLAP;

ORBIT_ORBIT_OVERLAP 	*orbit_overlap_region; 
ION_ORBIT_OVERLAP 	*ion_orbit_overlap_region_nl; 
ION_ORBIT_OVERLAP 	*ion_orbit_overlap_region_loc; 

REAL *rho_global;
REAL *kbpsi, *kbpsi_comm;
REAL *partial_kbpsi_x, *partial_kbpsi_y, *partial_kbpsi_z; 
int *num_nonlocal_ion;
int *ionidx_allproc;
int max_ion_nonlocal;

int *state_begin;
int *state_end;

int *send_to, *recv_from, num_sendrecv_loop;
int *send_to1, *recv_from1, num_sendrecv_loop1;


