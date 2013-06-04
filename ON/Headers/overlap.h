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
typedef struct ION_ORBIT_OVERLAP ION_ORBIT_OVERLAP;


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
typedef struct ORBIT_ORBIT_OVERLAP ORBIT_ORBIT_OVERLAP;

ORBIT_ORBIT_OVERLAP *orbit_overlap_region;
ION_ORBIT_OVERLAP *ion_orbit_overlap_region_nl;

REAL *rho_global;
REAL *wave_global;

REAL *kbpsi, *kbpsi_comm, *partial_kbpsi_x, *partial_kbpsi_y, *partial_kbpsi_z;
int *num_nonlocal_ion;
int *ionidx_allproc;
int max_ion_nonlocal;

int *state_begin;
int *state_end;


int *send_to, *recv_from, num_sendrecv_loop;
int *send_to1, *recv_from1, num_sendrecv_loop1;


REAL dot_product_orbit_nl (STATE *st1, int ion2, REAL * psi, REAL * prjptr);

void non_zero_pairs ();
void non_zero_pairs1 ();

void init_nl_xyz ();

void theta_phi_new (int st1, int st2, REAL theta_ion, REAL * st2_psi,
                    REAL * state1_psi, int mode, STATE * states);

void print_status (STATE *, REAL *, REAL *, REAL *, REAL *, char *);
void print_state_projections (STATE *, char);
void print_global_function (REAL *, char, char *);
void print_state_sum (STATE * states);
void print_state (STATE * state);
void print_sum (int size, double *data, char *msg);
void print_sum_square (int size, double *data, char *msg);


void init_comm (STATE * states);

void init_comm_res (STATE * states);


void get_orbit_overlap_region (STATE * states);

void get_ion_orbit_overlap_nl (STATE * states);

void duplicate_states_info (STATE * states, STATE * states1);


void get_ion_ion_overlap_region_orbit ();


void is_state_overlap (STATE * states);
