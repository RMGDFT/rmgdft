//#include "typedefs.h"
void GetNewRhoLocal (STATE * states_distribute, double *rho, double *mat_local, double *rho_matrix);
void InitStatedistribute ();
void MapOrbitalToProcess(int st2, STATE *states, STATE *states_distribute, double *psi_whole);
void MatrixToLocal (STATE *states_distribute, double * A_glob, double * A_local);


#ifdef __cplusplus
extern "C" {
#endif

void read_rhomatrix(char *, double *);
void init_TDDFT(STATE *states, STATE *states1);
void get_phi_xyz_phi(STATE *states, double *, double *, double *);
void dipole_calculation(double *, double *);
void update_TDDFT(double *);
void Cpdgemr2d(int m, int n,
                double* a, int ia, int ja, int* desca,
                double* b, int ib, int jb, int* descb,
                int gcontext);

void mat_dist_to_global (double *dist_mat, double *global_mat, int *desca);
void mat_global_to_dist (double *global_mat, double *dist_mat, int *desca);

void eldyn_(int *num_states, double *, double *, double *, double *, int *, int*);
FILE *my_fopen_increment(char *name);
void write_data(char *name, double *vh, double *vxc, double *vh_old,
        double *vxc_old, double *rho, STATE * states);


#ifdef __cplusplus
}
#endif
