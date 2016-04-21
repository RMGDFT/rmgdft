//#include "typedefs.h"
void GetNewRhoLocal (STATE * states_distribute, double *rho, double *mat_local, double *rho_matrix);
void InitStatedistribute (STATE *states_distrubte);

void MapOrbitalToProcess(int st2, STATE *states, STATE *states_distribute, double *psi_whole);
void MatrixToLocal (STATE *states_distribute, double * A_glob, double * A_local);
void MatrixToGlobal (STATE *states_distribute, double * A_local, double * A_glob);

void HijUpdateNCpp (STATE * states_distribute, double *vtot_c, double *mat_local, double *mat_global);


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


#ifdef __cplusplus
}
#endif


template <typename OrbitalType> 
void RmgTddft (double * vxc, double *, double * vnuc, double * rho, double * rho_oppo,
     double * rhocore, double * rhoc, Kpoint<OrbitalType> **Kptr);
template <typename KpointType>
void HmatrixUpdate (Kpoint<KpointType> *kptr, double *vtot_eig, KpointType *Aij);
void VhDriver(double *, double *, double *);
template <typename KpointType>
void HSmatrix (Kpoint<KpointType> *kptr, double *vtot_eig, KpointType *Aij, KpointType *Sij);
void GetNewRho_rmgtddft(double *, double *, double *, double *, int);
void ReadData_rmgtddft (char *filename, double * vh, double * vxc, 
        double *vh_corr, double *Pn0, double *Hmatrix, double *Smatrix, int *tot_steps);
void WriteData_rmgtddft (char *filename, double * vh, double * vxc, 
        double *vh_corr, double *Pn0, double *Hmatrix, double *Smatrix, int tot_steps);

