#pragma once
//#include "typedefs.h"
//#define eldyn_ eldyn
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
void eldyn(int *num_states, double *, double *, double *, double *, int *, int*);
FILE *my_fopen_increment(char *name);


#ifdef __cplusplus
}
#endif

template <typename KpointType>
void XyzMatrix (Kpoint<KpointType> *kptr, KpointType *Aij, int n, int m, int l);

template <typename OrbitalType, typename MatrixType> 
    void RmgTddft ( spinobj<double> &vxc,
                   fgobj<double> &vh,
                   fgobj<double> &vnuc,
                   spinobj<double> &rho,
                   fgobj<double> &rhocore,
                   fgobj<double> &rhoc,Kpoint<OrbitalType> **Kptr);

template <typename OrbitalType> 
    void VecPHmatrix (Kpoint<OrbitalType> *kptr, double *efield_tddft, int *desca, int tddft_start_state, int num_states);
template <typename OrbitalType> void CurrentNlpp (Kpoint<OrbitalType> *kptr, int *desca, int tddft_start_state, int num_states);
void CurrentOperator (Kpoint<std::complex<double>> *kptr, int *desca, int tddft_start_state);
void CurrentOperator (Kpoint<double> *kptr, int *desca, int tddft_start_state);

template <typename KpointType, typename CalType, typename MatrixType>
void HmatrixUpdate (Kpoint<KpointType> *kptr, wfobj<double> vtot_psi, wf_spinobj<double> vxc_psi, MatrixType *Aij, int tddft_start_state, int num_states, int *desca);
template <typename KpointType>
void HSmatrix (Kpoint<KpointType> *kptr, double *vtot_eig, double *vxc_psi,  KpointType *Aij, KpointType *Sij);
void ReadData_rmgtddft (const char *filename, double * vh, double * vxc, 
        double *vh_corr, double *Pn0, double *Hmatrix, double *H0, double *H1,int *tot_steps, int n2, int n2_C,int numst);
void WriteData_rmgtddft (const char *filename, double * vh, double * vxc, 
        double *vh_corr, double *Pn0, double *Hmatrix, double *H0, double *H1,int tot_steps, int n2,int n2_C,int numst);
void ReadData_rmgtddft_on (char *filename, double * vh, double * vxc, 
        double *vh_corr, double *Pn0, double *Hmatrix, double *Smatrix, double *Cmat, double *H0, double *H1,int *tot_steps, int n2);
void WriteData_rmgtddft_on (char *filename, double * vh, double * vxc, 
        double *vh_corr, double *Pn0, double *Hmatrix, double *Smatrix, double *Cmat, double *H0, double *H1,int tot_steps, int n2);

template <typename OrbitalType, typename CalType, typename MatrixType> 
void GetNewRho_rmgtddft (Kpoint<OrbitalType> *kptr, spinobj<double> &rho, MatrixType *rho_matrix, int numst, int tddft_start_state, Scalapack &);

template <typename MatrixType>
    void MatDiagSet (MatrixType *mat,  std::vector<double> diag_elem, double beta, int numst, Scalapack &SP);

void  magnus( double *H0, double *H1, double p_time_step , double *Hdt, int ldim);
void tst_conv_matrix  (double * p_err , int * p_ij_err ,   double *H0, double *H1,  int ldim, MPI_Comm comm);
void extrapolate_Hmatrix   (double  *Hm1, double *H0, double *H1,  int ldim);
template <typename OrbitalType> void TddftEnergyInit ( spinobj<double> &vxc,
        fgobj<double> &vh, fgobj<double> &vnuc,
        spinobj<double> &rho_ground,
        fgobj<double> &rhocore, fgobj<double> &rhoc,
        Kpoint<OrbitalType> **Kptr,
        Scalapack &SP, int Mdim, int Ndim, std::vector<double> &Eterms);
template <typename OrbitalType, typename MatrixType> void TddftEnergy (fgobj<double> &vh, spinobj<double> &rho, fgobj<double> &rhoc,
        Kpoint<OrbitalType> **Kptr, int Mdim, int Ndim, std::vector<double> &Eterms, MPI_Comm mat_comm);


