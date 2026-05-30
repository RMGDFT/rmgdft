#pragma once
#ifndef RMG_rmg_tddft_h
#define RMG_rmg_tddft_h

#include "GridObject.h"
#include "transition.h"
#include "const.h"
#include "State.h"
#include "Kpoint.h"
#include "TradeImages.h"
#include "RmgTimer.h"
#include "rmg_gemm.h"


namespace rmg
{
    template <typename OrbitalType, typename MatrixType> class tddft {
    public:
        tddft(spinobj<double> &vxc_in,
              fgobj<double> &vh_in,
              fgobj<double> &vnuc_in,
              spinobj<double> &rho_in,
              fgobj<double> &rhocore_in,
              fgobj<double> &rhoc_in,
              Kpoint<OrbitalType> **Kptr_in);

        void tddft_md(void);

        void tddft_energy_init (
             spinobj<double> &vxc,
             fgobj<double> &vh,
             fgobj<double> &vnuc,
             spinobj<double> &rho_ground,
             fgobj<double> &rhocore,
             fgobj<double> &rhoc,
             Kpoint<OrbitalType> **Kptr,
             Scalapack &SP,
             int Mdim,
             int Ndim,
             std::vector<double> &Eterms);

        void tddft_energy (
             fgobj<double> &vh,
             spinobj<double> &rho,
             fgobj<double> &rhoc,
             Kpoint<OrbitalType> **Kptr,
             int Mdim,
             int Ndim,
             std::vector<double> &Eterms,
             MPI_Comm mat_comm);

        void  tstconv(
              double *C,
              int *p_M,
              double *p_thrs,
              int *p_ierr,
              double *p_err,
              bool *p_tconv,
              MPI_Comm comm);

        void tstconv(
             float *C,
             int *p_M,
             double *p_thrs,
             int *p_ierr,
             double *p_err,
             bool *p_tconv,
             MPI_Comm comm);

        void gather_rho_matrix(OrbitalType *rho_matrix_global, MatrixType *rho_matrix);

        ~tddft(void);

    private:
        int Mdim, Ndim;
        spinobj<double> &vxc;
        fgobj<double> &vh;
        fgobj<double> &vnuc;
        spinobj<double> &rho;
        fgobj<double> &rhocore;
        fgobj<double> &rhoc;
        fgobj<double> vh_old, vh_dipole, vh_dipole_old;
        spinobj<double> vxc_old, rho_k, rho_ksum, vxc_diff;
        wfobj<double> vtot_psi;
        wf_spinobj<double> vxc_psi;
        fgobj<double> vtot;
        spinobj<double> rho_ground;
        double time_step = ct.tddft_time_step;

        Kpoint<OrbitalType> **Kptr;
        double dipole_tot[3];
        FILE *dfi = NULL, *efi = NULL, *current_fi = NULL, *dbp_fi = NULL;
        std::string filename;
        int n2, n22, n2_C, numst, ione =1;
        int tot_steps = 0;
        int pre_steps = 0;
        int Ieldyn = 1;    // BCH  
                        //int Ieldyn = 2;    // Diagev
        int iprint;
        MPI_Comm eldyn_comm;
        int FP0_BASIS;
        int pbasis;
        int pbasis_noncoll;
        int scalapack_groups = 1;
        Scalapack *Sp = NULL;
        int *desca;
        size_t matrix_size;
        double   thrs_dHmat =1.0e-5 ;
        std::vector<double> Eterms_ground = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        std::vector<double> Eterms_1step = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        std::vector<double> Eterms = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        double efactor;
        const char *eunits;
        double tot_bp_pol;
        double current[3], current0[3];
                
        OrbitalType *psi_dev_pool;
        MatrixType *Pn1         = NULL;
        MatrixType *Hmatrix_1   = NULL;
        MatrixType *Hmatrix     = NULL;
        MatrixType *Pn0         = NULL;
        MatrixType *Hmatrix_m1  = NULL;
        MatrixType *Hmatrix_0   = NULL;

    };
}


#endif

