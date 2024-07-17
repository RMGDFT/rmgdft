/*
 *
 * Copyright 2019 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/


#ifndef RMG_prototypes_rmg_H
#define RMG_prototypes_rmg_H 1

#include "Kpoint.h"
#include "PulayMixing.h"
#include "LaplacianCoeff.h"
#include "Prolong.h"
#include "transition.h"
#include "GridObject.h"


extern PulayMixing *Pulay_rho;
extern PulayMixing *Pulay_orbital;
extern BaseGrid *Rmg_G;
extern TradeImages *Rmg_T;
extern Lattice Rmg_L;
extern MpiQueue *Rmg_Q;


template <typename KpointType>
void MolecularDynamics (Kpoint<KpointType> **Kptr, double * vxc, double * vh, double * vnuc,
             double * rho, double * rho_oppo, double * rhoc, double * rhocore);

template <typename OrbitalType> void GetNewRhoOne(State<OrbitalType> *psi, Prolong *P, double *work, double scale);
template <typename OrbitalType> void GetNewRho(Kpoint<OrbitalType> **Kpts, double *rho);
template <typename OrbitalType> void GetNewRhoPre(Kpoint<OrbitalType> **Kpts, double *rho);
template <typename OrbitalType> void GetNewRhoPost(Kpoint<OrbitalType> **Kpts, double *rho);
template <typename OrbitalType> void GetAugRho(Kpoint<OrbitalType> **Kpts, double *rho);
#if HIP_ENABLED
template <typename OrbitalType> void GetNewRhoGpu(Kpoint<OrbitalType> **Kpts, double *rho);
template <typename OrbitalType> void GetNewRhoGpuOne(State<OrbitalType> *sp, Prolong *P, double scale);
#endif
template <typename OrbitalType> void Init (double * vh, double * rho, double * rho_oppo, double * rhocore, double * rhoc,
           double * vnuc, double * vxc, Kpoint<OrbitalType> **Kptr);
template <typename OrbitalType> void Reinit (double * vh, double * rho, double * rho_oppo, double * rhocore, double * rhoc,
           double * vnuc, double * vxc, Kpoint<OrbitalType> **Kptr);

template <typename OrbitalType> void Relax (
                   int steps,
                   spinobj<double> &vxc,
                   fgobj<double> &vh,
                   fgobj<double> &vnuc,
                   spinobj<double> &rho,
                   fgobj<double> &rhocore,
                   fgobj<double> &rhoc,
                   Kpoint<OrbitalType> **Kptr);

template <typename OrbitalType> void NEB_relax (int steps, double * vxc, double * vh, double * vnuc,
              double * rho, double * rho_oppo, double * rhocore, double * rhoc, Kpoint<OrbitalType> **Kptr);
template <typename OrbitalType> bool Quench (Kpoint<OrbitalType> **Kptr, bool compute_forces);
template <typename OrbitalType> bool Scf (spinobj<double> &vxc_in, fgobj<double> &vh_in,
                    double *vh_ext, int spin_flag, int boundaryflag, Kpoint<OrbitalType> **Kptr,
                    std::vector<double>& RMSdV);

template <typename KpointType> void AppNls(Kpoint<KpointType> *kpoint, KpointType *sintR,
            KpointType *psi, KpointType *nv, KpointType *ns, int first_state, int num_states);
template <typename KpointType> void AppS(Kpoint<KpointType> *kpoint, KpointType *sintR,
            KpointType *psi, KpointType *ns, int first_state, int num_states);
template <typename OrbitalType> double EnergyCorrection (Kpoint<OrbitalType> **Kptr,
          double *rho, double *new_rho, double *vh, double *vh_in);
template <typename OrbitalType> bool Nscf (double * vxc, double *vxc_correct, double * vh, double *vh_in, double *vh_ext,
          double * vnuc, double * rho, double * rho_oppo, double * rhocore, double * rhoc, int spin_flag,
          int boundaryflag, Kpoint<OrbitalType> **Kptr, std::vector<double>& RMSdV);
template <typename OrbitalType, typename CalcType> void MgEigState (Kpoint<OrbitalType> *kptr,
State<OrbitalType> * sp, double * vtot_psi, double *coarse_vtot, double *vxc_psi, OrbitalType *nv, OrbitalType *ns, int vcycle);

template <typename OrbitalType, typename CalcType>
void PotentialAcceleration(Kpoint<OrbitalType> *kptr, State<OrbitalType> *sp, double *vtot_psi, double *nvtot_psi, CalcType *tmp_psi_t, OrbitalType *saved_psi);
void PotentialAccelerationWait(int istate, int nstates, int skip);
void PotentialAccelerationReset(int skip);


template <typename KpointType>
void OutputDos (Kpoint<KpointType> **Kptr);
template <typename KpointType>
void OutputEigenvalues (Kpoint<KpointType> **Kptr, int ikbs, int iscf);
template <typename KpointType>
void ReadData (char *name, double * vh, double * rho, double * vxc, Kpoint<KpointType> ** Kptr);
template <typename KpointType>
void Read_nsocc (char *name, Kpoint<KpointType> * Kptr);
template <typename KpointType>
void Write_nsocc (char *name, Kpoint<KpointType> * Kptr);
template <typename KpointType>
void ExtrapolateOrbitals (char *name, Kpoint<KpointType> ** Kptr);
template <typename KpointType>
void GetOppositeEigvals (Kpoint<KpointType> **Kptr);
template <typename KpointType>
void GetOppositeOccupancies (Kpoint<KpointType> **Kptr);
template <typename StateType>
void LcaoGetAwave (StateType *psi, ION *iptr, int awave_idx, int l, int m, double coeff, double *kvec);
void LcaoGetRho (double * arho_f);
template <typename KpointType>
void GetTe (spinobj<double> &rho, fgobj<double> &rhocore, fgobj<double> &rhoc, fgobj<double> &vh, spinobj<double> &vxc, Kpoint<KpointType> ** Kptr , int ii_flag);
template <typename KpointType>
void WriteRestart (char *name, double * vh, double * rho, double * rho_oppo, double * vxc, Kpoint<KpointType> ** Kptr);
template <typename KpointType>
void WriteBGW_Wfng (int kpt, Kpoint<KpointType> * Kptr);
void WriteBGW_Rhog (double * rho, double * rho_oppo);
template <typename KpointType>
void WriteBGW_VxcEig (int kpt, double * vxc, Kpoint<KpointType> * Kptr);

template <typename KpointType>
void WriteData (int fhand, double * vh, double * rho, double * vxc, Kpoint<KpointType> ** Kptr);
template <typename KpointType>
void WriteSerialData (std::string& name, double * vh, double * rho, double * vxc, Kpoint<KpointType> ** Kptr);
template <typename KpointType>
void ReadSerialData (std::string& name, double * vh, double * rho, double * vxc, Kpoint<KpointType> ** Kptr);
template <typename KpointType>
double Fill (Kpoint<KpointType> **Kptr, double width, double nel, double mix, int num_st, int occ_flag, int mp_order);
template <typename KpointType>
double FillTetra(Kpoint<KpointType> **Kptr);

template <typename OrbitalType>
    void Force (double * rho, double * rho_oppo, double * rhoc, double * vh, double *vh_in,
        double * vxc, double *vxc_in, double * vnuc, Kpoint<OrbitalType> **Kptr);
template <typename OrbitalType>
    void Nlforce (double *, double *, Kpoint<OrbitalType> **Kptr, double *force);

template <typename OrbitalType> void nlforce_par_gamma (OrbitalType * par_gamma, int ion, int nh, double *force);
template <typename OrbitalType> void nlforce_par_omega (OrbitalType * par_omega, int ion, int nh, double *force);


template <typename OrbitalType> void PartialBetaxpsi (int ion, fftw_plan p2, double * newsintR_x, double * newsintR_y,
                       double * newsintR_z, double * newsintI_x, double * newsintI_y, double * newsintI_z,
                       ION * iptr,Kpoint<OrbitalType> **Kptr);

template <typename OrbitalType> void GetGamma (OrbitalType *gammaR, int ion, int nh , Kpoint<OrbitalType> **Kptr);


template <typename OrbitalType> void PartialGamma (int kpt,
                    int ion, OrbitalType * par_gammaR, OrbitalType * par_omegaR, int nion, int nh,
                    Kpoint<OrbitalType> **kptr, int state_start, int state_end,
                    OrbitalType *sint_derx, OrbitalType *sint_dery, OrbitalType *sint_derz);

template <typename OrbitalType> void AssignDerweight (Kpoint<OrbitalType> *kptr, SPECIES * sp, int ion, fftw_complex * beptr, OrbitalType
*rtptr);

template <typename KpointType>
void BandStructure(Kpoint<KpointType> ** Kptr, double *vh, double *vxc, double *vnuc, std::vector<bool> exclude_bands);

template <typename T> void nlforce_par_Q (double *veff, double *vxc, T *gamma_allions, double *forces, int num_owned_ions, int *owned_ion_list);

#endif
