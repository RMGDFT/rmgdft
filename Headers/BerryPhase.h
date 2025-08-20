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


#ifndef RMG_BerryPhase_H
#define RMG_BerryPhase_H 1


class BerryPhase
{
private:
    
public:
    int num_kort, num_kort_pe, num_kpp;
    int kort_start;
    int BerryPhase_dir;
    double efield_mag=0.0;
    double eps = 1.0e-10;
    double vel;
    std::complex<double> vel_C;
    double eai; //efield dot lattice vector in berryphase_dri
    double pol_elec, pol_ion, pol_tot;
    int nband_occ;
    int pbasis, pbasis_noncoll;
    size_t wfc_size;
    std::vector<double> kweight_string;
    std::complex<double> *mat, *psi_k0 = NULL;
    BerryPhase();

    ~BerryPhase(void);

    void init(Kpoint<std::complex<double>> **Kptr);;
    void init(Kpoint<double> **Kptr);;
    void Calc_Gnk(Kpoint<std::complex<double>> **Kptr);
    void Calc_Gnk(Kpoint<double> **Kptr);
    void Apply_BP_Hpsi(Kpoint<std::complex<double>> *kptr, int num_states, std::complex<double> *psi, std::complex<double> *h_psi);
    void Apply_BP_Hpsi(Kpoint<double> *kptr, int num_states, double *psi, double *h_psi);
    void CalcBP(Kpoint<std::complex<double>> **Kptr);
    void CalcBP(Kpoint<double> **Kptr);
    void psi_x_phase(std::complex<double> *psi_k0, double gr[3], int nband);
    void CalcBP_Skk1(Kpoint<std::complex<double>> **Kptr, int tddft_start_states, std::complex<double> *mat_glob, Scalapack &);
    void CalcBP_Skk1(Kpoint<double> **Kptr, int tddft_start_states, double *mat_glob, Scalapack &);
    void CalcBP_tddft(Kpoint<double> **Kptr, double &tot_bp_pol, double *, Scalapack &);
    void CalcBP_tddft(Kpoint<std::complex<double>> **Kptr, double &tot_bp_pol, std::complex<double> *mat_glob, Scalapack &);
    void tddft_Xml(Kpoint<double> **Kptr, int tddft_start_states, double *mat_glob, Scalapack &);
    void tddft_Xml(Kpoint<std::complex<double>> **Kptr, int tddft_start_states, std::complex<double> *mat_glob, Scalapack &);

};

#endif
