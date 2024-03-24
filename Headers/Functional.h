/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
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


#ifndef RMG_Functional_H
#define RMG_Functional_H 1

#include <string>
#include "BaseGrid.h"
#include "Lattice.h"
#include "TradeImages.h"
#include "vdW.h"
#include "rmg_mangling.h"

#if USE_LIBXC
#include "xc.h"
#endif

class Functional {

private:
    // BaseGrid class
    BaseGrid *Grid;

    // TradeImages object to use
    TradeImages *T;

    // Lattice object
    Lattice *L;

    int pbasis;     // Grid points on this processing node

    // Total number of grid points
    int N;

    int fd_order;

    bool gammaflag;
    static bool dft_set;
    static bool exx_started;
    double gau_scrlen;


    void gradcorr(double *rho, double *rho_core, double &etxc, double &vtxc, double *v);
    void gradcorr_spin(double *rho_up, double *rho_down, double *rho_core, double &etxc, double &vtxc, double *v_up, double *v_down);

public:
    Functional (BaseGrid &G, 
                Lattice &L, 
                TradeImages &T, 
                bool gamma_flag);

    ~Functional(void);

    double *vxc2, *v2cud;

    void set_dft_from_name_rmg(char *newdft_name);
    static const std::string & get_dft_name_rmg(void);
    void set_dft_from_name_rmg(std::string newdft_name);
    bool dft_is_gradient_rmg(void);
    bool dft_is_meta_rmg(void);
    bool dft_is_hybrid_rmg(void);
    bool igcc_is_lyp_rmg(void);
    bool dft_is_nonlocc_rmg(void);
    bool dft_has_finite_size_correction_rmg(void);
    void v_xc(double *rho, double *rho_core, double &etxc, double &vtxc, double *v, int nspin);
    void nlc_rmg(double *rho, double *rho_core, double &etxc, double &vtxc, double *v);
    static void start_exx_rmg(void);
    static void stop_exx_rmg(void);
    static bool is_exx_active(void);
    double get_exx_fraction_rmg(void);
    void set_exx_fraction_rmg(double);
    static double get_gau_parameter_rmg(void);
    static void set_gau_parameter_rmg(double p);
    static double get_screening_parameter_rmg(void);
    static void set_screening_parameter_rmg(double p);

    static std::string saved_dft_name;
    void stress_vdW_DF (double *rho, double *rho_core, int nspin, double *sigma);

};

#endif
