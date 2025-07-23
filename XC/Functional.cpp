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

#include <math.h>
#include <float.h>
#include <complex>
#include <iostream>
#include <fstream>
#include <fcntl.h>
#include <sys/stat.h>
#include <algorithm>

#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "Functional.h"
#include "FiniteDiff.h"
#include "RmgException.h"
#include "RmgTimer.h"
#include "FiniteDiff.h"
#include "RmgSumAll.h"
#include "RmgParallelFft.h"
#include "transition.h"
#if USE_LIBXC
#include "xc.h"
#endif
#include "rmg_mangling.h"

#define set_dft_from_name       RMG_FC_MODULE(funct,set_dft_from_name,mod_FUNCT,SET_DFT_FROM_NAME)
#define get_dft_name            RMG_FC_MODULE(funct,get_dft_name,mod_FUNCT,GET_DFT_NAME)
#define start_exx               RMG_FC_MODULE(dft_setting_routines,start_exx,mod_FUNCT,START_EXX)
#define stop_exx                RMG_FC_MODULE(dft_setting_routines,stop_exx,mod_FUNCT,STOP_EXX)
#define dft_is_gradient         RMG_FC_MODULE(dft_setting_routines,dft_is_gradient,mod_FUNCT,DFT_IS_GRADIENT)
#define dft_is_meta             RMG_FC_MODULE(dft_setting_routines,dft_is_meta,mod_FUNCT,DFT_IS_META)
#define dft_is_hybrid           RMG_FC_MODULE(dft_setting_routines,dft_is_hybrid,mod_FUNCT,DFT_IS_HYBRID)
#define get_exx_fraction        RMG_FC_MODULE(dft_setting_routines,xclib_get_exx_fraction,mod_FUNCT,GET_EXX_FRACTION)
#define set_exx_fraction        RMG_FC_MODULE(dft_setting_routines,xclib_set_exx_fraction,mod_FUNCT,SET_EXX_FRACTION)
#define set_rmg_epsg_guard      RMG_FC_MODULE(dft_setting_routines,set_rmg_epsg_guard,mod_FUNCT,SET_RMG_EPSG_GUARD)
#define igcc_is_lyp             RMG_FC_MODULE(dft_setting_routines,igcc_is_lyp,mod_FUNCT,IGCC_IS_LYP)
#define dft_has_finite_size_correction RMG_FC_MODULE(dft_setting_routines,dft_has_finite_size_correction,mod_FUNCT,DFT_HAS_FINITE_SIZE_CORRECTION)
#define dft_is_nonlocc          RMG_FC_MODULE(funct,dft_is_nonlocc,mod_FUNCT,DFT_IS_NONLOCC)
#define xc_spin                 RMG_FC_MODULE(funct,xc_spin,mod_FUNCT,XC)
#define nlc                     RMG_FC_MODULE(funct,nlc,mod_FUNCT,NLC)
#define xc_gcx                  RMG_FC_GLOBAL(xc_gcx, XC_GCX)
#define xc                      RMG_FC_GLOBAL(xc, XC)
#define xc_metagcx              RMG_FC_GLOBAL(xc_metagcx, XC_METAGCX)
#define xc_lsda                 RMG_FC_MODULE(qe_drivers_lda_lsda,xc_lsda,mod_FUNCT,XC_LSDA)
#define gcx_spin                RMG_FC_MODULE(qe_drivers_gga,gcx_spin,mod_FUNCT,GCX_SPIN)
#define gcc_spin_more           RMG_FC_MODULE(qe_drivers_gga,gcc_spin_more,mod_FUNCT,GCC_SPIN_MORE)
#define gcc_spin                RMG_FC_MODULE(qe_drivers_gga,gcc_spin,mod_FUNCT,GCC_SPIN)
#define tau_xc                  RMG_FC_MODULE(qe_drivers_mgga,tau_xc,mod_FUNCT,TAU_XC)
#define tau_xc_spin             RMG_FC_MODULE(qe_drivers_mgga,tau_xc_spin,mod_FUNCT,TAU_XC_SPIN)
#define get_inlc                RMG_FC_MODULE(funct,get_inlc,mod_FUNCT,GET_INLC)
#define get_gau_parameter       RMG_FC_MODULE(dft_setting_routines,get_gau_parameter,mod_FUNCT,GET_GAU_PARAMETER)
#define set_gau_parameter       RMG_FC_MODULE(dft_setting_routines,set_gau_parameter,mod_FUNCT,SET_GAU_PARAMETER)
#define get_screening_parameter       RMG_FC_MODULE(dft_setting_routines,get_screening_parameter,mod_FUNCT,GET_SCREENING_PARAMETER)
#define set_screening_parameter       RMG_FC_MODULE(dft_setting_routines,set_screening_parameter,mod_FUNCT,SET_SCREENING_PARAMETER)

double *Functional::ke_density;
double *Functional::ke_taur;
double *Functional::ke_taur_wf;

extern "C" void set_dft_from_name( const char *name, std::size_t len );
extern "C" char *get_dft_name(void);
extern "C" void start_exx(void);
extern "C" void stop_exx(void);
extern "C" bool dft_is_gradient(void);
extern "C" bool dft_is_meta(void);
extern "C" bool dft_is_hybrid(void);
extern "C" double get_exx_fraction(void);
extern "C" void set_exx_fraction(double *frac);
extern "C" void set_rmg_epsg_guard(double *epsg_guard);
extern "C" bool igcc_is_lyp(void);
extern "C" bool dft_has_finite_size_correction(void);
extern "C" bool dft_is_nonlocc(void);
extern "C" void xc (int *length, int *i1, int *i2, double *rho, double *ex, double *ec, double *vx, double *vc, bool *gargs);
extern "C" void xc_lsda (int *length, double *rho, double *zeta, double *ex, double *ec, double *vx, double *vc);
extern "C" void tau_xc (int *length, double *arho, double *grho2, double *atau, double *ex,
                        double *ec, double *v1x, double *v2x, double *v3x,
                        double *v1c, double *v2c, double *v3c);
extern "C" void xc_metagcx( int *length, int *ione, int *np, double *rho, double *grhof, 
                   double *ked, double *ex, double *ec, double *v1x, double *v2x, double *v3x,
                   double *v1c, double *v2c, double *v3c, bool *gargs );

extern "C" void xc_spin (double *rho, double *zeta, double *ex, double *ec, double *vxup, double *vxdw, double *vcup, double *vcdw);

extern "C" void nlc (double *rho_valence, double *rho_core, int *nspin, double *ec, double *vx, double *vc);
extern "C" void xc_gcx (int *length, int*nspin, double *rho, double *grho, double *sx, double *sc, 
                                  double *v1x, double *v2x, double *v1c, double *v2c, double *v2cud, bool *gargs);
extern "C" void gcx_spin(double *rhoup, double *rhodown, double *grhoup, double *grhodown, double *sx,
                                  double *v1xup, double *v1xdw, double *v2xup, double *v2xdw);
extern "C" void gcc_spin_more( double *arho_up, double *arho_down,  double *grhoup, double *grhodw, double *grhoud,
                                  double *sc, double *v1cup, double *v1cdw, double *v2cup, double *v2cdw, double *v2cud );
extern "C" void gcc_spin( double *arho, double *zeta, double *grh2, double *sc, double *v1cup, double *v1cdw, double *v2c );
extern "C" void tb09cxc(double *rho, double *grho, double *tau, double *ex, double *ec, 
            double *v1x, double *v2x,double *v3x, double *v1c, double *v2c, double *v3c);
extern "C" int get_inlc(void);
extern "C" double get_gau_parameter(void);
extern "C" void set_gau_parameter(double *);
extern "C" double get_screening_parameter(void);
extern "C" void set_screening_parameter(double *);
#if __LIBXC
#define xclib_init_libxc        RMG_FC_MODULE(dft_setting_routines,xclib_init_libxc,mod_FUNCT,XCLIB_INIT_LIBXC)
extern "C" void xclib_init_libxc(int *nspin, bool *domag);
#endif

bool Functional::dft_set=false;
bool Functional::exx_started=false;
std::string Functional::saved_dft_name;

double SMALL_CHARGE=1.0e-10;
double SMALL_MAG=1.0e-20;

void CToF_2d(int n, double *a, double *b)
{
    int kk=0;
    for(int i=0;i<3;i++)
    {
      for(int j=0;j < n;j++)b[i + j*3] = a[kk++];
    }
}

// Utility function to convert a c-string to a fixed length fortran string
void CstrToFortran(char* fstring, std::size_t flen, const char* cstring)
{
    std::size_t clen = std::strlen(cstring);
    std::size_t copylen = std::min(clen, flen);

    if (clen > flen)
    {
        // truncation error
    }

    std::copy(cstring, cstring + copylen, fstring);
    std::fill(fstring + copylen, fstring + flen, ' ');
}



Functional::Functional (
            BaseGrid &G,
            Lattice &L,
            TradeImages &T,
            bool gamma_flag)
{
    RmgTimer RT0("5-Functional");
#if __LIBXC
    static bool initialized;
    bool domag = false;
    if(!initialized)
    {
        xclib_init_libxc(&ct.nspin, &domag);
        initialized = false;
    }
#endif
    this->Grid = &G;
    this->T = &T;
    this->L = &L;
    this->gammaflag = gamma_flag;

    this->fd_order = ct.kohn_sham_fd_order;
    this->pbasis = G.get_P0_BASIS(G.default_FG_RATIO);
    this->N = G.get_GLOBAL_BASIS(G.default_FG_RATIO);

    if(dft_set && dft_is_hybrid() && exx_started) start_exx();
    if(dft_set && dft_is_hybrid() && !exx_started) stop_exx();

    if(this->dft_is_gradient_rmg()) 
    {
        vxc2 = new double[2 * this->pbasis]();
        v2cud = new double[this->pbasis]();
    }
    else
    {
        vxc2 = NULL;
        v2cud = NULL;
    }
    
}




Functional::~Functional(void)
{

    if(this->dft_is_gradient_rmg()) 
    {
        if(vxc2) delete [] vxc2;
        if(v2cud) delete [] v2cud;
        vxc2 = NULL;
        v2cud = NULL;
    }

}


void Functional::set_dft_from_name_rmg(char *newdft_name) 
{
    if(!this->dft_set) {
        set_dft_from_name(newdft_name, std::strlen(newdft_name) );
        saved_dft_name = newdft_name;
    }
    else {
        std::cout << "Warning! You have attempted to reset the dft type which is a programming error. Ignoring." << std::endl;
    }
    this->dft_set = true;
}

const std::string &Functional::get_dft_name_rmg(void)
{
    if(!dft_set) {
        throw RmgFatalException() << "Error! get_dft_name called before dft type was set." << " in " << __FILE__ << " at line " << __LINE__ << "\n";
    }
    return saved_dft_name;
}

void Functional::start_exx_rmg(void)
{
    exx_started = true;
    start_exx();
}

void Functional::stop_exx_rmg(void)
{
    exx_started = false;
    stop_exx();
}

bool Functional::is_exx_active(void)
{
    return exx_started;
}

void Functional::set_dft_from_name_rmg(std::string newdft_name)
{
    if(!this->dft_set) {
        set_dft_from_name(newdft_name.c_str(), std::strlen(newdft_name.c_str()) );
        saved_dft_name = newdft_name;
    }
    else {
        std::cout << "Warning! You have attempted to reset the dft type which is a programming error. Ignoring." << std::endl;
    }
    this->dft_set = true;
}

bool Functional::dft_is_gradient_rmg(void)
{
    return dft_is_gradient();
}

bool Functional::dft_is_meta_rmg(void)
{
    return dft_is_meta();
}

bool Functional::dft_is_hybrid_rmg(void)
{
    return dft_is_hybrid();
}

double Functional::get_exx_fraction_rmg(void)
{
    return get_exx_fraction();
}

void Functional::set_exx_fraction_rmg(double frac)
{
    set_exx_fraction(&frac);
}

void Functional::set_epsg_guard(double epsg_guard)
{
    set_rmg_epsg_guard(&epsg_guard);
}
double Functional::get_gau_parameter_rmg(void)
{
    return get_gau_parameter();
}

void Functional::set_gau_parameter_rmg(double p)
{
    set_gau_parameter(&p);
}

double Functional::get_screening_parameter_rmg(void)
{
    return get_screening_parameter();
}

void Functional::set_screening_parameter_rmg(double p)
{
    set_screening_parameter(&p);
}

bool Functional::igcc_is_lyp_rmg(void)
{
    return igcc_is_lyp();
}

bool Functional::dft_has_finite_size_correction_rmg(void)
{
    return dft_has_finite_size_correction();
}

bool Functional::dft_is_nonlocc_rmg(void)
{
    return dft_is_nonlocc();
}

extern Kpoint<double> **Kptr_g;
void Functional::v_xc(double *rho_in, double *rho_core, double &etxc, double &vtxc, double *v_out, int nspin)
{
    RmgTimer RT0("5-Functional");
    RmgTimer RT1("5-Functional: vxc");
    int ione = 1;
    int itwo = 2;
    int ifour = 4;
    const double epsr=1.0e-6;
    bool gargs=false;

    etxc = 0.0;
    vtxc = 0.0;
    for(int ix = 0;ix < nspin*this->pbasis;ix++) v_out[ix] = 0.0;

    if(dft_is_meta())
    {
        int wf_pbasis = Rmg_G->get_P0_BASIS(1);
        wfobj<double> kdetau_c;
        fgobj<double> kdetau_f;
        kdetau_c.set(0.0);
        if(ct.scf_steps >= 0)
        {
            for(int ik = 0; ik < ct.num_kpts_pe; ik++) Kptr_g[ik]->KineticEnergyDensity(kdetau_c.data());
//            FftInterpolation(*Rmg_G, kdetau_c.data(), kdetau_f.data(), 2, false);
            int ratio = Rmg_G->default_FG_RATIO;
            Prolong P(2, ct.prolong_order, 0.0, *Rmg_T,  Rmg_L, *Rmg_G);
            int dimx = Rmg_G->get_PX0_GRID(ratio);
            int dimy = Rmg_G->get_PY0_GRID(ratio);
            int dimz = Rmg_G->get_PZ0_GRID(ratio);
            int half_dimx = Rmg_G->get_PX0_GRID(1);
            int half_dimy = Rmg_G->get_PY0_GRID(1);
            int half_dimz = Rmg_G->get_PZ0_GRID(1);
            P.prolong(kdetau_f.data(), kdetau_c.data(), dimx, dimy, dimz, half_dimx, half_dimy, half_dimz);

        }
//  Need to update GetNewRho if we want to compute this on fine grid but not clear if it's necessary
        else
        {
            for(int ix=0;ix < this->pbasis;ix++) kdetau_f[ix] = this->ke_density[ix];
        }
        v_xc_meta(rho_in, rho_core, etxc, vtxc, v_out, kdetau_f.data(), nspin);
        return;
    }

    double rhoneg[2]{0.0,0.0};
    double *rho_up=NULL, *rho_down=NULL;
    double *v_up=NULL, *v_down=NULL;



    // First get the local exchange and correlation
    RmgTimer *RT2 = new RmgTimer("5-Functional: vxc local");
    if(nspin==1) {

        double etxcl=0.0, vtxcl=0.0, rhonegl=0.0;
        // spin unpolarized  
        fgobj<double> ex, ec, vx, vc;
        int length = ex.pbasis;
        fgobj<double> trho;
        for(int ix=0;ix < this->pbasis;ix++) trho[ix] = rho_in[ix] + rho_core[ix];
        xc ( &length, &ione, &ione, trho.data(), ex.data(), ec.data(), vx.data(), vc.data(), &gargs);
        for(int ir=0;ir < ex.pbasis;ir++)
        {
            v_out[ir] = vx[ir] + vc[ir];
            etxc = etxc + ( ex[ir] + ec[ir])*trho[ir];
            vtxc = vtxc + v_out[ir]*rho_in[ir];
            if(rho_in[ir] < 0.0) rhonegl = rhonegl-rho_in[ir];
        }
    } 
    else if(nspin == 2)
    {
        double etxcl=0.0, vtxcl=0.0;
        // spin polarized
        fgobj<double> ex, ec;
        spinobj<double> vx, vc, trho;

        vx.set(0.0);
        vc.set(0.0);
        ex.set(0.0);
        ec.set(0.0);

        // QE routines expect total charge + magnetization stored sequentially.
        for(int ix=0;ix < this->pbasis;ix++)
            trho.up[ix] = rho_in[ix] + rho_in[ix+this->pbasis] + rho_core[ix];

        // for collinear case, spin up and down are in different processor groups.
        if(pct.spinpe == 0) {
            for(int ix=0;ix < this->pbasis;ix++) trho.dw[ix] = rho_in[ix] - rho_in[ix+this->pbasis];
            rho_up = rho_in;
            rho_down = &rho_in[this->pbasis];
            v_up = v_out;
            v_down = &v_out[this->pbasis];
        }
        else {
            for(int ix=0;ix < this->pbasis;ix++) trho.dw[ix] = rho_in[ix+this->pbasis] - rho_in[ix];
            rho_down = rho_in;
            rho_up = &rho_in[this->pbasis];
            v_down = v_out;
            v_up = &v_out[this->pbasis];
        }

        int length = this->pbasis;
        xc ( &length, &itwo, &itwo, trho.data(), ex.data(), ec.data(), vx.data(), vc.data(), &gargs);
        for(int ix=0;ix < this->pbasis;ix++)
        {
            double atrho = fabs(trho[ix]);
            if(atrho > SMALL_CHARGE)
            {
                v_up[ix] = vx.up[ix] + vc.up[ix];
                v_down[ix] = vx.dw[ix] + vc.dw[ix];
                etxcl += (ex[ix] + ec[ix])*trho.up[ix];
                trho.up[ix] -= rho_core[ix];
                vtxcl += ((v_up[ix] + v_down[ix])*trho.up[ix] +
                          (v_up[ix] - v_down[ix])*trho.dw[ix] )*0.5;
            }
        }
        etxc += etxcl;
        vtxc += vtxcl;

    }
    else if(nspin == 4)
    {
        double etxcl=0.0, vtxcl=0.0;
        // spin-orbit
        fgobj<double> ex, ec;
        // wrap rho_in and v_out storage with spinobj for convenience
        spinobj<double> vx, vc, rho(rho_in), v(v_out);

        for(int ix=0;ix < this->pbasis;ix++) rho_in[ix] += rho_core[ix];
        xc(&this->pbasis, &ifour, &itwo, rho_in, ex.data(), ec.data(), vx.data(), vc.data(), &gargs);

        for(int ix=0;ix < this->pbasis;ix++)
        {
            double arho = std::abs(rho_in[ix]);
            if(arho < SMALL_CHARGE)
            {
                v.c0[ix] = 0.0;
                v.cx[ix] = 0.0;
                v.cy[ix] = 0.0;
                v.cz[ix] = 0.0;
                continue;
            }

            double vs = 0.5*( vx.up[ix] + vc.up[ix] - vx.dw[ix] - vc.dw[ix]);
            v.c0[ix] = (0.5*( vx.up[ix] + vc.up[ix] + vx.dw[ix] + vc.dw[ix]));
            double vtxc24;
            double amag = sqrt(rho.cx[ix]*rho.cx[ix] + rho.cy[ix]*rho.cy[ix] + rho.cz[ix]*rho.cz[ix]);
            if(amag > SMALL_MAG)
            {
                v.cx[ix] = vs * rho.cx[ix] / amag;
                v.cy[ix] = vs * rho.cy[ix] / amag;
                v.cz[ix] = vs * rho.cz[ix] / amag;
                vtxc24 = v[ix +   this->pbasis] * rho.cx[ix] +
                         v[ix + 2*this->pbasis] * rho.cy[ix] +
                         v[ix + 3*this->pbasis] * rho.cz[ix];
            }
            else
            {
                v.cx[ix] = 0.0;
                v.cy[ix] = 0.0;
                v.cz[ix] = 0.0;
                vtxc24 = 0.0;
            }
            etxcl += arho * (ex[ix] + ec[ix] );
            rho_in[ix]  = rho_in[ix] - rho_core[ix];
            //IF (rho%of_r(ir,1) < 0.D0 )  rhoneg1 = rhoneg1 - rho%of_r(ir,1)
            //IF ( amag / arho  > 1.D0 )  rhoneg2 = rhoneg2 + 1.D0/omega
            vtxcl += vtxc24 + v.c0[ix] * rho_in[ix];
        }
        etxc += etxcl;
        vtxc += vtxcl;
    }
    delete RT2;


    // Next add in any gradient corrections
    RmgTimer *RT3 = new RmgTimer("5-Functional: vxc grad");
    if(nspin == 1) {
        this->gradcorr(rho_in, rho_core, etxc, vtxc, v_out);
    }
    else if(nspin == 2) {
        this->gradcorr_spin(rho_up, rho_down, rho_core, etxc, vtxc, v_up, v_down);
    }
    delete RT3;

    // And finally any non-local corrections
    RmgTimer *RT4 = new RmgTimer("5-Functional: vxc nonlocal");
    if(this->dft_is_nonlocc_rmg()) {
        if(nspin == 4) 
        {
            throw RmgFatalException() << "vdw with noncollinear not programed. " << " in " << __FILE__ << " at line " << __LINE__ << "\n";

        }
        double netxc=0.0, nvtxc=0.0;
        this->nlc_rmg(rho_in, rho_core, netxc, nvtxc, v_out);
        vtxc += nvtxc;
        etxc += netxc;
    }
    delete RT4;

    vtxc = vtxc * L->omega / (double)this->N;
    etxc = etxc * L->omega / (double)this->N;

    vtxc = RmgSumAll(vtxc, this->T->get_MPI_comm());
    etxc = RmgSumAll(etxc, this->T->get_MPI_comm());

    if(Rmg_G->default_FG_RATIO > 1)
    {
        for(int is = 0; is < nspin; is++)
            FftFilter(&v_out[is*pbasis], *fine_pwaves, *coarse_pwaves, LOW_PASS);
    }

}

void Functional::v_xc_meta(double *rho_in, double *rho_core, double &etxc, double &vtxc, double *v, double *ked, int nspin)
{
    double eps8 = 1.0e-8, eps12 = 1.0e-12;
    double tpiba = 2.0 * PI / Rmg_L.celldm[0];
    etxc = 0.0;
    vtxc = 0.0;
    double rhoneg[2];
    int np = 1;
    int ione = 1;
    int itwo = 2;
    bool gargs = false;
    if(nspin == 2) np=3;

    if(nspin == 1)
    {
        fgobj<double> rho, grho2, lrho, d2rho;
        fgobj<double> hx, hy, hz, dh1, dh2, dh3;
        fgobj<double> ex, ec, v1x, v2x, v3x, v1c, v2c, v3c;
        double *gx = new double[3*this->pbasis];
        double *gy = gx + this->pbasis;
        double *gz = gy + this->pbasis;
        double *grhof = new double[3*this->pbasis];
        for(int ix=0;ix < rho.pbasis;ix++)rho[ix] = rho_in[ix] + rho_core[ix];
        ApplyGradient (rho.data(), gx, gy, gz, fd_order, "Fine");
        // Have to convert 2D array to Fortran order for QE routine.
        CToF_2d(this->pbasis, gx, grhof);

        for(int ix=0;ix < this->pbasis;ix++) grho2[ix] = gx[ix]*gx[ix] + gy[ix]*gy[ix] + gz[ix]*gz[ix];
        //ApplyLaplacian (rho.data(), lrho.data(), fd_order, "Fine");

        int length = this->pbasis;
        xc_metagcx( &length, &ione, &np, rho.data(), grhof, ked, 
                   ex.data(), ec.data(), v1x.data(), v2x.data(), v3x.data(),
                   v1c.data(), v2c.data(), v3c.data(), &gargs );

        for(int ix=0;ix < this->pbasis;ix++)
        {
            double arho = std::abs(rho[ix]);
            double atau = ked[ix];
            if ((arho > eps8) && (grho2[ix] > eps12) && (std::abs(atau) > eps8))
//if(1)
            {
                v[ix] +=  (v1x[ix] + v1c[ix]);
                // h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
                hx[ix] =  (v2c[ix] + v2x[ix])*gx[ix];
                hy[ix] =  (v2c[ix] + v2x[ix])*gy[ix];
                hz[ix] =  (v2c[ix] + v2x[ix])*gz[ix];
                ke_taur[ix] =  (v3x[ix] + v3c[ix]) * 0.5;
                etxc = etxc +  (ex[ix] + ec[ix]); // * segno
                vtxc = vtxc +  0.5*(v1x[ix]+v1c[ix])*rho_in[ix];
            }
            else
            {
                v[ix] = 0.0;
                hx[ix] = 0.0;
                hy[ix] = 0.0;
                hz[ix] = 0.0;
                ke_taur[ix] = 0.0;
            }
        }

#if 1
        ApplyGradient (hx.data(), dh1.data(), dh2.data(), dh3.data(), fd_order, "Fine");
        for(int ix=0;ix < rho.pbasis;ix++)
        {
             double dh = hx[ix]*dh1[ix];
dh = dh1[ix];
             v[ix] -= dh;
             vtxc -= dh * rho_in[ix];
        }
        ApplyGradient (hy.data(), dh1.data(), dh2.data(), dh3.data(), fd_order, "Fine");
        for(int ix=0;ix < rho.pbasis;ix++)
        {
             double dh = hy[ix]*dh2[ix];
dh = dh2[ix];
             v[ix] -= dh;
             vtxc -= dh * rho_in[ix];
        }
        ApplyGradient (hz.data(), dh1.data(), dh2.data(), dh3.data(), fd_order, "Fine");
        for(int ix=0;ix < rho.pbasis;ix++)
        {
             double dh = hz[ix]*dh3[ix];
dh = dh3[ix];
             v[ix] -= dh;
             vtxc -= dh * rho_in[ix];
        }
#endif
        delete [] grhof;
        delete [] gx;
    }
    else if(nspin == 2)
    {
#if 0
        fgobj<double> ex, ec;
        spinobj<double> v1x, v2x, v3x, v1c, v2c, v3c;
        spinobj<double> rho, grho2;
        double *rho_up, *rho_dw, *v_up, *v_dw;

        double *gx_up = new double[6*this->pbasis];
        double *gy_up = gx_up + this->pbasis;
        double *gz_up = gy_up + this->pbasis;
        double *gx_dw = gz_up + this->pbasis;
        double *gy_dw = gx_dw + this->pbasis;
        double *gz_dw = gy_dw + this->pbasis;
        double *grhof = new double[6*this->pbasis];
        if(pct.spinpe == 0) {
            rho_up = rho_in;
            rho_dw = &rho_in[this->pbasis];
            v_up = v;
            v_dw = &v[this->pbasis];
        }
        else {
            rho_dw = rho_in;
            rho_up = &rho_in[this->pbasis];
            v_dw = v;
            v_up = &v[this->pbasis];
        }
        for(int ix=0;ix < this->pbasis;ix++)
        {
            rho.up[ix] = rho_up[ix] + 0.5*rho_core[ix];
            rho.dw[ix] = rho_dw[ix] + 0.5*rho_core[ix];
        }
        ApplyGradient (rho.up.data(), gx_up, gy_up, gz_up, fd_order, "Fine");
        ApplyGradient (rho.dw.data(), gx_dw, gy_dw, gz_dw, fd_order, "Fine");
        CToF_2d(this->pbasis, gx_up, grhof);
        CToF_2d(this->pbasis, gx_dw, grhof+3*this->pbasis);
        int length = this->pbasis;
        xc_metagcx( &length, &itwo, &np, rho.data(), grhof, ked, 
                   ex.data(), ec.data(), v1x.data(), v2x.data(), v3x.data(),
                   v1c.data(), v2c.data(), v3c.data(), &gargs );

        for(int ix=0;ix < this->pbasis;ix++)
        {
            v_up[ix] = v1x.up[ix] + v1c.up[ix];
            v_dw[ix] = v1x.dw[ix] + v1c.dw[ix];

                // h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
                hx_up[ix] =  (v2x.up[ix] + v2c.up[ix])*gx_up[ix];
                hy_up[ix] =  (v2x.up[ix] + v2c.up[ix])*gy_up[ix];
                hz_up[ix] =  (v2x.up[ix] + v2c.up[ix])*gz_up[ix];
                hx_dw[ix] =  (v2x.dw[ix] + v2c.dw[ix])*gx_dw[ix];
                hy_dw[ix] =  (v2x.dw[ix] + v2c.dw[ix])*gy_dw[ix];
                hz_dw[ix] =  (v2x.dw[ix] + v2c.dw[ix])*gz_dw[ix];
                ke_taur[ix] =  (v3x.up[ix] + v3c.up[ix]) * 0.5;
                ke_taur[ix+this->pbasis] = (v3x.dw[ix] + v3c.dw[ix]) * 0.5;
                etxc = etxc +  (ex[ix] + ec[ix]);
                vtxc = vtxc +  (v1x.up[ix]+v1c.up[ix])*rho_in[ix] +
                               (v1x.dw[ix]+v1c.dw[ix])*rho_in[ix];

        }
#endif
    } 

    vtxc = vtxc * L->omega / (double)this->N;
    etxc = etxc * L->omega / (double)this->N;

    vtxc = RmgSumAll(vtxc, this->T->get_MPI_comm());
    etxc = RmgSumAll(etxc, this->T->get_MPI_comm());

    if(Rmg_G->default_FG_RATIO > 1)
    {
        int wf_pbasis = Rmg_G->get_P0_BASIS(1);
        for(int is = 0; is < nspin; is++)
        {
            FftFilter(&v[is*pbasis], *fine_pwaves, *coarse_pwaves, LOW_PASS);
            GetVtotPsi (&ke_taur_wf[is*wf_pbasis], 
                        &ke_taur[is*pbasis], Rmg_G->default_FG_RATIO);
        }
    }
}

// Applies non-local corrections for the correlation
void Functional::nlc_rmg(double *rho, double *rho_core, double &etxc, double &vtxc, double *v)
{
    int inlc = get_inlc();

    // No non-local correction just return
    if(inlc == 0) return;

    // inlc == 1 corresponds to vdW-DF1 and is the only one programmed currently
    if(inlc == 1) {
        Vdw *vdw = new Vdw (*this->Grid, *this->L, *this->T, 1, rho, rho_core, etxc, vtxc, v, this->gammaflag);
        delete vdw;
    }
    else {
        throw RmgFatalException() << "Non-local correlation correction type not programmed" << " in " << __FILE__ << " at line " << __LINE__ << "\n";

    }

}

// Applies any gradient corrections
void Functional::gradcorr(double *rho, double *rho_core, double &etxc, double &vtxc, double *v)
{
    if(!this->dft_is_gradient_rmg()) return;

    const double epsr=1.0e-6;
    const double epsg = 1.0e-10;
    double epsg_guard = ct.epsg_guard;
    int ione = 1;
    double etxcgc = 0.0;
    double vtxcgc = 0.0;
    bool gargs = false;

    fgobj<double> rhoout, sx, sc, v1x, v2x, v1c, v2c, d2rho;
    double *grho = new double[3*this->pbasis]();
    double *gx = grho;
    double *gy = gx + this->pbasis;
    double *gz = gy + this->pbasis;
    double *grhof = new double[3*this->pbasis]();
    double *vxc2 = this->vxc2;

    // Get rho plus rhocore
    for(int ix=0;ix < this->pbasis;ix++) rhoout[ix] = rho[ix] + rho_core[ix];

    // calculate the gradient of rho + rho_core
    RmgTimer *RT2 = new RmgTimer("5-Functional: apply gradient");
    ApplyGradient (rhoout.data(), gx, gy, gz, fd_order, "Fine");
    //FftGradientFine(rhoout.data(), gx, gy, gz);

    // and the Laplacian
    RmgTimer *RT3 = new RmgTimer("5-Functional: apply laplacian");
    //FftLaplacianFine(rhoout.data(), d2rho.data());
    ApplyLaplacian (rhoout.data(), d2rho.data(), fd_order, "Fine");
    delete RT3;

    // Have to convert 2D array to Fortran order for QE routine.
    CToF_2d(this->pbasis, grho, grhof);
    double *v2dummy=NULL;
    xc_gcx( &this->pbasis, &ione, rhoout.data(), grhof, sx.data(), sc.data(), v1x.data(),
            v2x.data(), v1c.data(), v2c.data(), v2dummy, &gargs);

    for(int k=0;k < this->pbasis;k++)
    {
        // ... first term of the gradient correction : D(rho*Exc)/D(rho)
        v[k] += v1x[k] + v1c[k];
        //  used later for second term of the gradient correction
        vxc2[k] = ( v2x[k] + v2c[k] );
        vtxcgc += (v1x[k] + v1c[k]) * (rhoout[k] - rho_core[k]);
        etxcgc += sx[k] + sc[k];
    } 

    // 
    // ... second term of the gradient correction :
    // ... \sum_alpha (D / D r_alpha) ( D(rho*Exc)/D(grad_alpha rho) )
    // 
    RmgTimer *RT5 = new RmgTimer("5-Functional: apply gradient");
    double *h = grhof;
    ApplyGradient (vxc2, h, &h[this->pbasis], &h[2*this->pbasis], fd_order, "Fine");
    delete RT5;

    double vtxcgc_1 = 0.0;
//#pragma omp parallel for reduction(+:vtxcgc_1)
    for(int ix=0;ix < this->pbasis;ix++) {
        double arho = fabs(rhoout[ix]);
        double grho2 = gx[ix]*gx[ix] + gy[ix]*gy[ix] + gz[ix]*gz[ix];
        if(arho > epsr && grho2 > epsg)
        {
            double gdot =  ( h[ix] * gx[ix] +
                    h[ix+this->pbasis] * gy[ix] + 
                    h[ix+2*this->pbasis] * gz[ix] ) ;
            v[ix] -= gdot;
            v[ix] -= vxc2[ix] * d2rho[ix];
            vtxcgc_1 -= rho[ix]*(gdot + vxc2[ix] * d2rho[ix]);
        }
    }
    vtxc = vtxc + vtxcgc + vtxcgc_1;
    etxc = etxc + etxcgc;

    delete [] grhof;
    delete [] grho;

}

// Applies gradient corrections for spin case
void Functional::gradcorr_spin(double *rho_up, double *rho_down, double *rho_core, double &etxc, double &vtxc, double *v_up, double *v_down)
{
    if(!this->dft_is_gradient_rmg()) return;

    int itwo = 2;
    bool gargs = false;
    double etxcgc = 0.0;
    double vtxcgc = 0.0;

    const double epsr=1.0e-6;
    const double epsg = 1.0e-10;
    double epsg_guard = ct.epsg_guard;

    double *grho_up = new double[6*this->pbasis];
    double *grho_down = grho_up + 3*this->pbasis;
    double *vxc2_up = this->vxc2;
    double *vxc2_down = vxc2_up + this->pbasis;
    double *v2cud = this->v2cud;
    double *rhoout_up = new double[this->pbasis];
    double *rhoout_down = new double[this->pbasis];

    double *gx_up = grho_up;
    double *gy_up = gx_up + this->pbasis;
    double *gz_up = gy_up + this->pbasis;
    double *gx_down = grho_down;
    double *gy_down = gx_down + this->pbasis;
    double *gz_down = gy_down + this->pbasis;

    // Get rho plus rhocore
    for(int ix=0;ix < this->pbasis;ix++) rhoout_up[ix] = rho_up[ix] + 0.5*rho_core[ix];
    for(int ix=0;ix < this->pbasis;ix++) rhoout_down[ix] = rho_down[ix] + 0.5*rho_core[ix];

    // calculate the gradient of rho + rho_core up
    RmgTimer *RT2 = new RmgTimer("5-Functional: apply gradient");
    ApplyGradient (rhoout_up, gx_up, gy_up, gz_up, fd_order, "Fine");
    ApplyGradient (rhoout_down, gx_down, gy_down, gz_down, fd_order, "Fine");

    delete RT2;

    fgobj<double> ex, ec;
    spinobj<double> v1x, v2x, v1c, v2c, trho;

    // QE routines expect total charge + magnetization stored sequentially.
    for(int ix=0;ix < this->pbasis;ix++) trho.up[ix] = rhoout_up[ix];
    for(int ix=0;ix < this->pbasis;ix++) trho.dw[ix] = rhoout_down[ix];

    double *grhof = new double[6*this->pbasis]();
    CToF_2d(this->pbasis, gx_up, grhof);
    CToF_2d(this->pbasis, gx_down, grhof + 3*this->pbasis);
    xc_gcx(&this->pbasis, &itwo, trho.data(), grhof, ex.data(), ec.data(),
            v1x.data(), v2x.data(), v1c.data(), v2c.data(), v2cud, &gargs);
    delete [] grhof;

    RmgTimer *RT4 = new RmgTimer("5-Functional: libxc");
///#pragma omp parallel for reduction(+:etxcgc,vtxcgc)
    for(int k=0;k < this->pbasis;k++) {
        double arho_up = fabs(rhoout_up[k]);
        double arho_down = fabs(rhoout_down[k]);
        double arho = arho_up + arho_down;
        double grho2[2];
        grho2[0] = gx_up[k]*gx_up[k] + gy_up[k]*gy_up[k] + gz_up[k]*gz_up[k];
        grho2[1] = gx_down[k]*gx_down[k] + gy_down[k]*gy_down[k] + gz_down[k]*gz_down[k];

//        if((arho_up > epsr) && (arho_down > epsr) && (grho2[0] > epsg) && (grho2[1] > epsg))
        {

            // first term of the gradient correction : D(rho*Exc)/D(rho)
            v_up[k] += (v1x.up[k] + v1c.up[k]);
            v_down[k] += ( v1x.dw[k] + v1c.dw[k]);

            vtxcgc += (v1x.up[k] + v1c.up[k]) * ( rhoout_up[k] - 0.5*rho_core[k]);
            vtxcgc += (v1x.dw[k] + v1c.dw[k]) * ( rhoout_down[k] - 0.5*rho_core[k]);
            etxcgc += (ex[k] + ec[k]);

            //  used later for second term of the gradient correction
            vxc2_up[k] = (v2x.up[k] + v2c.up[k]);
            vxc2_down[k] = (v2x.dw[k] + v2c.dw[k]);
        }

    }

    delete [] rhoout_down;
    delete [] rhoout_up;

    double *h = new double[6*this->pbasis]();
#if 1

    double *hx_up = h;
    double *hy_up = h + this->pbasis;
    double *hz_up = h + 2* this->pbasis;
    double *hx_dw = h + 3* this->pbasis;
    double *hy_dw = h + 4* this->pbasis;
    double *hz_dw = h + 5* this->pbasis;

#pragma omp parallel for 
    for(int k=0;k < this->pbasis;k++) {
        hx_up[k] = vxc2_up[k] * gx_up[k] + v2cud[k] * gx_down[k];
        hy_up[k] = vxc2_up[k] * gy_up[k] + v2cud[k] * gy_down[k];
        hz_up[k] = vxc2_up[k] * gz_up[k] + v2cud[k] * gz_down[k];
        hx_dw[k] = vxc2_down[k] * gx_down[k] + v2cud[k] * gx_up[k];
        hy_dw[k] = vxc2_down[k] * gy_down[k] + v2cud[k] * gy_up[k];
        hz_dw[k] = vxc2_down[k] * gz_down[k] + v2cud[k] * gz_up[k];
    }

    // second term of the gradient correction
    RmgTimer *RT5 = new RmgTimer("5-Functional: apply gradient");
    ApplyGradient (hx_up, gx_up, gy_up, gz_up, fd_order, "Fine");
    ApplyGradient (hx_dw, gx_down, gy_down, gz_down, fd_order, "Fine");
#pragma omp parallel for 
    for(int k=0;k < this->pbasis;k++) {
        v_up[k] -= gx_up[k];
        v_down[k] -= gx_down[k];
    }


    ApplyGradient (hy_up, gx_up, gy_up, gz_up, fd_order, "Fine");
    ApplyGradient (hy_dw, gx_down, gy_down, gz_down, fd_order, "Fine");
#pragma omp parallel for 
    for(int k=0;k < this->pbasis;k++) {
        v_up[k] -= gy_up[k];
        v_down[k] -= gy_down[k];
    }

    ApplyGradient (hz_up, gx_up, gy_up, gz_up, fd_order, "Fine");
    ApplyGradient (hz_dw, gx_down, gy_down, gz_down, fd_order, "Fine");
#pragma omp parallel for 
    for(int k=0;k < this->pbasis;k++) {
        v_up[k] -= gz_up[k];
        v_down[k] -= gz_down[k];
    }
    delete RT5;

#endif
    vtxc = vtxc + vtxcgc;
    etxc = etxc + etxcgc;
    delete RT4;

    delete [] h;
    delete [] grho_up;
}

void Functional::stress_vdW_DF (double *rho, double *rho_core, int nspin, double *sigma)
{
   int inlc = get_inlc();
        
    // No non-local correction just return
    if(inlc == 0) return;
    
    // inlc == 1 corresponds to vdW-DF1 and is the only one programmed currently
    if(inlc == 1) {
        double etxc, vtxc;
        double *v = new double[this->pbasis];
        Vdw *vdw = new Vdw (*this->Grid, *this->L, *this->T, 1, rho, rho_core, etxc, vtxc, v, this->gammaflag);
        vdw->stress_vdW_DF(rho, rho_core, nspin, sigma);
        delete vdw;
    }
    else { 
        throw RmgFatalException() << "Non-local correlation correction type not programmed" << " in " << __FILE__ << " at line " << __LINE__ << "\n";
    
    }

}

extern "C" const char *c_get_dft_name(void)
{
    return Functional::saved_dft_name.c_str();
}
