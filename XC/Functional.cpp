
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
#define start_exx               RMG_FC_MODULE(funct,start_exx,mod_FUNCT,START_EXX)
#define stop_exx                RMG_FC_MODULE(funct,stop_exx,mod_FUNCT,STOP_EXX)
#define dft_is_gradient         RMG_FC_MODULE(funct,dft_is_gradient,mod_FUNCT,DFT_IS_GRADIENT)
#define dft_is_meta             RMG_FC_MODULE(funct,dft_is_meta,mod_FUNCT,DFT_IS_META)
#define dft_is_hybrid           RMG_FC_MODULE(funct,dft_is_hybrid,mod_FUNCT,DFT_IS_HYBRID)
#define get_exx_fraction        RMG_FC_MODULE(funct,get_exx_fraction,mod_FUNCT,GET_EXX_FRACTION)
#define set_exx_fraction        RMG_FC_MODULE(funct,set_exx_fraction,mod_FUNCT,SET_EXX_FRACTION)
#define igcc_is_lyp             RMG_FC_MODULE(funct,igcc_is_lyp,mod_FUNCT,IGCC_IS_LYP)
#define dft_has_finite_size_correction RMG_FC_MODULE(funct,dft_has_finite_size_correction,mod_FUNCT,DFT_HAS_FINITE_SIZE_CORRECTION)
#define dft_is_nonlocc          RMG_FC_MODULE(funct,dft_is_nonlocc,mod_FUNCT,DFT_IS_NONLOCC)
#define xc                      RMG_FC_MODULE(funct,xc,mod_FUNCT,XC)
#define xc_spin                 RMG_FC_MODULE(funct,xc_spin,mod_FUNCT,XC_SPIN)
#define nlc                     RMG_FC_MODULE(funct,nlc,mod_FUNCT,NLC)
#define gcxc                    RMG_FC_MODULE(funct,gcxc,mod_FUNCT,GCXC)
#define gcx_spin                RMG_FC_MODULE(funct,gcx_spin,mod_FUNCT,GCX_SPIN)
#define gcc_spin_more           RMG_FC_MODULE(funct,gcc_spin_more,mod_FUNCT,GCC_SPIN_MORE)
#define gcc_spin                RMG_FC_MODULE(funct,gcc_spin,mod_FUNCT,GCC_SPIN)
#define get_inlc                RMG_FC_MODULE(funct,get_inlc,mod_FUNCT,GET_INLC)
#define get_gau_parameter       RMG_FC_MODULE(funct,get_gau_parameter,mod_FUNCT,GET_GAU_PARAMETER)
#define set_gau_parameter       RMG_FC_MODULE(funct,set_gau_parameter,mod_FUNCT,SET_GAU_PARAMETER)
#define get_screening_parameter       RMG_FC_MODULE(funct,get_screening_parameter,mod_FUNCT,GET_SCREENING_PARAMETER)
#define set_screening_parameter       RMG_FC_MODULE(funct,set_screening_parameter,mod_FUNCT,SET_SCREENING_PARAMETER)

extern "C" void set_dft_from_name( const char *name, std::size_t len );
extern "C" char *get_dft_name(void);
extern "C" void start_exx(void);
extern "C" void stop_exx(void);
extern "C" bool dft_is_gradient(void);
extern "C" bool dft_is_meta(void);
extern "C" bool dft_is_hybrid(void);
extern "C" double get_exx_fraction(void);
extern "C" void set_exx_fraction(double *frac);
extern "C" bool igcc_is_lyp(void);
extern "C" bool dft_has_finite_size_correction(void);
extern "C" bool dft_is_nonlocc(void);
extern "C" void xc (double *rho, double *ex, double *ec, double *vx, double *vc);
extern "C" void xc_spin (double *rho, double *zeta, double *ex, double *ec, double *vxup, double *vxdw, double *vcup, double *vcdw);

extern "C" void nlc (double *rho_valence, double *rho_core, int *nspin, double *ec, double *vx, double *vc);
extern "C" void gcxc (double *rho, double *grho, double *sx, double *sc, 
                                  double *v1x, double *v2x, double *v1c, double *v2c);
extern "C" void gcx_spin(double *rhoup, double *rhodown, double *grhoup, double *grhodown, double *sx,
                                  double *v1xup, double *v1xdw, double *v2xup, double *v2xdw);
extern "C" void gcc_spin_more( double *arho_up, double *arho_down,  double *grhoup, double *grhodw, double *grhoud,
                                  double *sc, double *v1cup, double *v1cdw, double *v2cup, double *v2cdw, double *v2cud );
extern "C" void gcc_spin( double *arho, double *zeta, double *grh2, double *sc, double *v1cup, double *v1cdw, double *v2c );

extern "C" int get_inlc(void);
extern "C" double get_gau_parameter(void);
extern "C" void set_gau_parameter(double *);
extern "C" double get_screening_parameter(void);
extern "C" void set_screening_parameter(double *);

bool Functional::dft_set=false;
bool Functional::exx_started=false;
std::string Functional::saved_dft_name;

double SMALL_CHARGE=1.0e-10;

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
    this->Grid = &G;
    this->T = &T;
    this->L = &L;
    this->gammaflag = gamma_flag;

    this->fd_order = ct.kohn_sham_fd_order;
    this->pbasis = G.get_P0_BASIS(G.default_FG_RATIO);
    this->N = G.get_NX_GRID(G.default_FG_RATIO) *
              G.get_NY_GRID(G.default_FG_RATIO) *
              G.get_NZ_GRID(G.default_FG_RATIO);

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

void Functional::v_xc(double *rho_in, double *rho_core, double &etxc, double &vtxc, double *v, int nspin)
{

    RmgTimer RT0("5-Functional");
    RmgTimer RT1("5-Functional: vxc");
    const double epsr=1.0e-6;
    double rhoneg[2]{0.0,0.0};
    double *rho_up=NULL, *rho_down=NULL;
    double *v_up=NULL, *v_down=NULL;
    double *rho = new double[nspin*this->pbasis];

    for(int ix=0;ix < this->pbasis;ix++)rho[ix] = rho_in[ix];
    // for collinear case, spin up and down are in different processor groups.
    if(nspin==2)
    {
        for(int ix=0;ix < this->pbasis;ix++)rho[ix+this->pbasis] = rho_in[ix+this->pbasis];

        if(pct.spinpe == 0) {
            rho_up = rho;
            rho_down = &rho[this->pbasis];
            v_up = v;
            v_down = &v[this->pbasis];
        }
        else {
            rho_down = rho;
            rho_up = &rho[this->pbasis];
            v_down = v;
            v_up = &v[this->pbasis];
        }
    }
    else if(nspin == 4)
    {
        rho_up = rho;
        rho_down = &rho[this->pbasis];
        v_up = v;
        v_down = &v[this->pbasis];
        double mrho;
        
        for(int idx = 0; idx < this->pbasis; idx++)
        {
            mrho = rho_in[idx + this->pbasis] *rho_in[idx + this->pbasis];
            mrho += rho_in[idx + 2*this->pbasis] *rho_in[idx + 2*this->pbasis];
            mrho += rho_in[idx + 3*this->pbasis] *rho_in[idx + 3*this->pbasis];
            mrho = std::sqrt(mrho);
            rho_up[idx] = 0.5 * (rho_in[idx] + mrho);
            rho_down[idx] = 0.5* (rho_in[idx] - mrho);
        }

    }

    etxc = 0.0;
    vtxc = 0.0;
    for(int ix = 0;ix < this->pbasis;ix++) v[ix] = 0.0;
    for(int ispin = 1; ispin < nspin; ispin++)
        for(int ix = 0;ix < this->pbasis;ix++) v[ix + ispin * this->pbasis] = 0.0;


    // First get the local exchange and correlation
    RmgTimer *RT2 = new RmgTimer("5-Functional: vxc local");
    if(nspin==1) {

        double etxcl=0.0, vtxcl=0.0, rhonegl=0.0;
#pragma omp parallel for reduction(+:etxcl), reduction(+:vtxcl), reduction(-:rhonegl)
        // spin unpolarized  
        for(int ix=0;ix < this->pbasis;ix++) {

            double trho = rho[ix] + rho_core[ix];
            double atrho = fabs(trho);
            double ex, ec, vx, vc;
            if(atrho > SMALL_CHARGE) {

                xc( &trho, &ex, &ec, &vx, &vc);
                v[ix] = vx + vc;
                etxcl = etxcl + ( ex + ec ) * trho;
                vtxcl = vtxcl + v[ix] * rho[ix];

            }

            else {
                double rhotem = SMALL_CHARGE * (1.0 + SMALL_CHARGE);
                xc( &rhotem, &ex, &ec, &vx, &vc );
                double frac = std::cbrt(atrho/SMALL_CHARGE);
                v[ix] = (vx + vc) * frac;
                etxcl = etxcl + ( ex + ec ) * trho * frac;
                vtxcl = vtxcl + v[ix] * rho[ix];

            }

            if(rho[ix] < 0.0) rhonegl = rhonegl - rho[ix];

        } 
        etxc += etxcl;
        vtxc += vtxcl;
        rhoneg[0] += rhonegl;

    } 
    else {

        // spin polarized
        double etxcl=0.0, vtxcl=0.0;
#pragma omp parallel for reduction(+:etxcl,vtxcl)
        for(int ix=0;ix < this->pbasis;ix++) {

            double trho = rho_up[ix] + rho_down[ix] + rho_core[ix];
            double atrho = fabs(trho);
            double vx0, vx1, vc0, vc1;
            double ex, ec;
            if(atrho > SMALL_CHARGE) {

                double zeta = (rho_up[ix] - rho_down[ix]) / atrho;
                if( fabs( zeta ) > 1.0 ) {
                    double tzeta = 1.0;
                    if(zeta < 0.0) tzeta = -1.0;
                    zeta = tzeta;
                }
                xc_spin( &trho, &zeta, &ex, &ec, &vx0, &vx1, &vc0, &vc1 );
                v_up[ix] = vx0 + vc0;
                v_down[ix] = vx1 + vc1;
                etxcl = etxcl + ( ex + ec ) * trho;
                vtxcl = vtxcl + v_up[ix] * rho_up[ix] + v_down[ix] * rho_down[ix];

            }
            else {

            }

        }
        etxc += etxcl;
        vtxc += vtxcl;

    }
    delete RT2;


    // Next add in any gradient corrections
    RmgTimer *RT3 = new RmgTimer("5-Functional: vxc grad");
    if(nspin == 1) {
        this->gradcorr(rho, rho_core, etxc, vtxc, v);
    }
    else {
        this->gradcorr_spin(rho_up, rho_down, rho_core, etxc, vtxc, v_up, v_down);
    }
    delete RT3;

    if(nspin == 4)
    {
        double vt, vd, mrho;
        for(int idx = 0; idx < this->pbasis; idx++)
        {
            vt = 0.5 * (v_up[idx] + v_down[idx]);
            vd = 0.5 * (v_up[idx] - v_down[idx]);

            mrho = rho_in[idx + this->pbasis] *rho_in[idx + this->pbasis];
            mrho += rho_in[idx + 2*this->pbasis] *rho_in[idx + 2*this->pbasis];
            mrho += rho_in[idx + 3*this->pbasis] *rho_in[idx + 3*this->pbasis];
            mrho = std::sqrt(mrho);
            
            v[idx] = vt;
            
            if(mrho > epsr)
            {
                v[idx +   this->pbasis] = rho_in[idx +   this->pbasis] /mrho * vd;
                v[idx + 2*this->pbasis] = rho_in[idx + 2*this->pbasis] /mrho * vd;
                v[idx + 3*this->pbasis] = rho_in[idx + 3*this->pbasis] /mrho * vd;
            }
            else
            {
                v[idx +   this->pbasis] = 0.0;
                v[idx + 2*this->pbasis] = 0.0;
                v[idx + 3*this->pbasis] = 0.0;
            }
        }    
    }
    // And finally any non-local corrections
    RmgTimer *RT4 = new RmgTimer("5-Functional: vxc nonlocal");
    if(this->dft_is_nonlocc_rmg()) {
        if(nspin == 4) 
        {
            throw RmgFatalException() << "vdw with noncollinear not programed. " << " in " << __FILE__ << " at line " << __LINE__ << "\n";

        }
        double netxc=0.0, nvtxc=0.0;
        this->nlc_rmg(rho, rho_core, netxc, nvtxc, v);
        vtxc += nvtxc;
        etxc += netxc;
    }
    delete RT4;

    vtxc = vtxc * L->omega / (double)this->N;
    etxc = etxc * L->omega / (double)this->N;

    vtxc = RmgSumAll(vtxc, this->T->get_MPI_comm());
    etxc = RmgSumAll(etxc, this->T->get_MPI_comm());

    delete [] rho;
    if(Rmg_G->default_FG_RATIO > 1)
    {
        for(int is = 0; is < nspin; is++)
            FftFilter(&v[is*pbasis], *fine_pwaves, *coarse_pwaves, LOW_PASS);
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

    double etxcgc = 0.0;
    double vtxcgc = 0.0;
    const double epsr=1.0e-6;
    const double epsg = 1.0e-10;
    double epsg_guard = ct.epsg_guard;

    double *rhoout = new double[this->pbasis];
    double *grho = new double[3*this->pbasis];
    double *vxc2 = this->vxc2;
    double *d2rho = new double[this->pbasis];
    double *gx = grho;
    double *gy = gx + this->pbasis;
    double *gz = gy + this->pbasis;
    double *h = new double[3*this->pbasis]();

    // Get rho plus rhocore
    for(int ix=0;ix < this->pbasis;ix++) rhoout[ix] = rho[ix] + rho_core[ix];

    // calculate the gradient of rho + rho_core
    RmgTimer *RT2 = new RmgTimer("5-Functional: apply gradient");
    ApplyGradient (rhoout, gx, gy, gz, fd_order, "Fine");
    //FftGradientFine(rhoout, gx, gy, gz);
    delete RT2;


    // and the Laplacian
    RmgTimer *RT3 = new RmgTimer("5-Functional: apply laplacian");
    //FftLaplacianFine(rhoout, d2rho);
    ApplyLaplacian (rhoout, d2rho, fd_order, "Fine");

    delete RT3;



    RmgTimer *RT4 = new RmgTimer("5-Functional: libxc");

#pragma omp parallel for reduction(+:etxcgc,vtxcgc)
    for(int k=0;k < this->pbasis;k++) {

        double arho = fabs(rhoout[k]);
        double grho2[2];
        if(arho > epsr) {

            grho2[0] = gx[k]*gx[k] + gy[k]*gy[k] + gz[k]*gz[k];

            if(grho2[0] > epsg) {
                double sx, sc, v1x, v2x, v1c, v2c;
                double segno = 1.0;
                if(rhoout[k] < 0.0)segno = -1.0; 

                double pgrho2 = grho2[0] + epsg_guard;
                gcxc( &arho, &pgrho2, &sx, &sc, &v1x, &v2x, &v1c, &v2c );
                //
                // first term of the gradient correction : D(rho*Exc)/D(rho)
                v[k] = v[k] +( v1x + v1c );
                //
                //  used later for second term of the gradient correction
                vxc2[k] = ( v2x + v2c );
                // 
                vtxcgc = vtxcgc+( v1x + v1c ) * ( rhoout[k] - rho_core[k] );
                etxcgc = etxcgc+( sx + sc ) * segno;

            }
        }
    }
    delete RT4;

    for(int ix=0;ix < this->pbasis;ix++) rhoout[ix] = rhoout[ix] - rho_core[ix];

    // 
    // ... second term of the gradient correction :
    // ... \sum_alpha (D / D r_alpha) ( D(rho*Exc)/D(grad_alpha rho) )
    // 

    RmgTimer *RT5 = new RmgTimer("5-Functional: apply gradient");
    ApplyGradient (vxc2, h, &h[this->pbasis], &h[2*this->pbasis], fd_order, "Fine");
    delete RT5;

    double vtxcgc_1 = 0.0;
#pragma omp parallel for reduction(+:vtxcgc_1)
    for(int ix=0;ix < this->pbasis;ix++) {
        double arho = fabs(rhoout[ix]);
        if(arho > epsr) {
            double gdot =  ( h[ix] * gx[ix] +
                    h[ix+this->pbasis] * gy[ix] + 
                    h[ix+2*this->pbasis] * gz[ix] ) ;
            v[ix] -= gdot;
            v[ix] -= vxc2[ix] * d2rho[ix];
            vtxcgc_1 -= rhoout[ix]*(gdot + vxc2[ix] * d2rho[ix]);
        }
    }

    vtxc = vtxc + vtxcgc + vtxcgc_1;
    etxc = etxc + etxcgc;

    delete [] h;
    delete [] d2rho;
    delete [] grho;
    delete [] rhoout;


}

// Applies gradient corrections for spin case
void Functional::gradcorr_spin(double *rho_up, double *rho_down, double *rho_core, double &etxc, double &vtxc, double *v_up, double *v_down)
{
    if(!this->dft_is_gradient_rmg()) return;

    double etxcgc = 0.0;
    double vtxcgc = 0.0;

    const double epsr=1.0e-6;
    const double epsg = 1.0e-10;
    double epsg_guard = ct.epsg_guard;

    double *grho_up = new double[3*this->pbasis];
    double *grho_down = new double[3*this->pbasis];
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


    RmgTimer *RT4 = new RmgTimer("5-Functional: libxc");
#pragma omp parallel for reduction(+:etxcgc,vtxcgc)
    for(int k=0;k < this->pbasis;k++) {
        double arho_up = fabs(rhoout_up[k]);
        double arho_down = fabs(rhoout_down[k]);
        double arho = arho_up + arho_down;
        double grho2[2];
        grho2[0] = gx_up[k]*gx_up[k] + gy_up[k]*gy_up[k] + gz_up[k]*gz_up[k];
        grho2[1] = gx_down[k]*gx_down[k] + gy_down[k]*gy_down[k] + gz_down[k]*gz_down[k];

        double pgrho2_up = grho2[0] + epsg_guard;
        double pgrho2_down = grho2[1] + epsg_guard;
        double v1xup, v1xdw, v2xup, v2xdw, sx;

        gcx_spin( &arho_up, &arho_down, &pgrho2_up,
                &pgrho2_down, &sx, &v1xup, &v1xdw, &v2xup, &v2xdw );

        double sc    = 0.0;
        double v1cup = 0.0;
        double v1cdw = 0.0;
        double v2c   = 0.0;
        double v2cup = 0.0;
        double v2cdw = 0.0;
        v2cud[k] = 0.0;

        if(arho > epsr && grho2[0] > epsg && grho2[1] > epsg) {

            if(igcc_is_lyp()) {

                double grhoup = gx_up[k]*gx_up[k] + gy_up[k]*gy_up[k] + gz_up[k]*gz_up[k];
                double grhodw = gx_down[k]*gx_down[k] + gy_down[k]*gy_down[k] + gz_down[k]*gz_down[k];
                double grhoud = gx_up[k]*gx_down[k] + gy_up[k]*gy_down[k] + gz_up[k]*gz_down[k];

                gcc_spin_more( &arho_up, &arho_down, &grhoup, &grhodw, &grhoud,
                        &sc, &v1cup, &v1cdw, &v2cup, &v2cdw, &v2cud[k] );

            }
            else {

                double zeta = ( rhoout_up[k] - rhoout_down[k]) / arho;
                zeta = (arho_up - arho_down) / arho;
                double grh2 = (gx_up[k] + gx_down[k]) * (gx_up[k] + gx_down[k]) +
                    (gy_up[k] + gy_down[k]) * (gy_up[k] + gy_down[k]) +
                    (gz_up[k] + gz_down[k]) * (gz_up[k] + gz_down[k]);

                grh2 += epsg_guard;
                gcc_spin( &arho, &zeta, &grh2, &sc, &v1cup, &v1cdw, &v2c );
                v2cup = v2c;
                v2cdw = v2c;
                v2cud[k] = v2c;


                // first term of the gradient correction : D(rho*Exc)/D(rho)
                v_up[k] = v_up[k] + ( v1xup + v1cup );
                v_down[k] = v_down[k] + ( v1xdw + v1cdw );

                vtxcgc = vtxcgc +
                    ( v1xup + v1cup ) * ( rhoout_up[k] - 0.5*rho_core[k]);
                vtxcgc = vtxcgc + 
                    ( v1xdw + v1cdw ) * ( rhoout_down[k] - 0.5*rho_core[k]);
                etxcgc = etxcgc + ( sx + sc );

                //  used later for second term of the gradient correction
                vxc2_up[k] = ( v2xup + v2cup );
                vxc2_down[k] = ( v2xdw + v2cdw );

            }

        }

    }

    delete [] rhoout_down;
    delete [] rhoout_up;

    double *h = new double[6*this->pbasis]();

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

    vtxc = vtxc + vtxcgc;
    etxc = etxc + etxcgc;

    delete RT4;

    delete [] h;
    delete [] grho_down;
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
