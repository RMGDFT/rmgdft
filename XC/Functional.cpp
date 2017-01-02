
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

extern "C" void __funct_MOD_set_dft_from_name( const char *name, std::size_t len );
extern "C" char *__funct_MOD_get_dft_name(void);
extern "C" bool __funct_MOD_dft_is_gradient(void);
extern "C" bool __funct_MOD_dft_is_meta(void);
extern "C" bool __funct_MOD_dft_is_hybrid(void);
extern "C" bool __funct_MOD_igcc_is_lyp(void);
extern "C" bool __funct_MOD_dft_has_finite_size_correction(void);
extern "C" bool __funct_MOD_dft_is_nonlocc(void);
extern "C" void __funct_MOD_xc (double *rho, double *ex, double *ec, double *vx, double *vc);
extern "C" void __funct_MOD_xc_spin (double *rho, double *zeta, double *ex, double *ec, double *vxup, double *vxdw, double *vcup, double *vcdw);

extern "C" void __funct_MOD_nlc (double *rho_valence, double *rho_core, int *nspin, double *ec, double *vx, double *vc);
extern "C" void __funct_MOD_gcxc (double *rho, double *grho, double *sx, double *sc, 
                                  double *v1x, double *v2x, double *v1c, double *v2c);
extern "C" void __funct_MOD_gcx_spin(double *rhoup, double *rhodown, double *grhoup, double *grhodown, double *sx,
                                  double *v1xup, double *v1xdw, double *v2xup, double *v2xdw);
extern "C" void __funct_MOD_gcc_spin_more( double *arho_up, double *arho_down,  double *grhoup, double *grhodw, double *grhoud,
                                  double *sc, double *v1cup, double *v1cdw, double *v2cup, double *v2cdw, double *v2cud );
extern "C" void __funct_MOD_gcc_spin( double *arho, double *zeta, double *grh2, double *sc, double *v1cup, double *v1cdw, double *v2c );

extern "C" int __funct_MOD_get_inlc(void);

bool Functional::dft_set=false;
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

    this->pbasis = G.get_P0_BASIS(G.default_FG_RATIO);
    this->N = G.get_NX_GRID(G.default_FG_RATIO) *
              G.get_NY_GRID(G.default_FG_RATIO) *
              G.get_NZ_GRID(G.default_FG_RATIO);

}




Functional::~Functional(void)
{


}


void Functional::set_dft_from_name(char *newdft_name) 
{
    if(!this->dft_set) {
        __funct_MOD_set_dft_from_name(newdft_name, std::strlen(newdft_name) );
        saved_dft_name = newdft_name;
    }
    else {
        std::cout << "Warning! You have attempted to reset the dft type which is a programming error. Ignoring." << std::endl;
    }
    this->dft_set = true;
}

const char *Functional::get_dft_name(void)
{
    if(!dft_set) {
        throw RmgFatalException() << "Error! get_dft_name called before dft type was set." << " in " << __FILE__ << " at line " << __LINE__ << "\n";
    }
    return saved_dft_name.c_str();
}

void Functional::set_dft_from_name(std::string newdft_name)
{
    if(!this->dft_set) {
        __funct_MOD_set_dft_from_name(newdft_name.c_str(), std::strlen(newdft_name.c_str()) );
        saved_dft_name = newdft_name;
    }
    else {
        std::cout << "Warning! You have attempted to reset the dft type which is a programming error. Ignoring." << std::endl;
    }
    this->dft_set = true;
}

bool Functional::dft_is_gradient(void)
{
    return __funct_MOD_dft_is_gradient();
}

bool Functional::dft_is_meta(void)
{
    return __funct_MOD_dft_is_meta();
}

bool Functional::dft_is_hybrid(void)
{
    return __funct_MOD_dft_is_hybrid();
}

bool Functional::igcc_is_lyp(void)
{
    return __funct_MOD_igcc_is_lyp();
}

bool Functional::dft_has_finite_size_correction(void)
{
    return __funct_MOD_dft_has_finite_size_correction();
}

bool Functional::dft_is_nonlocc(void)
{
    return __funct_MOD_dft_is_nonlocc();
}

void Functional::v_xc(double *rho, double *rho_core, double &etxc, double &vtxc, double *v, int spinflag)
{

   RmgTimer RT0("5-Functional");
   RmgTimer RT1("5-Functional: vxc");
   double rhoneg[2]{0.0,0.0};
   double *rho_up, *rho_down;
   double *v_up, *v_down;

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

   etxc = 0.0;
   vtxc = 0.0;
   for(int ix = 0;ix < this->pbasis;ix++) v[ix] = 0.0;
   if(spinflag) for(int ix = 0;ix < this->pbasis;ix++) v[ix + this->pbasis] = 0.0;


   // First get the local exchange and correlation
   RmgTimer *RT2 = new RmgTimer("5-Functional: vxc local");
   if(!spinflag) {

       double etxcl=0.0, vtxcl=0.0, rhonegl=0.0;
#pragma omp parallel for reduction(+:etxcl,vtxcl), reduction(-:rhonegl)
       // spin unpolarized  
       for(int ix=0;ix < this->pbasis;ix++) {

           double trho = rho[ix] + rho_core[ix];
           double atrho = fabs(trho);
           double ex, ec, vx, vc;
           if(atrho > SMALL_CHARGE) {

               __funct_MOD_xc( &trho, &ex, &ec, &vx, &vc);
               v[ix] = vx + vc;
               etxcl = etxcl + ( ex + ec ) * trho;
               vtxcl = vtxcl + v[ix] * rho[ix];

           }

           else {
               double rhotem = SMALL_CHARGE * (1.0 + SMALL_CHARGE);
               __funct_MOD_xc( &rhotem, &ex, &ec, &vx, &vc );
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
               __funct_MOD_xc_spin( &trho, &zeta, &ex, &ec, &vx0, &vx1, &vc0, &vc1 );
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

   vtxc = vtxc * L->omega / (double)this->N;
   etxc = etxc * L->omega / (double)this->N;

   // Next add in any gradient corrections
   RmgTimer *RT3 = new RmgTimer("5-Functional: vxc grad");
   if(!spinflag) {
       this->gradcorr(rho, rho_core, etxc, vtxc, v);
   }
   else {
       this->gradcorr_spin(rho, rho_core, etxc, vtxc, v);
   }
   delete RT3;

   // And finally any non-local corrections
   RmgTimer *RT4 = new RmgTimer("5-Functional: vxc nonlocal");
   if(this->dft_is_nonlocc()) {
       this->nlc(rho, rho_core, etxc, vtxc, v, spinflag);
       //__funct_MOD_nlc( rho, rho_core, &nspin, &etxc, &vtxc, v );
   }
   delete RT4;

   vtxc = RmgSumAll(vtxc, this->T->get_MPI_comm());
   etxc = RmgSumAll(etxc, this->T->get_MPI_comm());
   //printf("GGGGGGGG  %20.12f  %20.12f\n",vtxc,etxc);

}

// Applies non-local corrections for the correlation
void Functional::nlc(double *rho, double *rho_core, double &etxc, double &vtxc, double *v, int spinflag)
{
    int inlc = __funct_MOD_get_inlc();

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

    if(!this->dft_is_gradient()) return;

    double etxcgc = 0.0;
    double vtxcgc = 0.0;
    double grho2[2];
    const double epsr=1.0e-10;
    const double epsg = 1.0e-16;
    double epsg_guard = sqrt(epsg);


    double *rhoout = new double[this->pbasis];
    double *grho = new double[3*this->pbasis];
    double *vxc2 = new double[this->pbasis]();
    double *d2rho = new double[this->pbasis];
    double *gx = grho;
    double *gy = gx + this->pbasis;
    double *gz = gy + this->pbasis;
    double *h = new double[3*this->pbasis]();

    // Get rho plus rhocore
    for(int ix=0;ix < this->pbasis;ix++) rhoout[ix] = rho[ix] + rho_core[ix];

    // calculate the gradient of rho + rho_core
    RmgTimer *RT2 = new RmgTimer("5-Functional: apply gradient");
    ApplyGradient (rhoout, gx, gy, gz, APP_CI_EIGHT, "Fine");
    delete RT2;


    // and the Laplacian
    RmgTimer *RT3 = new RmgTimer("5-Functional: apply laplacian");
    ApplyLaplacian (rhoout, d2rho, APP_CI_EIGHT, "Fine");
    delete RT3;



    RmgTimer *RT4 = new RmgTimer("5-Functional: libxc");

#pragma omp parallel for reduction(+:etxcgc,vtxcgc)
    for(int k=0;k < this->pbasis;k++) {

        double arho = fabs(rhoout[k]);
        if(arho > epsr) {

            grho2[0] = gx[k]*gx[k] + gy[k]*gy[k] + gz[k]*gz[k];

            if(grho2[0] > epsg) {
                double sx, sc, v1x, v2x, v1c, v2c;
                double segno = 1.0;
                if(rhoout[k] < 0.0)segno = -1.0; 

                double pgrho2 = grho2[0] + epsg_guard;
                __funct_MOD_gcxc( &arho, &pgrho2, &sx, &sc, &v1x, &v2x, &v1c, &v2c );
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
    ApplyGradient (vxc2, h, &h[this->pbasis], &h[2*this->pbasis], APP_CI_EIGHT, "Fine");
    delete RT5;

#pragma omp parallel for
    for(int ix=0;ix < this->pbasis;ix++) {

        v[ix] -= ( h[ix] * gx[ix] +
                h[ix+this->pbasis] * gy[ix] + 
                h[ix+2*this->pbasis] * gz[ix] ) ;
        v[ix] -= vxc2[ix] * d2rho[ix];

    }

    //printf("VTXC1 = %18.12f  ETXC1 = %18.12f\n", 
    //2.0*L->omega*RmgSumAll(vtxcgc, this->T->get_MPI_comm())/(double)this->N,
    //2.0*L->omega*RmgSumAll(etxcgc, this->T->get_MPI_comm())/(double)this->N);

    vtxc = vtxc + L->omega * vtxcgc / (double)this->N;
    etxc = etxc + L->omega * etxcgc / (double)this->N;


    delete [] h;
    delete [] vxc2;
    delete [] d2rho;
    delete [] grho;
    delete [] rhoout;


}

// Applies gradient corrections for spin case
void Functional::gradcorr_spin(double *rho, double *rho_core, double &etxc, double &vtxc, double *v)
{
    if(!this->dft_is_gradient()) return;

    double *rho_up, *rho_down;
    double *v_up, *v_down;

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

    double etxcgc = 0.0;
    double vtxcgc = 0.0;
 
    double grho2[2];
    const double epsr=1.0e-10;
    const double epsg = 1.0e-16;
    double epsg_guard = sqrt(epsg);

    double *grho_up = new double[3*this->pbasis];
    double *grho_down = new double[3*this->pbasis];
    double *d2rho_up = new double[this->pbasis];
    double *d2rho_down = new double[this->pbasis];
    double *vxc2_up = new double[this->pbasis]();
    double *vxc2_down = new double[this->pbasis]();
    double *v2cud = new double[this->pbasis]();
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
    ApplyGradient (rhoout_up, gx_up, gy_up, gz_up, APP_CI_EIGHT, "Fine");
    ApplyGradient (rhoout_down, gx_down, gy_down, gz_down, APP_CI_EIGHT, "Fine");
    delete RT2;

    // and the Laplacian
    RmgTimer *RT3 = new RmgTimer("5-Functional: apply laplacian");
    ApplyLaplacian (rhoout_up, d2rho_up, APP_CI_EIGHT, "Fine");
    ApplyLaplacian (rhoout_down, d2rho_down, APP_CI_EIGHT, "Fine");
    delete RT3;


    RmgTimer *RT4 = new RmgTimer("5-Functional: libxc");
#pragma omp parallel for reduction(+:etxcgc,vtxcgc)
    for(int k=0;k < this->pbasis;k++) {
        double arho_up = fabs(rhoout_up[k]);
        double arho_down = fabs(rhoout_down[k]);
        double arho = arho_up + arho_down;

        grho2[0] = gx_up[k]*gx_up[k] + gy_up[k]*gy_up[k] + gz_up[k]*gz_up[k];
        grho2[1] = gx_down[k]*gx_down[k] + gy_down[k]*gy_down[k] + gz_down[k]*gz_down[k];

        double pgrho2_up = grho2[0] + epsg_guard;
        double pgrho2_down = grho2[1] + epsg_guard;
        double v1xup, v1xdw, v2xup, v2xdw, sx;

        __funct_MOD_gcx_spin( &arho_up, &arho_down, &pgrho2_up,
                        &pgrho2_down, &sx, &v1xup, &v1xdw, &v2xup, &v2xdw );

         double sc    = 0.0;
         double v1cup = 0.0;
         double v1cdw = 0.0;
         double v2c   = 0.0;
         double v2cup = 0.0;
         double v2cdw = 0.0;
         v2cud[k] = 0.0;

        if(arho > epsr) {

            if(__funct_MOD_igcc_is_lyp()) {

                double grhoup = gx_up[k]*gx_up[k] + gy_up[k]*gy_up[k] + gz_up[k]*gz_up[k];
                double grhodw = gx_down[k]*gx_down[k] + gy_down[k]*gy_down[k] + gz_down[k]*gz_down[k];
                double grhoud = gx_up[k]*gx_down[k] + gy_up[k]*gy_down[k] + gz_up[k]*gz_down[k];

                __funct_MOD_gcc_spin_more( &arho_up, &arho_down, &grhoup, &grhodw, &grhoud,
                                  &sc, &v1cup, &v1cdw, &v2cup, &v2cdw, &v2cud[k] );

            }
            else {

                double zeta = ( rhoout_up[k] - rhoout_down[k]) / arho;
                zeta = (arho_up - arho_down) / arho;
                double grh2 = (gx_up[k] + gx_down[k]) * (gx_up[k] + gx_down[k]) +
                              (gy_up[k] + gy_down[k]) * (gy_up[k] + gy_down[k]) +
                              (gz_up[k] + gz_down[k]) * (gz_up[k] + gz_down[k]);
                
                __funct_MOD_gcc_spin( &arho, &zeta, &grh2, &sc, &v1cup, &v1cdw, &v2c );
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

    double *h = new double[3*this->pbasis]();

    // second term of the gradient correction
    RmgTimer *RT5 = new RmgTimer("5-Functional: apply gradient");
    ApplyGradient (vxc2_up, h, &h[this->pbasis], &h[2*this->pbasis], APP_CI_EIGHT, "Fine");
    delete RT5;

#pragma omp parallel for 
    for(int ix=0;ix < this->pbasis;ix++) {

        v_up[ix] -= ( h[ix] * gx_up[ix] +
                h[ix+this->pbasis] * gy_up[ix] +
                h[ix+2*this->pbasis] * gz_up[ix] ) ;

        v_up[ix] -= vxc2_up[ix] * d2rho_up[ix];

    }

    RmgTimer *RT6 = new RmgTimer("5-Functional: apply gradient");
    ApplyGradient (vxc2_down, h, &h[this->pbasis], &h[2*this->pbasis], APP_CI_EIGHT, "Fine");
    delete RT6;

#pragma omp parallel for 
    for(int ix=0;ix < this->pbasis;ix++) {

        v_down[ix] -= ( h[ix] * gx_down[ix] +
                h[ix+this->pbasis] * gy_down[ix] +
                h[ix+2*this->pbasis] * gz_down[ix] ) ;

        v_down[ix] -= vxc2_down[ix] * d2rho_down[ix];

    }

    ApplyGradient (v2cud, h, &h[this->pbasis], &h[2*this->pbasis], APP_CI_EIGHT, "Fine");
#pragma omp parallel for 
    for(int ix=0;ix < this->pbasis;ix++) {
        v_up[ix] -= ( h[ix] * gx_down[ix] +
                h[ix+this->pbasis] * gy_down[ix] +
                h[ix+2*this->pbasis] * gz_down[ix] ) ;
        v_up[ix] -= v2cud[ix] * d2rho_down[ix];
    }
#pragma omp parallel for 
    for(int ix=0;ix < this->pbasis;ix++) {
        v_down[ix] -= ( h[ix] * gx_up[ix] +
                h[ix+this->pbasis] * gy_up[ix] +
                h[ix+2*this->pbasis] * gz_up[ix] ) ;
        v_down[ix] -= v2cud[ix] * d2rho_up[ix];
    }

    vtxc = vtxc + L->omega * vtxcgc / (double)this->N;
    etxc = etxc + L->omega * etxcgc / (double)this->N;


    delete RT4;

    delete [] h;
    delete [] v2cud;
    delete [] vxc2_down;
    delete [] vxc2_up;
    delete [] d2rho_up;
    delete [] d2rho_down;
    delete [] grho_down;
    delete [] grho_up;
}


extern "C" const char *c_get_dft_name(void)
{
    return Functional::saved_dft_name.c_str();
}
