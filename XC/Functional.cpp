
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
#include "FiniteDiff.h"
#include "xc.h"
#include "RmgSumAll.h"
#include "transition.h"

extern "C" void __funct_MOD_set_dft_from_name( const char *name, std::size_t len );
extern "C" char *__funct_MOD_get_dft_name(void);
extern "C" bool __funct_MOD_dft_is_gradient(void);
extern "C" bool __funct_MOD_dft_is_meta(void);
extern "C" bool __funct_MOD_dft_is_hybrid(void);
extern "C" bool __funct_MOD_igcc_is_lyp(void);
extern "C" bool __funct_MOD_dft_has_finite_size_correction(void);
extern "C" bool __funct_MOD_dft_is_nonlocc(void);
extern "C" void __funct_MOD_xc (double *rho, double *ex, double *ec, double *vx, double *vc);
extern "C" void __funct_MOD_nlc (double *rho_valence, double *rho_core, int *nspin, double *ec, double *vx, double *vc);
extern "C" void __funct_MOD_gcxc (double *rho, double *grho, double *sx, double *sc, 
                                  double *v1x, double *v2x, double *v1c, double *v2c);
extern "C" int __funct_MOD_get_inlc(void);

bool Functional::dft_set=false;
std::string Functional::saved_dft_name;

#define SMALL_CHARGE 1.0e-10

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
    this->Grid = &G;
    this->T = &T;
    this->L = &L;
    this->gammaflag = gamma_flag;

    this->pbasis = G.get_P0_BASIS(G.default_FG_RATIO);
    this->hxgrid = G.get_hxgrid(G.default_FG_RATIO);
    this->hygrid = G.get_hygrid(G.default_FG_RATIO);
    this->hzgrid = G.get_hzgrid(G.default_FG_RATIO);
    this->dimx = G.get_PX0_GRID(G.default_FG_RATIO);
    this->dimy = G.get_PY0_GRID(G.default_FG_RATIO);
    this->dimz = G.get_PZ0_GRID(G.default_FG_RATIO);
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

   double vx[2]{0.0,0.0}, vc[2]{0.0,0.0}, rhoneg[2]{0.0,0.0};
   double ex=0.0, ec=0.0;

   etxc = 0.0;
   vtxc = 0.0;
   for(int ix = 0;ix < this->pbasis;ix++) v[ix] = 0.0;


   // First get the local exchange and correlation
   if(!spinflag) {

       // spin unpolarized  
       for(int ix=0;ix < this->pbasis;ix++) {

           double trho = rho[ix] + rho_core[ix];
           double atrho = fabs(trho);
           if(atrho > SMALL_CHARGE) {

               __funct_MOD_xc( &trho, &ex, &ec, &vx[0], &vc[0] );
               v[ix] = vx[0] + vc[0];
               etxc = etxc + ( ex + ec ) * trho;
               vtxc = vtxc + v[ix] * rho[ix];
 
           }

           if(rho[ix] < 0.0) rhoneg[0] = rhoneg[0] - rho[ix];

       } 
       
       
   } 

   vtxc = vtxc * L->omega / (double)this->N;
   etxc = etxc * L->omega / (double)this->N;

   // Next add in any gradient corrections
   this->gradcorr(rho, rho_core, etxc, vtxc, v );

   // And finally any non-local corrections
   if(this->dft_is_nonlocc()) {
       this->nlc(rho, rho_core, etxc, vtxc, v, spinflag);
       //__funct_MOD_nlc( rho, rho_core, &nspin, &etxc, &vtxc, v );
   }

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
    const double epsr=1.0e-6;
    const double epsg = 1.0e-10;


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
    CPP_app_grad_driver (this->L, this->T, rhoout, gx, gy, gz,
                         this->dimx, this->dimy, this->dimz, 
                         this->hxgrid, this->hygrid, this->hzgrid, 
                         APP_CI_EIGHT);


    // and the Laplacian
    //app6_del2(rhoout, d2rho, this->dimx, this->dimy, this->dimz, this->hxgrid, this->hygrid, this->hzgrid);
    CPP_app_del2_driver (this->L, this->T, rhoout, d2rho, this->dimx, this->dimy, this->dimz, this->hxgrid, this->hygrid, this->hzgrid, APP_CI_EIGHT);



    for(int k=0;k < this->pbasis;k++) {

        double arho = fabs(rhoout[k]);
        if(arho > epsr) {

            grho2[0] = gx[k]*gx[k] + gy[k]*gy[k] + gz[k]*gz[k];

            if(grho2[0] > epsg) {
                double sx, sc, v1x, v2x, v1c, v2c;
                double segno = 1.0;
                if(rhoout[k] < 0.0)segno = -1.0; 

                __funct_MOD_gcxc( &arho, &grho2[0], &sx, &sc, &v1x, &v2x, &v1c, &v2c );
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

    for(int ix=0;ix < this->pbasis;ix++) rhoout[ix] = rhoout[ix] - rho_core[ix];

    // 
    // ... second term of the gradient correction :
    // ... \sum_alpha (D / D r_alpha) ( D(rho*Exc)/D(grad_alpha rho) )
    // 
    double *dh = new double[this->pbasis]();

    CPP_app_grad_driver (this->L, this->T, vxc2, 
                         h, &h[this->pbasis], &h[2*this->pbasis],
                         this->dimx, this->dimy, this->dimz, 
                         this->hxgrid, this->hygrid, this->hzgrid, 
                         APP_CI_EIGHT);

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


    delete [] dh;
    delete [] h;
    delete [] vxc2;
    delete [] d2rho;
    delete [] grho;
    delete [] rhoout;


}


extern "C" const char *c_get_dft_name(void)
{
    return Functional::saved_dft_name.c_str();
}
