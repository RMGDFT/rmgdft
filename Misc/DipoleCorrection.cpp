
/************************** SVN Revision Information **************************
 **    $Id: $    **
 ******************************************************************************/


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "RmgTimer.h"
#include "transition.h"
#include "blas.h"

#include "fft3d.h"
#include "RmgParallelFft.h"

void VhcorrPeriodicPart(double *vh_x, double *vh_y, double *vh_z, double alpha, double *r0);
void VhcorrDipoleInit(double *vh_x, double *vh_y, double *vh_z);
double gaussintegral(double r, int n);

void DipoleCorrection(double *dipole,  double *vh_dipole)
{
    static double *vh_x=NULL, *vh_y=NULL, *vh_z=NULL;
    int nfp0 = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);

    if(vh_x == NULL)
    {
        vh_x = new double[nfp0];
        vh_y = new double[nfp0];
        vh_z = new double[nfp0];
        VhcorrDipoleInit(vh_x, vh_y, vh_z);
    }

    if(ct.dipole_corr[0])
    {
        for(int idx = 0; idx < nfp0; idx++) 
        {
            vh_dipole[idx]  = vh_x[idx] *dipole[0];
        }
    }
    if(ct.dipole_corr[1])
    {
        for(int idx = 0; idx < nfp0; idx++) 
        {
            vh_dipole[idx] += vh_y[idx] *dipole[1];
        }
    }
    if(ct.dipole_corr[2])
    {
        for(int idx = 0; idx < nfp0; idx++) 
        {
            vh_dipole[idx] += vh_z[idx] *dipole[2];
        }
    }
}



void VhcorrDipoleInit(double *vh_x, double *vh_y, double *vh_z)
{
    int FPX0_GRID = Rmg_G->get_PX0_GRID(Rmg_G->default_FG_RATIO);
    int FPY0_GRID = Rmg_G->get_PY0_GRID(Rmg_G->default_FG_RATIO);
    int FPZ0_GRID = Rmg_G->get_PZ0_GRID(Rmg_G->default_FG_RATIO);
    int FPX_OFFSET = get_FPX_OFFSET();
    int FPY_OFFSET = get_FPY_OFFSET();
    int FPZ_OFFSET = get_FPZ_OFFSET();

    double xside = Rmg_L.get_xside();
    double yside = Rmg_L.get_yside();
    double zside = Rmg_L.get_zside();
    double hxxgrid = Rmg_G->get_hxgrid(Rmg_G->default_FG_RATIO);
    double hyygrid = Rmg_G->get_hygrid(Rmg_G->default_FG_RATIO);
    double hzzgrid = Rmg_G->get_hzgrid(Rmg_G->default_FG_RATIO);

    double r, x, y, z, alpha, fac2l1, vcorr;
    int idx, i,j,k, lpole;

    double dipole_center[3], dipole_ion[3];



    for (i = 0; i < 3; i++) 
    {
        dipole_ion[i] = 0.0;
    }
    for (int ion = 0; ion < ct.num_ions; ion++)
    {

        /*Ionic pointer and ionic charge */
        ION *iptr = &Atoms[ion];
        double icharge = Species[iptr->species].zvalence;

        /*Ionic contribution to dipole moment */
        dipole_ion[0] += icharge * iptr->crds[0];
        dipole_ion[1] += icharge * iptr->crds[1];
        dipole_ion[2] += icharge * iptr->crds[2];

    }                           /*end for (ion=0; ion<c->num_ions; ion++) */


    lpole = 1;
    fac2l1 = 1.0;
    for (i = 1; i <= 2*lpole+1; i+=2) fac2l1 *=i;

    alpha = xside /9.0;  // small enough so exp[-(r/alpha)^] decays to eps at boundary
    if(yside /9.0 < alpha) alpha = yside/9.0;
    if(zside /9.0 < alpha) alpha = zside/9.0;


    // use the charge center as the dipole center  
    for (i = 0; i < 3; i++) 
    {
        dipole_center[i] = dipole_ion[i]/ct.nel;
    }

    //  make the center on the grid
//    idx = (int)(dipole_center[0]/(hxxgrid * xside));
//    dipole_center[0] = idx * hxxgrid * xside;
//    idx = (int)(dipole_center[1]/(hyygrid * yside));
//    dipole_center[1] = idx * hyygrid * yside;
//    idx = (int)(dipole_center[2]/(hzzgrid * zside));
//    dipole_center[2] = idx * hzzgrid * zside;


    rmg_printf("\n dipole center at %f %f %f", dipole_center[0], dipole_center[1], dipole_center[2]);

    VhcorrPeriodicPart(vh_x, vh_y, vh_z, alpha, dipole_center);
    for(i = 0; i < FPX0_GRID; i++)
    {
        x = (FPX_OFFSET + i)*hxxgrid * xside - dipole_center[0];
        if(x > xside * 0.5) x = x- xside;
        if(x < -xside * 0.5) x = x + xside;


        for(j = 0; j < FPY0_GRID; j++)
        {
            y = (FPY_OFFSET + j)*hyygrid * yside - dipole_center[1];
            if(y > yside * 0.5) y = y- yside;
            if(y < -yside * 0.5) y = y + yside;

            for(k = 0; k < FPZ0_GRID; k++)
            {
                z = (FPZ_OFFSET + k)*hzzgrid * zside - dipole_center[2];
                if(z > zside * 0.5) z = z- zside;
                if(z < -zside * 0.5) z = z + zside;



                r = sqrt(x*x + y*y + z*z);
                if(r < 1.e-5)
                {
                    vcorr = 0.0;
                }
                else
                {
                    vcorr = gaussintegral(r/alpha, lpole);  //  for dipole only, l=1

                    vcorr *= sqrt(PI) *pow(2, lpole + 3)/fac2l1/pow(r, 2*lpole +1);

                    vcorr *= 3.0/(4.0 *PI);  // normalize Ylm = root(3/(4Pi))x for l =1 
                }

                idx = i * FPY0_GRID * FPZ0_GRID + j*FPZ0_GRID + k;
                vh_x[idx] = vcorr * x - vh_x[idx];
                vh_y[idx] = vcorr * y - vh_y[idx];
                vh_z[idx] = vcorr * z - vh_z[idx];


            }

        }

    }

}

double gaussintegral(double r, int n)
{

    // int(t^(2n)exp(-t*t)dt = -exp(-t*t) *{sum_j=0 to n-1 [(2n-1)!!/(2j+1)!!/2^(n-j) t^(2j+1)]}
    //                         +(2n-1)!! root(pi) /2^(n+1) erf(t)   
    //  when t = 0.0, function = 0.0

    int i, j;
    double fac2n1, fac2j1;
    double rootpi = sqrt(PI);

    double vcorr, sumj;

    if(n < 0)
    {
        rmg_printf("\n  n < 0 in gaussintegral \n");
        fflush(NULL);
        exit(0);
    }

    fac2n1 = 1.0;
    for(i = 1; i<=2*n-1; i+=2) fac2n1 *= i;

    vcorr = fac2n1 *rootpi /(pow(2, n+1)) * erf(r);

    sumj = 0.0;
    for(j = 0; j <=n-1; j++)
    {
        fac2j1 = 1;
        for(i = 1; i<=2*j+1; i+=2) fac2j1 *= i;
        sumj += pow(r, 2*j+1)/pow(2, n-j) *fac2n1/fac2j1;
    }

    vcorr = vcorr - exp(-r*r) * sumj;

    return vcorr;
}


#include "fft3d.h"
#include "RmgParallelFft.h"

void VhcorrPeriodicPart(double *vh_x, double *vh_y, double *vh_z, double alpha, double *r0)
{

    double gsquare, gx, g_r0;
    std::complex<double> phase_r0;

    //double g2cut = (sqrt(fine_pwaves->gmax))*(sqrt(fine_pwaves->gmax));
    //int global_basis = fine_pwaves->global_basis;
    int pbasis = fine_pwaves->pbasis;

    std::complex<double> *crho = new std::complex<double>[pbasis];

    //  for(int ig=0;ig < pbasis;ig++) {
    //      if(pwaves.gmags[ig] > g2cut) {
    //          crho[ig] = std::complex<double>(0.0, 0.0);
    //      }
    //  }


    double tpiba = 2.0 * PI / Rmg_L.celldm[0];
    double tpiba2 = tpiba * tpiba;

    for(int ig=0;ig < pbasis;ig++) {
        if(fine_pwaves->gmags[ig] > 1.0e-6) 
        {   
            gsquare = fine_pwaves->gmags[ig] * tpiba2;
            g_r0  = fine_pwaves->g[ig].a[0] *r0[0];
            g_r0 += fine_pwaves->g[ig].a[1] *r0[1];
            g_r0 += fine_pwaves->g[ig].a[2] *r0[2];
            g_r0 *= tpiba;
            phase_r0 = exp(std::complex<double>(0.0, -g_r0));

            gx  = fine_pwaves->g[ig].a[0] * tpiba;
            crho[ig] = exp(-alpha * alpha * gsquare/4.0)/gsquare * gx *phase_r0;
        }
        else
        {
            crho[ig] = 0.0;
        }
        
    }

    fine_pwaves->FftInverse(crho, crho);

    for(int i = 0;i < pbasis;i++) vh_x[i] = std::imag(crho[i]) * 4 * PI/Rmg_L.get_omega();


    for(int ig=0;ig < pbasis;ig++) {
        if(fine_pwaves->gmags[ig] > 1.0e-6) 
        {   
            gsquare = fine_pwaves->gmags[ig] * tpiba2;
            g_r0  = fine_pwaves->g[ig].a[0] *r0[0];
            g_r0 += fine_pwaves->g[ig].a[1] *r0[1];
            g_r0 += fine_pwaves->g[ig].a[2] *r0[2];
            g_r0 *= tpiba;
            phase_r0 = exp(std::complex<double>(0.0, -g_r0));

            gx  = fine_pwaves->g[ig].a[1] * tpiba;
            crho[ig] = exp(-alpha * alpha * gsquare/4.0)/gsquare * gx *phase_r0;
        }
        else
        {
            crho[ig] = 0.0;
        }
        
    }

    fine_pwaves->FftInverse(crho, crho);
    for(int i = 0;i < pbasis;i++) vh_y[i] = std::imag(crho[i]) * 4 * PI/Rmg_L.get_omega();

    for(int ig=0;ig < pbasis;ig++) {
        if(fine_pwaves->gmags[ig] > 1.0e-6) 
        {   
            gsquare = fine_pwaves->gmags[ig] * tpiba2;
            g_r0  = fine_pwaves->g[ig].a[0] *r0[0];
            g_r0 += fine_pwaves->g[ig].a[1] *r0[1];
            g_r0 += fine_pwaves->g[ig].a[2] *r0[2];
            g_r0 *= tpiba;
            phase_r0 = exp(std::complex<double>(0.0, -g_r0));

            //rmg_printf("\n  iggg %d %e %e %e %e %e", ig, fine_pwaves->g[ig].a[2],g_r0, g_r0/PI, phase_r0)

            gx  = fine_pwaves->g[ig].a[2] * tpiba;
            crho[ig] = exp(-alpha * alpha * gsquare/4.0)/gsquare * gx *phase_r0;
        }
        else
        {
            crho[ig] = 0.0;
        }
        
    }

    fine_pwaves->FftInverse(crho, crho);
    for(int i = 0;i < pbasis;i++) vh_z[i] = std::imag(crho[i]) * 4 * PI/Rmg_L.get_omega();


    delete [] crho;

}



