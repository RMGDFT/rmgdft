
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


#include "prototypes_on.h"
#include "init_var.h"

void VhcorrDipoleInit(double *vh_x, double *vh_y, double *vh_z, double *rhoc)
{
    int nfp0 = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);
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

    lpole = 1;
    fac2l1 = 1.0;
    for (i = 1; i <= 2*lpole+1; i+=2) fac2l1 *=i;

    alpha = xside /9.0;  // small enough so exp[-(r/alpha)^] decays to eps at boundary
    if(yside /9.0 < alpha) alpha = yside/9.0;
    if(zside /9.0 < alpha) alpha = zside/9.0;


    for(i = 0; i < FPX0_GRID; i++)
    {
        x = (FPX_OFFSET + i)*hxxgrid * xside;
        for(j = 0; j < FPY0_GRID; j++)
        {
            y = (FPY_OFFSET + j)*hyygrid * yside;

            for(k = 0; k < FPZ0_GRID; k++)
            {
                z = (FPZ_OFFSET + k)*hzzgrid * zside;

                idx = i * FPY0_GRID * FPZ0_GRID + j*FPZ0_GRID + k;

                dipole_ion[0] += x * rhoc[idx];
                dipole_ion[1] += y * rhoc[idx];
                dipole_ion[2] += z * rhoc[idx];
            }
        }
    }

    idx = 3;
    global_sums (dipole_ion, &idx, pct.grid_comm);

    // use the charge center as the dipole center  
    for (i = 0; i < 3; i++) 
    {
        dipole_ion[i] *= get_vel_f();
        dipole_center[i] = dipole_ion[i]/ct.nel;
    }

    //  make the center on the grid
    idx = (int)(dipole_center[0]/(hxxgrid * xside));
    dipole_center[0] = idx * hxxgrid * xside;
    idx = (int)(dipole_center[1]/(hyygrid * yside));
    dipole_center[1] = idx * hyygrid * yside;
    idx = (int)(dipole_center[2]/(hzzgrid * zside));
    dipole_center[2] = idx * hzzgrid * zside;



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
        printf("\n  n < 0 in gaussintegral \n");
        fflush(NULL);
        exit(0);
    }

    fac2n1 = 1.0;
    for(i = 1; i<=2*n-1; i+=2) fac2n1 *= i;

    vcorr = fac2n1 *rootpi /(pow(2, n+1)) * erf(r);

    sumj = 0.0;
    for(j = 0; j <n-1; j++)
    {
        fac2j1 = 1;
        for(i = 1; i<=2*j+1; i+=2) fac2j1 *= i;
        sumj += pow(r, 2*j+1)/pow(2, n-j) *fac2n1/fac2j1;
    }

    vcorr = vcorr - exp(-r*r) * sumj;

    return vcorr;
}







