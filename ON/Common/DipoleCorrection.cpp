
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

void DipoleCorrection(double *rho_tot, double *rhoc, double *vh_dipole, double *vh_x, double *vh_y, double *vh_z)
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

    double x, y, z;
    int idx, i,j,k;

    double dipole_ele[3], dipole_ion[3];

    for (i = 0; i < 3; i++) 
    {
        dipole_ele[i] = 0.0;
        dipole_ion[i] = 0.0;
    }

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

                dipole_ele[0] += x * rho_tot[idx];
                dipole_ele[1] += y * rho_tot[idx];
                dipole_ele[2] += z * rho_tot[idx];
                dipole_ion[0] += x * rhoc[idx];
                dipole_ion[1] += y * rhoc[idx];
                dipole_ion[2] += z * rhoc[idx];
            }
        }
    }

    idx = 3;
    global_sums (dipole_ele, &idx, pct.grid_comm);
    global_sums (dipole_ion, &idx, pct.grid_comm);

    for (i = 0; i < 3; i++) 
    {
        dipole_ele[i] *= get_vel_f();
        dipole_ion[i] *= get_vel_f();
    }


    for(idx = 0; idx < nfp0; idx++) 
    {
        vh_dipole[idx]  = vh_x[idx] *(dipole_ele[0] - dipole_ion[0]);
        vh_dipole[idx] += vh_y[idx] *(dipole_ele[1] - dipole_ion[1]);
        vh_dipole[idx] += vh_z[idx] *(dipole_ele[2] - dipole_ion[2]);
    }
}

