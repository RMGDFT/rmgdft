/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include "main.h"
#include "prototypes_on.h"
#include "init_var.h"
//#include "svnrev.h"




void dipole_calculation(double *rhooo, double *dipole)
{



    int i, j, k, idx;

    double x, y, z;

    dipole[0] = 0.0;
    dipole[1] = 0.0;
    dipole[2] = 0.0;

    for(i = 0; i < get_FPX0_GRID(); i++)
    {
        x = (get_FPX_OFFSET() + i)*get_hxxgrid() * get_xside();
        for(j = 0; j < get_FPY0_GRID(); j++)
        {
            y = (get_FPY_OFFSET() + j)*get_hyygrid() * get_yside();

            for(k = 0; k < get_FPZ0_GRID(); k++)
            {
                z = (get_FPZ_OFFSET() + k)*get_hzzgrid() * get_zside();

                idx = i * get_FPY0_GRID() * get_FPZ0_GRID() + j*get_FPZ0_GRID() + k;
                dipole[0] += x * rhooo[idx];
                dipole[1] += y * rhooo[idx];
                dipole[2] += z * rhooo[idx];
            }
        }
    }

    dipole[0] *= get_vel_f();
    dipole[1] *= get_vel_f();
    dipole[2] *= get_vel_f();
    
    idx = 3;
    global_sums (dipole, &idx, pct.grid_comm);


}
