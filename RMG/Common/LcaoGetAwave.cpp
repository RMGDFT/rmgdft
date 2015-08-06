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

#include <complex>
#include "const.h"
#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "State.h"
#include "Kpoint.h"
#include "transition.h"

#if RMG_FAST_MATH
   #include "fastonebigheader.h"
#else
   #define fasterlog log 
   #define fasterexp exp
#endif

/*This function calculates atomic wavefunctions using wavefunctions read from PP files
 * with angular part added. The result is in psi, which is assumed to be initialized to zero*/

template void LcaoGetAwave(double *, ION *, int, int, int, double);
template void LcaoGetAwave(std::complex<double>  *, ION *, int, int, int, double);

template <typename StateType>
void LcaoGetAwave (StateType *psi, ION *iptr, int awave_idx, int l, int m, double coeff)
{


    int  yindex;
    double r, vector[3];

    int ix, iy, iz, ixx, iyy, izz;
    int xstart, ystart, zstart, xend, yend, zend;
    int i_r, idx;
    int ilow, jlow, klow, ihi, jhi, khi;
    int dimx, dimy, dimz;
    int PX0_GRID, PY0_GRID, PZ0_GRID;
    int PX_OFFSET, PY_OFFSET, PZ_OFFSET;
    int NX_GRID, NY_GRID, NZ_GRID;


    double r1, r2, fradius, coef1, coef2;
    double x[3];
    double hxgrid, hygrid, hzgrid;
    double xside, yside, zside;
    double a, b, c;
    SPECIES *sp;



    hxgrid = get_hxgrid();
    hygrid = get_hygrid();
    hzgrid = get_hzgrid();
    xside = get_xside();
    yside = get_yside();
    zside = get_zside();

    PX0_GRID = get_PX0_GRID();
    PY0_GRID = get_PY0_GRID();
    PZ0_GRID = get_PZ0_GRID();
    PX_OFFSET = get_PX_OFFSET();
    PY_OFFSET = get_PY_OFFSET();
    PZ_OFFSET = get_PZ_OFFSET();
    NX_GRID = get_NX_GRID();
    NY_GRID = get_NY_GRID();
    NZ_GRID = get_NZ_GRID();


    ilow = PX_OFFSET;
    jlow = PY_OFFSET;
    klow = PZ_OFFSET;
    ihi = ilow + PX0_GRID;
    jhi = jlow + PY0_GRID;
    khi = klow + PZ0_GRID;




    /*Starting index for ylm function: Indexing is as follows: 0:s, 1:px, 2:py, 3:pz, 4:dxx, etc.*/
    yindex = l*l + m;

    /* Get species type */
    sp = &ct.sp[iptr->species];

    // r[i] = a exp[(i*b] -c
    // for norm-conserving in FHI, c= 0
    // for ultrasoft with Vanderbilt code, c != 0.0

    b = log((sp->r[2] - sp->r[1])/(sp->r[1] - sp->r[0]));
    c = (sp->r[0] * exp(b) - sp->r[1])/(1.0 -exp(b) );
    a = sp->r[0] + c;



    dimx =  sp->aradius/(hxgrid*xside);
    dimy =  sp->aradius/(hygrid*yside);
    dimz =  sp->aradius/(hzgrid*zside);

    dimx = dimx * 2 + 1;
    dimy = dimy * 2 + 1;
    dimz = dimz * 2 + 1;

    xstart = iptr->xtal[0] / hxgrid - dimx/2;
    xend = xstart + dimx;
    ystart = iptr->xtal[1] / hygrid - dimy/2;
    yend = ystart + dimy;
    zstart = iptr->xtal[2] / hzgrid - dimz/2;
    zend = zstart + dimz;

    for (ix = xstart; ix < xend; ix++)
    {
        x[0] = ix * hxgrid - iptr->xtal[0];

        // fold the grid into the unit cell
        ixx = (ix + 20 * NX_GRID) % NX_GRID;
        if(ixx >= ilow && ixx < ihi)
        {

            for (iy = ystart; iy < yend; iy++)
            {
                x[1] = iy * hygrid - iptr->xtal[1];

                // fold the grid into the unit cell
                iyy = (iy + 20 * NY_GRID) % NY_GRID;
                if(iyy >= jlow && iyy < jhi)
                {
                    for (iz = zstart; iz < zend; iz++)
                    {
                        x[2] = iz * hzgrid - iptr->xtal[2];

                        // fold the grid into the unit cell
                        izz = (iz + 20 * NZ_GRID) % NZ_GRID;
                        if(izz >= klow && izz < khi)
                        {

                            idx = (ixx-ilow) * PY0_GRID * PZ0_GRID + (iyy-jlow) * PZ0_GRID + izz-klow;
                            r = Rmg_L.metric(x);

                            Rmg_L.to_cartesian(x, vector);
                        
                            if(r < sp->r[0])
                            {
                                fradius = sp->atomic_wave[awave_idx][0]; 
                            }
                            else
                            {
                                i_r = (int)(fasterlog ( (r+c)/a) /b);

                                r1 = sp->r[i_r];
                                r2 = sp->r[i_r+1];
                                coef1 = (r2-r)/(r2-r1);
                                coef2 = (r-r1)/(r2-r1);

                                fradius = coef1 * sp->atomic_wave[awave_idx][i_r] 
                                    + coef2 * sp->atomic_wave[awave_idx][i_r+1];
                            }

                            psi[idx] += coeff * fradius * ylm(yindex, vector);


                        }   
                    }
                }
            }
        }
    }
}

/******/


