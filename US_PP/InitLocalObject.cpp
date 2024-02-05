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

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "const.h"
#include "RmgTimer.h"
#include "RmgException.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "InputKey.h"
#include "blas.h"
//#include "main.h"
#include "init_var.h"
#include "transition.h"
#include "prototypes_on.h"
#include "Kbpsi.h"
#include "rmgthreads.h"
#include "RmgGemm.h"
#include "AtomicInterpolate.h"
#include "Atomic.h"
#include "RmgException.h"
#include "transition.h"

// This is used to initialize 4 types of atomic data structures that live
// on the high density grid. Each object can have a sum representation
// where the sum is over ions and a local representation which separates
// the contributions of individual ions and is used to calculate forces. 
// The object types are.
//
// ATOMIC_LOCAL_PP  Local ionic pseudopotential
// ATOMIC_RHO       Atomic rho
// ATOMIC_RHOCOMP   Atomic compensating charges
// ATOMIC_RHOCORE   Atomic core charges for non-linear core corrections.
//

void LcaoGetAtomicRho(double *arho)
{
    InitLocalObject(arho, ATOMIC_RHO);
}

void InitLocalObject(double *sumobject, int object_type)
{
    double *dum_array = NULL;
    InitLocalObject(sumobject, dum_array, object_type, false);
}

void InitLocalObject (double *sumobject, double * &lobject, int object_type, bool compute_lobject)
{

    int ilow, jlow, klow, ihi, jhi, khi;
    int FP0_BASIS;
    int FPX0_GRID, FPY0_GRID, FPZ0_GRID;
    int FPX_OFFSET, FPY_OFFSET, FPZ_OFFSET;
    int FNX_GRID, FNY_GRID, FNZ_GRID;

    double hxxgrid, hyygrid, hzzgrid;
    double xside, yside, zside;

    hxxgrid = get_hxxgrid();
    hyygrid = get_hyygrid();
    hzzgrid = get_hzzgrid();
    xside = get_xside();
    yside = get_yside();
    zside = get_zside();

    FP0_BASIS = get_FP0_BASIS();
    FPX0_GRID = get_FPX0_GRID();
    FPY0_GRID = get_FPY0_GRID();
    FPZ0_GRID = get_FPZ0_GRID();
    FPX_OFFSET = get_FPX_OFFSET();
    FPY_OFFSET = get_FPY_OFFSET();
    FPZ_OFFSET = get_FPZ_OFFSET();
    FNX_GRID = get_FNX_GRID();
    FNY_GRID = get_FNY_GRID();
    FNZ_GRID = get_FNZ_GRID();


    ilow = FPX_OFFSET;
    jlow = FPY_OFFSET;
    klow = FPZ_OFFSET;
    ihi = ilow + FPX0_GRID;
    jhi = jlow + FPY0_GRID;
    khi = klow + FPZ0_GRID;

    /* Loop over ions determine num of ions whose local pp has overlap with this pe */

    pct.num_loc_ions = 0;
    for (int ion = 0; ion < ct.num_ions; ion++)
    {
        /* Generate ion pointer */
        ION *iptr = &Atoms[ion];

        /* Get species type */
        SPECIES *sp = &Species[iptr->species];

        int dimx =  sp->lradius/(hxxgrid*xside);
        int dimy =  sp->lradius/(hyygrid*yside);
        int dimz =  sp->lradius/(hzzgrid*zside);

        dimx = dimx * 2 + 1;
        dimy = dimy * 2 + 1;
        dimz = dimz * 2 + 1;
        

        int xstart = iptr->xtal[0] / hxxgrid - dimx/2;
        int xend = xstart + dimx;
        int ystart = iptr->xtal[1] / hyygrid - dimy/2;
        int yend = ystart + dimy;
        int zstart = iptr->xtal[2] / hzzgrid - dimz/2;
        int zend = zstart + dimz;


        bool map_x = false, map_y = false, map_z = false;
        for (int ix = xstart; ix < xend; ix++)
        {
            // fold the grid into the unit cell
            // maxium fold 20 times, local potential can extend to 20-1 unit cell
            int ixx = (ix + 20 * FNX_GRID) % FNX_GRID;
            if(ixx >= ilow && ixx < ihi)
            {
                map_x = true;
                break;
            }
        }

        if(!map_x) continue;

        for (int iy = ystart; iy < yend; iy++)
        {
            // fold the grid into the unit cell
            int iyy = (iy + 20 * FNY_GRID) % FNY_GRID;
            if(iyy >= jlow && iyy < jhi)
            {
                map_y = true;
                break;
            }
        }

        if(!map_y) continue;


        for (int iz = zstart; iz < zend; iz++)
        {
            // fold the grid into the unit cell
            int izz = (iz + 20 * FNZ_GRID) % FNZ_GRID;
            if(izz >= klow && izz < khi)
            {
                map_z = true;
                break;
            }
        }

        if(!map_z) continue;

        pct.loc_ions_list[pct.num_loc_ions] = ion;
        pct.num_loc_ions ++;
    }

    // Zero out sumobject

    int factor = 1;
    if(object_type == ATOMIC_RHOCORE_STRESS)
    {
        factor = 3;
    }

    if( (ct.nspin == 4) && (object_type == ATOMIC_RHO) && !compute_lobject)
    {
        factor = 4;
    }

    for(int idx = 0; idx < factor* FP0_BASIS; idx++) sumobject[idx] = 0.0;
    if(object_type == ATOMIC_RHOCOMP) {
        int npes = get_PE_X() * get_PE_Y() * get_PE_Z();
        double t1 = ct.background_charge / (double)FP0_BASIS / get_vel_f() / (double)npes;
        for (int idx = 0; idx < FP0_BASIS; idx++) sumobject[idx] = t1;

    }

    // Responsibility for freeing lobject lies in the calling routine
    size_t alloc = (size_t)pct.num_loc_ions * (size_t)FP0_BASIS + 128;
    if(compute_lobject) lobject = new double[alloc]();

#pragma omp parallel 
    {
        double *sumobj_omp = new double[factor * FP0_BASIS](); 
        double *lobject_omp = new double[alloc]();
#pragma omp barrier
#pragma omp for schedule(static, 1) nowait
        for (int ion1 = 0; ion1 < pct.num_loc_ions; ion1++)
        {
            int ion = pct.loc_ions_list[ion1];
            /* Generate ion pointer */
            ION *iptr = &Atoms[ion];

            /* Get species type */
            SPECIES *sp = &Species[iptr->species];
            double Zv = sp->zvalence;
            double rc = sp->rc;
            double rc2 = rc * rc;
            double rcnorm = rc * rc * rc * pow (PI, 1.5);
            rcnorm = 1.0 / rcnorm;

            int dimx =  sp->lradius/(hxxgrid*xside);
            int dimy =  sp->lradius/(hyygrid*yside);
            int dimz =  sp->lradius/(hzzgrid*zside);

            dimx = dimx * 2 + 1;
            dimy = dimy * 2 + 1;
            dimz = dimz * 2 + 1;


            int xstart = iptr->xtal[0] / hxxgrid - dimx/2;
            int xend = xstart + dimx;
            int ystart = iptr->xtal[1] / hyygrid - dimy/2;
            int yend = ystart + dimy;
            int zstart = iptr->xtal[2] / hzzgrid - dimz/2;
            int zend = zstart + dimz;

            for (int ix = xstart; ix < xend; ix++)
            {
                // fold the grid into the unit cell
                // maxium fold 20 times, local potential can extend to 20-1 unit cell
                int ixx = (ix + 20 * FNX_GRID) % FNX_GRID;
                if(ixx >= ilow && ixx < ihi)
                {

                    for (int iy = ystart; iy < yend; iy++)
                    {
                        // fold the grid into the unit cell
                        int iyy = (iy + 20 * FNY_GRID) % FNY_GRID;
                        if(iyy >= jlow && iyy < jhi)
                        {
                            for (int iz = zstart; iz < zend; iz++)
                            {
                                // fold the grid into the unit cell
                                int izz = (iz + 20 * FNZ_GRID) % FNZ_GRID;
                                if(izz >= klow && izz < khi)
                                {

                                    int idx = (ixx-ilow) * FPY0_GRID * FPZ0_GRID + (iyy-jlow) * FPZ0_GRID + izz-klow;
                                    double x[3], cx[3];
                                    x[0] = ix * hxxgrid - iptr->xtal[0];
                                    x[1] = iy * hyygrid - iptr->xtal[1];
                                    x[2] = iz * hzzgrid - iptr->xtal[2];
                                    double r = Rmg_L.metric (x);
                                    Rmg_L.to_cartesian(x, cx);

                                    if(r > sp->lradius) continue;
                                    double t1;

                                    switch(object_type) 
                                    {
                                        case ATOMIC_LOCAL_PP:
                                            t1= AtomicInterpolateInline (&sp->localig[0], r);
                                            break;

                                        case ATOMIC_RHO:
                                            t1= AtomicInterpolateInline (&sp->arho_lig[0], r);
                                            break;

                                        case ATOMIC_RHOCOMP:
                                            t1= Zv * exp (-r * r / rc2) * rcnorm;
                                            break;

                                        case ATOMIC_RHOCORE: 
                                            if(sp->nlccflag) 
                                            {
                                                t1 = AtomicInterpolateInline (&sp->rhocorelig[0], r);
                                            }
                                            else
                                            {
                                                t1 = 0.0;
                                            }
                                            break;
                                        case ATOMIC_RHOCORE_STRESS: 
                                            if(sp->nlccflag) 
                                            {
                                                t1 = AtomicInterpolateInline (&sp->rhocorelig[0], r);
                                            }
                                            else
                                            {
                                                t1 = 0.0;
                                            }
                                            break;

                                        default:
                                            throw RmgFatalException() << "Undefined local object type" << 
                                                " in " << __FILE__ << " at line " << __LINE__ << "\n";

                                    }

                                    if( (ct.nspin == 2) && (object_type == ATOMIC_RHO) && !compute_lobject)
                                    { 
                                        if (pct.spinpe == 0)
                                            sumobj_omp[idx] += t1 * (0.5 + iptr->init_spin_rho) ;
                                        else
                                            sumobj_omp[idx] += t1 * (0.5 - iptr->init_spin_rho) ;

                                    }
                                    else if( (ct.nspin == 4) && (object_type == ATOMIC_RHO) && !compute_lobject)
                                    { 
                                        sumobj_omp[idx] += t1;
                                        sumobj_omp[idx+FP0_BASIS] += t1 * iptr->init_spin_x   ;
                                        sumobj_omp[idx+2*FP0_BASIS] += t1 * iptr->init_spin_y ;
                                        sumobj_omp[idx+3*FP0_BASIS] += t1 * iptr->init_spin_z ;
                                    }
                                    else if(object_type == ATOMIC_RHOCORE_STRESS)
                                    {
                                        sumobj_omp[0*FP0_BASIS + idx] += t1* cx[0];
                                        sumobj_omp[1*FP0_BASIS + idx] += t1* cx[1];
                                        sumobj_omp[2*FP0_BASIS + idx] += t1* cx[2];
                                    }
                                    else
                                    {
                                        sumobj_omp[idx] += t1;
                                    }
                                    if(compute_lobject) lobject_omp[(size_t)ion1 * (size_t)FP0_BASIS + idx] += t1;

                                }                           /* end for */

                            }
                        }
                    }
                }

            }
        }

#pragma omp critical(init_lo1)
        for(int idx = 0; idx < factor* FP0_BASIS; idx++) sumobject[idx] += sumobj_omp[idx];

#pragma omp critical(init_lo2)
        if(compute_lobject) for(int idx = 0; idx < factor* FP0_BASIS; idx++) lobject[idx] += lobject_omp[idx];
        delete [] sumobj_omp;
        delete [] lobject_omp;

    }

    // Renormalize atomic rho
    if(object_type == ATOMIC_RHO) {
        double t2 = 0.0;
        for (int idx = 0; idx < FP0_BASIS; idx++) t2 += sumobject[idx];
        t2 = get_vel_f() *  real_sum_all (t2, pct.grid_comm);
        t2 = real_sum_all (t2, pct.spin_comm);
        if(ct.AFM) t2 *= 2.0;

        double t1 = ct.nel / t2;
        double difference = fabs(t1 - 1.0);
        if ((ct.verbose == 1) || (difference > 0.05))
        {
            if (pct.imgpe == 0)
                printf ("\n LCAO initialization: Normalization constant for initial atomic charge is %f\n", t1);
        }

        for(int idx = 0;idx < FP0_BASIS;idx++) sumobject[idx] *= t1;
        if(!compute_lobject)
            for(int is = 1; is < ct.noncoll_factor * ct.noncoll_factor; is++)
                for(int idx = 0;idx < FP0_BASIS;idx++) sumobject[idx + is * FP0_BASIS] *= t1;

    }

    // Core charges may have a small negative component because of filtering
    if(object_type == ATOMIC_RHOCORE) 
    {
        for (int idx = 0; idx < FP0_BASIS; idx++) 
        {
            if(sumobject[idx] < 0.0) sumobject[idx] = 0.0;
        }
    }

    if(object_type == ATOMIC_RHOCOMP)
    {

        /* Check compensating charges */
        ct.crho = 0.0;
        for (int idx = 0; idx < FP0_BASIS; idx++) ct.crho += sumobject[idx];

        ct.crho = ct.crho * get_vel_f();
        ct.crho = real_sum_all (ct.crho, pct.grid_comm);  /* sum over pct.grid_comm  */

        if (ct.verbose)
        {
            if (pct.imgpe==0)
                printf("\nCompensating charge is %.8e\n", ct.crho);
        }
    }

    if(ct.runflag == RESTART || ct.runflag == Restart_TDDFT || ct.forceflag == TDDFT) return;

    if(object_type == ATOMIC_LOCAL_PP) init_efield (sumobject);

}   // end InitLocalObject
