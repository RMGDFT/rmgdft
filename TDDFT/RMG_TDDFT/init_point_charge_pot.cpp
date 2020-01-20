/****f* QMD-MGDFT/init_efield.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5
 * COPYRIGHT
 *   Copyright (C) 2001  Vincent Meunier, Wenchang Lu, Jerzy Bernholc
 * FUNCTION
 *   void init_efield(double *vnuc)
 *   Set up the electric field (add it to vnuc already defined)
 *   sawtooth potential
 * INPUTS
 *   vnuc: local part of pseudopotential
 * OUTPUT
 *   vnuc: local part of pseudopotential + V_efield 
 * PARENTS
 *   init_nuc.c
 * CHILDREN
 *   nothing
 * SOURCE
 */



#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>     
#include <algorithm>    
#include <cfloat>
#include <math.h>       
#include <mpi.h>       
#include "RmgException.h"
#include "transition.h"
#include "GlobalSums.h"
#include "Atomic.h"
#include <boost/math/special_functions/erf.hpp>

void init_point_charge_pot (double * vpot, int density)
{

    if (pct.gridpe == 0)
    {
        printf ("\n TDDFT point charge at: (%12.8f, %12.8f, %12.8f)", 
                ct.tddft_qpos[0],ct.tddft_qpos[1],ct.tddft_qpos[2]);
        printf ("\n ====================================== \n");
    }


    int PX0_GRID = Rmg_G->get_PX0_GRID(density);
    int PY0_GRID = Rmg_G->get_PY0_GRID(density);
    int PZ0_GRID = Rmg_G->get_PZ0_GRID(density);
    int PX_OFFSET = Rmg_G->get_PX_OFFSET(density);
    int PY_OFFSET = Rmg_G->get_PY_OFFSET(density);
    int PZ_OFFSET = Rmg_G->get_PZ_OFFSET(density);
    double hx = Rmg_G->get_hxgrid(density);
    double hy = Rmg_G->get_hygrid(density);
    double hz = Rmg_G->get_hzgrid(density);
    double xc = PX_OFFSET * hx; 
    double yc = PY_OFFSET * hy;
    double zc = PZ_OFFSET * hz;

    double xtal[3], dist, xcry[3];
    double qxtal[3];
    Rmg_L.to_crystal(qxtal, ct.tddft_qpos);
    for(int ix = 0; ix < PX0_GRID; ix++)
        for(int iy = 0; iy < PY0_GRID; iy++)
            for(int iz = 0; iz < PZ0_GRID; iz++)
            {
                int idx = ix * PY0_GRID * PZ0_GRID + iy * PZ0_GRID + iz;
                xtal[0] = xc + ix * hx;
                xtal[1] = yc + iy * hy;
                xtal[2] = zc + iz * hz;

                xcry[0] = xtal[0] - qxtal[0];
                xcry[1] = xtal[1] - qxtal[1];
                xcry[2] = xtal[2] - qxtal[2];

                if(xcry[0] > 0.5) xcry[0] -= 1.;
                if(xcry[1] > 0.5) xcry[1] -= 1.0;
                if(xcry[2] > 0.5) xcry[2] -= 1.0;
                if(xcry[0] < -0.5) xcry[0] += 1.0;
                if(xcry[1] < -0.5) xcry[1] += 1.0;
                if(xcry[2] < -0.5) xcry[2] += 1.0;

                dist = Rmg_L.metric(xcry);

                if( dist < 1.0e-5)
                {
                    vpot[idx] = 2.0/std::sqrt(PI)/ct.tddft_qgau;
                }
                else
                {
                    vpot[idx] = boost::math::erf (dist/ct.tddft_qgau) / dist;
                }

            }

    for(int idx = 0; idx < PZ0_GRID; idx++)
        printf("\n %d %f  ppp", idx, vpot[idx]);
}

