/************************** SVN Revision Information **************************
 **    $Id: pack_stop.c 2093 2014-02-06 01:02:40Z ebriggs $    **
******************************************************************************/

/****f* QMD-MGDFT/pack_stop.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5
 * COPYRIGHT
 *   Copyright (C) 1995  Emil Briggs
 *   Copyright (C) 1998  Emil Briggs, Charles Brabec, Mark Wensell, 
 *                       Dan Sullivan, Chris Rapcewicz, Jerzy Bernholc
 *   Copyright (C) 2001  Emil Briggs, Wenchang Lu,
 *                       Marco Buongiorno Nardelli,Charles Brabec, 
 *                       Mark Wensell,Dan Sullivan, Chris Rapcewicz,
 *                       Jerzy Bernholc
 * FUNCTION
 *   void pack_stop(double *sg, double *pg, int dimx, int dimy, int dimz)
 *   pack smoothing grids sg to processor grids pg
 * INPUTS
 *   sg[(dimx+2)*(dimy+2)*(dimz+2)]: array with images
 *   dimx, dimy, dimz: dimensions of array 
 * OUTPUT
 *   pg[dimx*dimy*dimz]: processor grid, its value are copied from sg
 * PARENTS
 *   get_vh.c
 * CHILDREN
 *   nothing
 * SOURCE
 */


#include "packfuncs.h"



extern "C" void pack_ptos(double * sg, double * pg, int dimx, int dimy, int dimz)
{
    CPP_pack_ptos<double> (sg, pg, dimx, dimy, dimz);
}

extern "C" void pack_ptos_f(float * sg, float * pg, int dimx, int dimy, int dimz)
{
    CPP_pack_ptos<float> (sg, pg, dimx, dimy, dimz);
}

extern "C" void pack_stop (double * sg, double * pg, int dimx, int dimy, int dimz)
{
    CPP_pack_stop<double> (sg, pg, dimx, dimy, dimz);
}

extern "C" void pack_stop_f (float * sg, float * pg, int dimx, int dimy, int dimz)
{
    CPP_pack_stop<float> (sg, pg, dimx, dimy, dimz);
}
