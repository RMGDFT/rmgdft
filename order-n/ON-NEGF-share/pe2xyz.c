/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/pe2xyz.c *****
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
 *   void pe2xyz(int pe, int *x, int *y, int *z)
 *   Determines cartesian mapping for a given processor.
 *   Automatically wraps around if we are running with some processors
 *   duplicated.
 *   It is also usful for determine the corner coordinate for each processor
 * INPUTS
 *   pe: rank of the processor
 * OUTPUT
 *   x,y,z: coordinates of the processor in 3D array
 * PARENTS
 *   get_index.c getpoi_bc.c init_wf.c init_wflcao.c pack_vhdtos.c
 *   pack_vhstod.c init_pe.c set_bc.c symmetry.c write_avgd.c write_zstates.c
 * CHILDREN
 *   nothing  
 * SOURCE
 */


#include "md.h"



void pe2xyz(int pe, int *x, int *y, int *z)
{

    *x = pe;
    *z = *x % pct.pe_z;
    *x /= pct.pe_z;
    *y = *x % pct.pe_y;
    *x /= pct.pe_y;

    if (*x >= pct.pe_x)
        *x -= pct.pe_x;
    if (*x >= pct.pe_x)
        *x -= pct.pe_x;


}                               /* end pe2xyz */


/******/
