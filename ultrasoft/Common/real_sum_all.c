/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/real_sum_all.c *****
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
 *   REAL real_sum_all(REAL x)
 *   Performs a scalar sum over all processors.
 * INPUTS
 *   x: defined in each processor
 * OUTPUT
 *   sum over all processors is returned to each processor
 * PARENTS
 *   get_ke.c get_rho.c get_te.c get_vh.c getpoi_bc.c init_nuc.c
 *   lforce.c mg_eig_state.c norm_psi.c 
 * CHILDREN
 *   nothing
 * SOURCE
 */


#include <float.h>
#include <math.h>
#include "main.h"



REAL real_sum_all (REAL x)
{

    REAL inreg;
    REAL outreg;
#if MD_TIMERS
    REAL time0;

    time0 = my_crtc ();
#endif
	


    inreg = x;

    MPI_Allreduce (&inreg, &outreg, 1, MPI_DOUBLE, MPI_SUM, pct.grid_comm);

#if MD_TIMERS
    rmg_timings (REAL_SUM_ALL_TIME, my_crtc () - time0, 0);
#endif


    return outreg;



}                               /* end real_sum_all */




REAL real_sum_all_spin (REAL x)
{

    REAL inreg;
    REAL outreg;
#if MD_TIMERS
    REAL time0;

    time0 = my_crtc ();
#endif
	


    inreg = x; 
    
    MPI_Allreduce (&inreg, &outreg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#if MD_TIMERS
    rmg_timings (REAL_SUM_ALL_TIME, my_crtc () - time0, 0);
#endif


    return outreg;

}                               /* end real_sum_all_spin */




















