/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*
     	Just generates a random start.
*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "main.h"
#include "init_var.h"

static void init_wf_gamma(STATE * states);

void init_wf(STATE * states)
{

#if GAMMA_PT
    init_wf_gamma(states);
#else
    init_wf_complex(states);
#endif

}

static void init_wf_gamma(STATE * states)
{

    int idx, state, ix, iy, iz;
    STATE *sp;
    long idum;
    int ix1, iy1, iz1;
    int ixx, iyy, izz;
    rmg_double_t temp;
    int i;
    int idx1, idx2, idx3, idx4, idx5,idx6;


    if (pct.gridpe == 0)
        printf(" Begin init_wf ...\n");
    my_barrier();

    if (pct.gridpe == 0)
        printf(" Initialize random functions\n");
	
        dprintf(" Initialize random functions\n");

    /* Initialize the random number generator */

    /*  fill the borders of the orbit and smooth it */
    for (state = ct.state_begin; state < ct.state_end; state++)
    {
        idum = 3356 + state;  /*  seeds are different for different orbitals  */
        rand0(&idum);
        ixx = states[state].orbit_nx;
        iyy = states[state].orbit_ny;
        izz = states[state].orbit_nz;
        sp = &states1[state];
        for(ix = 0; ix < ixx; ix++)
        for(iy = 0; iy < iyy; iy++)
        for(iz = 0; iz < izz; iz++)
        {
            idx = ix * iyy * izz + iy * izz + iz;
            sp->psiR[idx] = 0.0;
        }

        for(ix = ixx/2 -3; ix < ixx/2+3; ix++)
        for(iy = iyy/2 -3; iy < iyy/2+3; iy++)
        for(iz = izz/2 -3; iz < izz/2+3; iz++)
        {
            idx = ix * iyy * izz + iy * izz + iz;
            sp->psiR[idx] = rand0(&idum);
        }

        for(ix = ixx/2 -4; ix < ixx/2+4; ix++)
        for(iy = iyy/2 -4; iy < iyy/2+4; iy++)
        for(iz = izz/2 -4; iz < izz/2+4; iz++)
        {
            idx = ix * iyy * izz + iy * izz + iz;
            idx1 = (ix-1) *iyy * izz + (iy+0) * izz + iz +0;
            idx2 = (ix+1) *iyy * izz + (iy+0) * izz + iz +0;
            idx3 = (ix+0) *iyy * izz + (iy-1) * izz + iz +0;
            idx4 = (ix+0) *iyy * izz + (iy+0) * izz + iz +0;
            idx5 = (ix+0) *iyy * izz + (iy+0) * izz + iz -1;
            idx6 = (ix+0) *iyy * izz + (iy+0) * izz + iz +1;

            states[state].psiR[idx] += (sp->psiR[idx1] +sp->psiR[idx2] +sp->psiR[idx3]
                   +sp->psiR[idx4] +sp->psiR[idx5] +sp->psiR[idx6])/6.0 ;
        }


    }



    normalize_orbits(states);
    /*
     *	ortho_norm_local(states); 
     */




    if (pct.gridpe == 0)
        printf(" init_wf done  \n");

#if     DEBUG
    print_state_sum(states);
    print_sum(pct.psi_size, states[ct.state_begin].psiR, "init_wf.c states sum");
#endif


}                               /* end init_wf */
