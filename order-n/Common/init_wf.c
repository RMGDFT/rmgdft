/************************** SVN Revision Information **************************
 **    $Id: init_wf.c 587 2006-08-18 22:57:56Z miro $    **
******************************************************************************/
 
/*
     	Just generates a random start.
*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "md.h"

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
    REAL *xrand, *yrand, *zrand;
    long idum;
    int ix1, iy1, iz1;
    int ixx, iyy, izz;
    REAL temp;
    int i;


    my_malloc_init( xrand, NX_GRID + NY_GRID + NZ_GRID, REAL );
    yrand = xrand + NX_GRID;
    zrand = yrand + NY_GRID;

    if (pct.thispe == 0)
        printf(" Begin init_wf ...\n");
    my_barrier();

    if (pct.thispe == 0)
        printf(" Initialize random functions\n");

    /* Initialize the random number generator */
    idum = 3356;
    rand0(&idum);


    for (state = 0; state < ct.num_states; state++)
    {
        ixx = states[state].ixmax - states[state].ixmin + 1;
        iyy = states[state].iymax - states[state].iymin + 1;
        izz = states[state].izmax - states[state].izmin + 1;

        assert(ixx <= states[state].orbit_nx);
        assert(iyy <= states[state].orbit_ny);
        assert(izz <= states[state].orbit_nz);

        /* Generate x, y, z random number sequences */
        for (idx = 0; idx < NX_GRID; idx++)
            xrand[idx] = rand0(&idum) - 0.5;
        for (idx = 0; idx < NY_GRID; idx++)
            yrand[idx] = rand0(&idum) - 0.5;
        for (idx = 0; idx < NZ_GRID; idx++)
            zrand[idx] = rand0(&idum) - 0.5;

        sp = &states[state];
        if (state >= ct.state_begin && state < ct.state_end)
            for (ix = 0; ix < ixx; ix++)
                for (iy = 0; iy < iyy; iy++)
                    for (iz = 0; iz < izz; iz++)
                    {
                        ix1 = states[state].ixmin + ix;
                        if (ix1 < 0)
                            ix1 += NX_GRID;
                        if (ix1 >= NX_GRID)
                            ix1 -= NX_GRID;

                        iy1 = states[state].iymin + iy;
                        if (iy1 < 0)
                            iy1 += NY_GRID;
                        if (iy1 >= NY_GRID)
                            iy1 -= NY_GRID;

                        iz1 = states[state].izmin + iz;
                        if (iz1 < 0)
                            iz1 += NZ_GRID;
                        if (iz1 >= NZ_GRID)
                            iz1 -= NZ_GRID;

                        idx = ix * iyy * izz + iy * izz + iz;
                        sp->psiR[idx] = xrand[ix1] * yrand[iy1] * zrand[iz1];
                        sp->psiR[idx] = sp->psiR[idx] * sp->psiR[idx];
                    }
    }                           /* end for */


    my_barrier();

    /*  fill the borders of the orbit and smooth it */
    for (state = ct.state_begin; state < ct.state_end; state++)
    {
        ixx = states[state].ixmax - states[state].ixmin + 1;
        iyy = states[state].iymax - states[state].iymin + 1;
        izz = states[state].izmax - states[state].izmin + 1;
        sp = &states[state];
        app_mask(state, sp->psiR, 0);
        pack_ptos(sg_orbit, sp->psiR, ixx, iyy, izz);
        fill_orbit_borders(sg_orbit, ixx, iyy, izz);
        app_cir(sg_orbit, sp->psiR, ixx, iyy, izz);
        app_mask(state, sp->psiR, 0);

    }


    my_free(xrand);

    normalize_orbits(states);
/*
 *	ortho_norm_local(states); 
*/




    if (pct.thispe == 0)
        printf(" init_wf done  \n");

#if     DEBUG
    print_state_sum(states);
    print_sum(pct.psi_size, states[ct.state_begin].psiR, "init_wf.c states sum");
#endif


}                               /* end init_wf */
