/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/*

  This routine interpolates the orbitals to a new fine-grid used in NEGF 
  Calculation. For a 3/4-probe system this interpolation is carried out 
  along both x and y directions but for a 2-probe system it is required 
  only along the x-direction.

*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "main.h"

/* which_part = 0 for conductors 
 *              1 for left lead
 *              2 for right lead
 *              3 for down lead 
 *              4 for up lead 
 */

void interpolation_orbit (STATE * states)
{


    int ix, iy, iz, subsystem, idx0;
    int idx1, incx, incy;
    int NX, NY, NZ;
    int ixmin_old, ixmin_new, iymin_old, iymin_new;
    int st, st1, st_new, max_orbit_nx_ny;

    REAL *psi_old, *psi_new; 
    REAL hx_old, hx_new, hy_old, hy_new;
    REAL x1_old, x1_new, y1_old, y1_new;

/*  when the transport calculation use different gird space hx because of the different length
 *  in lead and conductor, interpolate the readin potential and rho to the current grid
 *  numbers of grids for both are same, the hx is only slightly different
 */

    hx_new = ct.hxgrid * ct.xside;
    hy_new = ct.hygrid * ct.yside;

    max_orbit_nx_ny = max(ct.max_orbit_nx, ct.max_orbit_ny);
    my_malloc_init( psi_old, max_orbit_nx_ny, REAL );
    my_malloc_init( psi_new, max_orbit_nx_ny, REAL );

#if 	LDEBUG
    printf ("\n PE: %d  xside %f  %f %f  ", pct.gridpe, lcr[1].xside, lcr[0].xside, lcr[2].xside);
    printf ("\n PE: %d  x_shift %f  %f %f  ", pct.gridpe, lcr[1].x_shift, lcr[0].x_shift,
            lcr[2].x_shift);
    printf ("\n PE: %d  NX_GRID %d  %d %d  ", pct.gridpe, lcr[1].NX_GRID, lcr[0].NX_GRID,
            lcr[2].NX_GRID);
    printf ("\n PE: %d  num_states %d  %d %d  ", pct.gridpe, lcr[1].num_states, lcr[0].num_states,
            lcr[2].num_states);
#endif

    for (st = ct.state_begin; st < ct.state_end; st++)
    {
        NX = states[st].orbit_nx;
        NY = states[st].orbit_ny;
        NZ = states[st].orbit_nz;

        incx = NY * NZ;
        incy = NZ * NX;

        ixmin_new = states[st].ixmin;
        ixmin_old = states[st].ixmin_old;
        iymin_new = states[st].iymin;
        iymin_old = states[st].iymin_old;

       /* printf ("\n %d %d %d st ixmn\n", st, ixmin, ixmin_old);*/
        x1_new = ixmin_new * hx_new;
        y1_new = iymin_new * hy_new;

        st_new = 0;
        for (subsystem = 0; subsystem < cei.num_subsystem; subsystem++)
        {
            idx0 = cei.subsystem_idx[subsystem];
       
            for (st1 = lcr[idx0].state_begin; st1 < lcr[idx0].state_end; st1++)
            {
                if(st_new == st)
                {
                    if(cei.num_probe < 3)
                    {
                        hx_old = lcr[idx0].xside / lcr[idx0].NX_GRID;
                        x1_old = lcr[idx0].x_shift + ixmin_old  * hx_old;

                        /* ==== Interpolation along x-direction ==== */
                        for (iy = 0; iy < NY; iy++)
                        {
                            for (iz = 0; iz < NZ; iz++)
                            {
                                for (ix = 0; ix < NX; ix++)
                                {
                                    idx1 = iz + iy * NZ + ix * incx;
                                    psi_old[ix] = states[st].psiR[idx1];
                                }
                                diff_hx_interpolation (st, psi_new, psi_old, NX, hx_new, hx_old, x1_new, x1_old);

                                for (ix = 0; ix < NX; ix++)
                                {
                                    idx1 = iz + iy * NZ + ix * incx;
                                    states[st].psiR[idx1] = psi_new[ix];
                                }
                            }
                        }

                    }
                    else if(cei.num_probe < 5)
                    {
                        hx_old = lcr[idx0].xside / lcr[idx0].NX_GRID;
                        x1_old = lcr[idx0].x_shift + ixmin_old * hx_old;

                        /* ==== Interpolation along x-direction ==== */
                        for (iy = 0; iy < NY; iy++)
                        {
                            for (iz = 0; iz < NZ; iz++)
                            {
                                for (ix = 0; ix < NX; ix++)
                                {
                                    idx1 = iz + iy * NZ + ix * incx;
                                    psi_old[ix] = states[st].psiR[idx1];
                                }
                                diff_hx_interpolation (st, psi_new, psi_old, NX, hx_new, hx_old, x1_new, x1_old);

                                for (ix = 0; ix < NX; ix++)
                                {
                                    idx1 = iz + iy * NZ + ix * incx;
                                    states[st].psiR[idx1] = psi_new[ix];
                                }
                            }
                        }


                        hy_old = lcr[idx0].yside / lcr[idx0].NY_GRID;
                        y1_old = lcr[idx0].y_shift + iymin_old  * hy_old;
/*
                        y1_old = lcr[idx0].y_shift + iymin_old * hy_old;
                        y1_new = (iymin_new + lcr[idx0].y0) * hy_new;
*/
/*
                        if(idx0 == 0 || idx0 == 3 || idx0 == 4)
                        printf (" old_urgent %d %d %f %f %d %f %d %f\n", st, idx0, 
                                       hy_old, hy_new, iymin_old, y1_old, iymin_new, y1_new);
*/

                        /* ==== Interpolation along y-direction ==== */

                        for (iz = 0; iz < NZ; iz++)
                        {
                            for (ix = 0; ix < NX; ix++)
                            {
                                for (iy = 0; iy < NY; iy++)
                                {
                                    idx1 = iz + iy * NZ + ix * incx;
                                    psi_old[iy] = states[st].psiR[idx1];
                                }
/*
                                if(pct.gridpe ==0) 
                                printf (" urgent %d %d %f %f %f %f \n", idx0, NY, hy_new, hy_old, y1_new, y1_old);
*/
                                diff_hx_interpolation (st, psi_new, psi_old, NY, hy_new, hy_old, y1_new, y1_old); 

                                for (iy = 0; iy < NY; iy++)
                                {
                                    idx1 = iz + iy * NZ + ix * incx;
                                    states[st].psiR[idx1] = psi_new[iy];
/*
                                    if(pct.gridpe ==0) 
                                    printf (" urgent %f %f \n", psi_old[iy], psi_new[iy]);
*/
                                }
                            }
                        }

                    } 
                    else
                    {
                        error_handler (" WARNING: Code may need to modify for num_probe > 4 ");
                    }


                } /* if statement ends */
            st_new++;
            }   /* st1 loop ends */ 
        }   /* subsystem loop ends */ 


#if 	LDEBUG
        printf ("\n %d %f %f state x x  %d ixmin_old %d  hx_old %f", st, x1, x1_old, ixmin,
                ixmin_old, hx_old);
        fflush (NULL);
#endif


    }                           /* end for st */

    my_free(psi_old);
    my_free(psi_new);

    my_barrier ();
    if (pct.gridpe == 0) printf ("\n interpolation_orbit  is done! ");

}
