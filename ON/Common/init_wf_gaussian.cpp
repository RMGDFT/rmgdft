/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*
     	Just generates a random start.
*/

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "main.h"
#include "prototypes_on.h"


#define s      0
#define px     1
#define py     2
#define pz     3
#define dxy    4
#define dxz    5
#define dyz    6
#define dx2y2  7
#define dz2    8
#define ss     9
#define s1     10
#define px1    11
#define py1    12
#define pz1    13

void init_wf_gaussian(STATE * states)
{

    /*
     *      number of orbitals 
     *      at the same center     gaussian orbitals
     *       1                     s
     *       2                     s,s*
     *       3                     s,s*, s1
     *       4                     s,p
     *       5                     s,p,s*
     *       6                     s,p,s*,s1
     *       7                     s,p,p1
     *       8                     s,p,s1,p1
     *       9                     s,p,d
     *      10                     s,p,d,s*
     */

    int idx, state, ix, iy, iz;
    int ixx, iyy, izz;
    int ist ;

    int n_orbital_same_center;
    int gaussian;
    double crds[3], hx, hy, hz, x, y, z, r2;
    double radius, alpha1, alpha2;


    if (pct.gridpe == 0)
        printf("\n initial orbitals with gaussian functions \n");
    MPI_Barrier(pct.img_comm);


    hx = get_hxgrid() * get_xside();
    hy = get_hygrid() * get_yside();
    hz = get_hzgrid() * get_zside();

    for (state = ct.state_begin; state < ct.state_end; state++)
    {

        n_orbital_same_center = states[state].n_orbital_same_center;
        ist = states[state].gaussian_orbital_index;
        crds[0] = states[state].crds[0];
        crds[1] = states[state].crds[1];
        crds[2] = states[state].crds[2];
        radius = states[state].radius;
        ixx = states[state].ixmax - states[state].ixmin + 1;
        iyy = states[state].iymax - states[state].iymin + 1;
        izz = states[state].izmax - states[state].izmin + 1;

        alpha1 =  (6.0 /radius) * (6.0 /radius);
        alpha2 =  (7.0 /radius) * (7.0 /radius);

        gaussian = ist;
        if (n_orbital_same_center <4) 
        {
            if(ist == 1 ) gaussian = ss;
            if(ist == 2 ) gaussian = s1;
        }

        if (n_orbital_same_center <7) 
        {
            if(ist == 4 ) gaussian = ss;
            if(ist == 5 ) gaussian = s1;
        }

        if (n_orbital_same_center ==7) 
        {
            if(ist == 4 ) gaussian = px1;
            if(ist == 5 ) gaussian = py1;
            if(ist == 6 ) gaussian = pz1;
        }

        if (n_orbital_same_center ==8) 
        {
            if(ist == 4 ) gaussian = s1;
            if(ist == 5 ) gaussian = px1;
            if(ist == 6 ) gaussian = py1;
            if(ist == 7 ) gaussian = pz1;
        }


        switch(gaussian)
        {

            case s: /*   first gaussian always be the s orbital */

                for(ix = 0; ix < ixx; ix++ )
                {
                    x = (ix + states[state].ixmin) *hx - crds[0];
                    for(iy = 0; iy < iyy; iy++ )
                    {
                        y = (iy + states[state].iymin) *hy - crds[1];
                        for(iz = 0; iz < izz; iz++ )
                        {
                            z = (iz + states[state].izmin) *hz - crds[2];
                            r2 = x*x + y*y + z*z;
                            idx = ix * iyy * izz + iy * izz + iz;

                            states[state].psiR[idx] = exp( - alpha1 * r2);
                        }
                    }
                }
               break;

            case px: 

                for(ix = 0; ix < ixx; ix++ )
                {
                    x = (ix + states[state].ixmin) *hx - crds[0];
                    for(iy = 0; iy < iyy; iy++ )
                    {
                        y = (iy + states[state].iymin) *hy - crds[1];
                        for(iz = 0; iz < izz; iz++ )
                        {
                            z = (iz + states[state].izmin) *hz - crds[2];
                            r2 = x*x + y*y + z*z;
                            idx = ix * iyy * izz + iy * izz + iz;

                            states[state].psiR[idx] = x * exp( - alpha1 * r2);
                        }
                    }
                }
               break;

            case py: 

                for(ix = 0; ix < ixx; ix++ )
                {
                    x = (ix + states[state].ixmin) *hx - crds[0];
                    for(iy = 0; iy < iyy; iy++ )
                    {
                        y = (iy + states[state].iymin) *hy - crds[1];
                        for(iz = 0; iz < izz; iz++ )
                        {
                            z = (iz + states[state].izmin) *hz - crds[2];
                            r2 = x*x + y*y + z*z;
                            idx = ix * iyy * izz + iy * izz + iz;

                            states[state].psiR[idx] = y * exp( - alpha1 * r2);
                        }
                    }
                }
               break;

            case pz: 

                for(ix = 0; ix < ixx; ix++ )
                {
                    x = (ix + states[state].ixmin) *hx - crds[0];
                    for(iy = 0; iy < iyy; iy++ )
                    {
                        y = (iy + states[state].iymin) *hy - crds[1];
                        for(iz = 0; iz < izz; iz++ )
                        {
                            z = (iz + states[state].izmin) *hz - crds[2];
                            r2 = x*x + y*y + z*z;
                            idx = ix * iyy * izz + iy * izz + iz;

                            states[state].psiR[idx] = z * exp( - alpha1 * r2);
                        }
                    }
                }
               break;


            case dxy: 

                for(ix = 0; ix < ixx; ix++ )
                {
                    x = (ix + states[state].ixmin) *hx - crds[0];
                    for(iy = 0; iy < iyy; iy++ )
                    {
                        y = (iy + states[state].iymin) *hy - crds[1];
                        for(iz = 0; iz < izz; iz++ )
                        {
                            z = (iz + states[state].izmin) *hz - crds[2];
                            r2 = x*x + y*y + z*z;
                            idx = ix * iyy * izz + iy * izz + iz;

                            states[state].psiR[idx] = x * y * exp( - alpha1 * r2);
                        }
                    }
                }
               break;

            case dxz: 

                for(ix = 0; ix < ixx; ix++ )
                {
                    x = (ix + states[state].ixmin) *hx - crds[0];
                    for(iy = 0; iy < iyy; iy++ )
                    {
                        y = (iy + states[state].iymin) *hy - crds[1];
                        for(iz = 0; iz < izz; iz++ )
                        {
                            z = (iz + states[state].izmin) *hz - crds[2];
                            r2 = x*x + y*y + z*z;
                            idx = ix * iyy * izz + iy * izz + iz;

                            states[state].psiR[idx] = x * z * exp( - alpha1 * r2);
                        }
                    }
                }
               break;

            case dyz: 

                for(ix = 0; ix < ixx; ix++ )
                {
                    x = (ix + states[state].ixmin) *hx - crds[0];
                    for(iy = 0; iy < iyy; iy++ )
                    {
                        y = (iy + states[state].iymin) *hy - crds[1];
                        for(iz = 0; iz < izz; iz++ )
                        {
                            z = (iz + states[state].izmin) *hz - crds[2];
                            r2 = x*x + y*y + z*z;
                            idx = ix * iyy * izz + iy * izz + iz;

                            states[state].psiR[idx] = y * z * exp( - alpha1 * r2);
                        }
                    }
                }
               break;

            case dx2y2: 

                for(ix = 0; ix < ixx; ix++ )
                {
                    x = (ix + states[state].ixmin) *hx - crds[0];
                    for(iy = 0; iy < iyy; iy++ )
                    {
                        y = (iy + states[state].iymin) *hy - crds[1];
                        for(iz = 0; iz < izz; iz++ )
                        {
                            z = (iz + states[state].izmin) *hz - crds[2];
                            r2 = x*x + y*y + z*z;
                            idx = ix * iyy * izz + iy * izz + iz;

                            states[state].psiR[idx] = (x*x - y*y) * exp( - alpha1 * r2);
                        }
                    }
                }
               break;

            case dz2: 
            printf("\n gaussian  aas %d  %d  %d\n", gaussian, ist, n_orbital_same_center);

                for(ix = 0; ix < ixx; ix++ )
                {
                    x = (ix + states[state].ixmin) *hx - crds[0];
                    for(iy = 0; iy < iyy; iy++ )
                    {
                        y = (iy + states[state].iymin) *hy - crds[1];
                        for(iz = 0; iz < izz; iz++ )
                        {
        //        printf("\n %d %d %d ix \n", ix,iy,iz);
                            z = (iz + states[state].izmin) *hz - crds[2];
                            r2 = x*x + y*y + z*z;
                            idx = ix * iyy * izz + iy * izz + iz;

                            states[state].psiR[idx] = (r2 - 3.0*z*z) * exp( - alpha1 * r2);
                        }
                    }
                }
               break;


            case ss: 

                for(ix = 0; ix < ixx; ix++ )
                {
                    x = (ix + states[state].ixmin) *hx - crds[0];
                    for(iy = 0; iy < iyy; iy++ )
                    {
                        y = (iy + states[state].iymin) *hy - crds[1];
                        for(iz = 0; iz < izz; iz++ )
                        {
                            z = (iz + states[state].izmin) *hz - crds[2];
                            r2 = x*x + y*y + z*z;
                            idx = ix * iyy * izz + iy * izz + iz;

                            states[state].psiR[idx] = r2 * exp( - alpha1 * r2);
                        }
                    }
                }
               break;

            case s1: 
                for(ix = 0; ix < ixx; ix++ )
                {
                    x = (ix + states[state].ixmin) *hx - crds[0];
                    for(iy = 0; iy < iyy; iy++ )
                    {
                        y = (iy + states[state].iymin) *hy - crds[1];
                        for(iz = 0; iz < izz; iz++ )
                        {
                            z = (iz + states[state].izmin) *hz - crds[2];
                            r2 = x*x + y*y + z*z;
                            idx = ix * iyy * izz + iy * izz + iz;

                            states[state].psiR[idx] = exp( - alpha2 * r2);
                        }
                    }
                }

               break;
            case px1: 

                for(ix = 0; ix < ixx; ix++ )
                {
                    x = (ix + states[state].ixmin) *hx - crds[0];
                    for(iy = 0; iy < iyy; iy++ )
                    {
                        y = (iy + states[state].iymin) *hy - crds[1];
                        for(iz = 0; iz < izz; iz++ )
                        {
                            z = (iz + states[state].izmin) *hz - crds[2];
                            r2 = x*x + y*y + z*z;
                            idx = ix * iyy * izz + iy * izz + iz;

                            states[state].psiR[idx] = x * exp( - alpha2 * r2);
                        }
                    }
                }
               break;

            case py1: 

                for(ix = 0; ix < ixx; ix++ )
                {
                    x = (ix + states[state].ixmin) *hx - crds[0];
                    for(iy = 0; iy < iyy; iy++ )
                    {
                        y = (iy + states[state].iymin) *hy - crds[1];
                        for(iz = 0; iz < izz; iz++ )
                        {
                            z = (iz + states[state].izmin) *hz - crds[2];
                            r2 = x*x + y*y + z*z;
                            idx = ix * iyy * izz + iy * izz + iz;

                            states[state].psiR[idx] = y * exp( - alpha2 * r2);
                        }
                    }
                }
               break;

            case pz1: 

                for(ix = 0; ix < ixx; ix++ )
                {
                    x = (ix + states[state].ixmin) *hx - crds[0];
                    for(iy = 0; iy < iyy; iy++ )
                    {
                        y = (iy + states[state].iymin) *hy - crds[1];
                        for(iz = 0; iz < izz; iz++ )
                        {
                            z = (iz + states[state].izmin) *hz - crds[2];
                            r2 = x*x + y*y + z*z;
                            idx = ix * iyy * izz + iy * izz + iz;

                            states[state].psiR[idx] = z * exp( - alpha2 * r2);
                        }
                    }
                }
               break;
            default:
                error_handler("Gaussian orbitals > 10 ");
        }

    }

    //normalize_orbits(states);
    ortho_norm_local(states);

    if (pct.gridpe == 0)
        printf(" initial orbitals  done  \n");


}                               /* end init_wf_atom */
