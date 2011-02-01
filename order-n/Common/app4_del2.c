/************************** SVN Revision Information **************************
 **    $Id: app4_del2.c 587 2006-08-18 22:57:56Z miro $    **
******************************************************************************/
 
/*

                            app4_del2.c


    Functions to apply the sixth order del-squared operator on a
    cubic grid.


*/



#include "md.h"
#include <float.h>
#include <math.h>



void app4_del2(REAL * f, REAL * work, int dimx, int dimy, int dimz)
{

    int iz, ix, iy;
    REAL h2, t0, t1x, t2x;
    REAL t1y, t2y;
    REAL t1z, t2z;
    REAL *dum2;
    int incx = dimy * dimz;
    int incx2 = (dimy + 4) * (dimz + 4);
    int incy = dimz;
    int incy2 = dimz + 4;
    int ixs, ixs1, ixs2, ixs3, ixs4, ixs5;
    int iys, iys1, iys2, iys3, iys4, iys5;


#if MD_TIMERS
    REAL time1, time2;
    time1 = my_crtc();
#endif

    my_malloc_init( dum2, (dimx + 4) * (dimy + 4) * (dimz + 4), REAL );

    fill_orbit_borders2(dum2, f, dimx, dimy, dimz);

    h2 = ct.hxgrid * ct.hxgrid * ct.xside * ct.xside;
    t0 = -30.0 / (12.0 * h2);
    t1x = 16.0 / (12.0 * h2);
    t2x = -1.0 / (12.0 * h2);

    h2 = ct.hygrid * ct.hygrid * ct.yside * ct.yside;
    t0 -= 30.0 / (12.0 * h2);
    t1y = 16.0 / (12.0 * h2);
    t2y = -1.0 / (12.0 * h2);

    h2 = ct.hzgrid * ct.hzgrid * ct.zside * ct.zside;
    t0 -= 30.0 / (12.0 * h2);
    t1z = 16.0 / (12.0 * h2);
    t2z = -1.0 / (12.0 * h2);



    for (ix = 2; ix < dimx + 2; ix++)
    {
        ixs = (ix - 2) * incx;
        ixs1 = ix * incx2;
        ixs2 = (ix - 1) * incx2;
        ixs3 = (ix + 1) * incx2;
        ixs4 = (ix - 2) * incx2;
        ixs5 = (ix + 2) * incx2;

        for (iy = 2; iy < dimy + 2; iy++)
        {
            iys = (iy - 2) * incy;
            iys1 = iy * incy2;
            iys2 = (iy - 1) * incy2;
            iys3 = (iy + 1) * incy2;
            iys4 = (iy - 2) * incy2;
            iys5 = (iy + 2) * incy2;

            for (iz = 2; iz < dimz + 2; iz++)
            {

                work[ixs + iys + iz - 2] = t0 * dum2[ixs1 + iys1 + iz] +
                    t1x * dum2[ixs2 + iys1 + iz] +
                    t1x * dum2[ixs3 + iys1 + iz] +
                    t1y * dum2[ixs1 + iys2 + iz] +
                    t1y * dum2[ixs1 + iys3 + iz] +
                    t1z * dum2[ixs1 + iys1 + iz - 1] +
                    t1z * dum2[ixs1 + iys1 + iz + 1] +
                    t2x * dum2[ixs4 + iys1 + iz] +
                    t2x * dum2[ixs5 + iys1 + iz] +
                    t2y * dum2[ixs1 + iys4 + iz] +
                    t2y * dum2[ixs1 + iys5 + iz] +
                    t2z * dum2[ixs1 + iys1 + iz - 2] + t2z * dum2[ixs1 + iys1 + iz + 2];

            }                   /* end for */
        }                       /* end for */
    }                           /* end for */

    my_free(dum2);


}                               /* end app6_del2 */
