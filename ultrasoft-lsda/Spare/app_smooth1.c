/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include "main.h"
#include <float.h>
#include <math.h>


/**

 Smooths f and returns result in work
*/
void app_smooth1 (FS0_GRID * f, FS0_GRID * work)
{

    int iz, ix, iy;
    REAL scale, ec, fc, crn, cc;

    cc = 8.0;
    fc = 8.0;
    ec = 8.0;
    crn = 8.0;
    scale = 1.0 / 216.0;


    for (ix = 1; ix <= FPX0_GRID; ix++)
    {

        for (iy = 1; iy <= FPY0_GRID; iy++)
        {

            for (iz = 1; iz <= FPZ0_GRID; iz++)
            {


                work->s1.b[ix][iy][iz] = cc * f->s1.b[ix][iy][iz] +
                    fc * f->s1.b[ix - 1][iy][iz] +
                    fc * f->s1.b[ix + 1][iy][iz] +
                    fc * f->s1.b[ix][iy - 1][iz] +
                    fc * f->s1.b[ix][iy + 1][iz] +
                    fc * f->s1.b[ix][iy][iz - 1] + fc * f->s1.b[ix][iy][iz + 1];

                work->s1.b[ix][iy][iz] +=
                    ec * f->s1.b[ix - 1][iy][iz - 1] +
                    ec * f->s1.b[ix + 1][iy][iz - 1] +
                    ec * f->s1.b[ix][iy - 1][iz - 1] +
                    ec * f->s1.b[ix][iy + 1][iz - 1] +
                    ec * f->s1.b[ix - 1][iy - 1][iz] +
                    ec * f->s1.b[ix - 1][iy + 1][iz] +
                    ec * f->s1.b[ix + 1][iy - 1][iz] +
                    ec * f->s1.b[ix + 1][iy + 1][iz] +
                    ec * f->s1.b[ix - 1][iy][iz + 1] +
                    ec * f->s1.b[ix + 1][iy][iz + 1] +
                    ec * f->s1.b[ix][iy - 1][iz + 1] + ec * f->s1.b[ix][iy + 1][iz + 1];

                work->s1.b[ix][iy][iz] +=
                    crn * f->s1.b[ix - 1][iy - 1][iz - 1] +
                    crn * f->s1.b[ix - 1][iy - 1][iz + 1] +
                    crn * f->s1.b[ix - 1][iy + 1][iz - 1] +
                    crn * f->s1.b[ix - 1][iy + 1][iz + 1] +
                    crn * f->s1.b[ix + 1][iy - 1][iz - 1] +
                    crn * f->s1.b[ix + 1][iy - 1][iz + 1] +
                    crn * f->s1.b[ix + 1][iy + 1][iz - 1] + crn * f->s1.b[ix + 1][iy + 1][iz + 1];

                work->s1.b[ix][iy][iz] = scale * work->s1.b[ix][iy][iz];


            }                   /* end for */

        }                       /* end for */

    }                           /* end for */



}                               /* end app_smooth */






/******/
