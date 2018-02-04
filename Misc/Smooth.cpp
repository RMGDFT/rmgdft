#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
//#include "main.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "transition.h"
#include "packfuncs.h"
#include "math.h"



void Smooth(double *in, double *out, int dimx, int dimy, int dimz, double q)
{

    double factor = q*q*q / (q*q*q + 6.0*q*q + 12.0*q + 8.0);
    double q000 = 1.0;
    double q100 = 1.0 / q;
    double q110 = 1.0 / (q * q);
    double q111 = 1.0 / (q * q * q);

    int incx = (dimy + 2)*(dimz + 2);
    int incy = dimz + 2;
    double *buf = new double[(dimx+2)*(dimy+2)*(dimz+2)]();
    Rmg_T->trade_imagesx (in, buf, dimx, dimy, dimz, 1, FULL_TRADE);


    for(int ix=1;ix <= dimx;ix++)
    {
        for(int iy=1;iy <= dimy;iy++)
        {
            for(int iz=1;iz <= dimz;iz++)
            {
                out[(ix-1)*dimy*dimz + (iy-1)*dimz + iz - 1] =
                q111*buf[(ix-1)*incx + (iy-1)*incy + iz-1] +
                q110*buf[(ix-1)*incx + (iy-1)*incy + iz] +
                q111*buf[(ix-1)*incx + (iy-1)*incy + iz + 1] +

                q110*buf[(ix-1)*incx + iy*incy + iz-1] +
                q100*buf[(ix-1)*incx + iy*incy + iz] +
                q110*buf[(ix-1)*incx + iy*incy + iz + 1] +

                q111*buf[(ix-1)*incx + (iy+1)*incy + iz-1] +
                q110*buf[(ix-1)*incx + (iy+1)*incy + iz] +
                q111*buf[(ix-1)*incx + (iy+1)*incy + iz + 1] +

                q110*buf[ix*incx + (iy-1)*incy + iz-1] +
                q100*buf[ix*incx + (iy-1)*incy + iz] +
                q110*buf[ix*incx + (iy-1)*incy + iz + 1] +

                q100*buf[ix*incx + iy*incy + iz-1] +
                q000*buf[ix*incx + iy*incy + iz] +
                q100*buf[ix*incx + iy*incy + iz + 1] +

                q110*buf[ix*incx + (iy+1)*incy + iz-1] +
                q100*buf[ix*incx + (iy+1)*incy + iz] +
                q110*buf[ix*incx + (iy+1)*incy + iz + 1] +

                q111*buf[(ix+1)*incx + (iy-1)*incy + iz-1] +
                q110*buf[(ix+1)*incx + (iy-1)*incy + iz] +
                q111*buf[(ix+1)*incx + (iy-1)*incy + iz + 1] +

                q110*buf[(ix+1)*incx + iy*incy + iz-1] +
                q100*buf[(ix+1)*incx + iy*incy + iz] +
                q110*buf[(ix+1)*incx + iy*incy + iz + 1] +

                q111*buf[(ix+1)*incx + (iy+1)*incy + iz-1] +
                q110*buf[(ix+1)*incx + (iy+1)*incy + iz] +
                q111*buf[(ix+1)*incx + (iy+1)*incy + iz + 1];

                out[(ix-1)*dimy*dimz + (iy-1)*dimz + iz - 1] *= factor;
            }
        }
    }

    delete [] buf;

}
