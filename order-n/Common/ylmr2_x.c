/************************** SVN Revision Information **************************
 **    $Id: ylmr2_x.c 587 2006-08-18 22:57:56Z miro $    **
******************************************************************************/
 
#include <math.h>
#include <stdio.h>
#include "md.h"

void ylmr2_x(double *r, double *ylm_x)
{
    int l, i;
    REAL c, fpi, rmod, rhat[3];

    fpi = 4.0 * PI;
    rmod = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
    rmod = sqrt(rmod);

    rmod += 1.0e-9;
    for (i = 0; i < 3; i++)
        rhat[i] = r[i] / rmod;

    for (l = 0; l < 25; l++)
    {
        if (l == 0)
        {
            ylm_x[l] = 0;
        }
        else if (l == 1)
        {
            c = sqrt(3.0 / fpi);
            ylm_x[l] = c * (1.0 - rhat[0] * rhat[0]) / rmod;
            if (rmod < 1.0e-8)
                ylm_x[l] = 0.0;
        }
        else if (l == 2)
        {
            c = sqrt(3.0 / fpi);
            ylm_x[l] = -c * rhat[0] * rhat[2] / rmod;
        }
        else if (l == 3)
        {
            c = sqrt(3.0 / fpi);
            ylm_x[l] = -c * rhat[0] * rhat[1] / rmod;
        }
        else if (l == 4)
        {
            c = sqrt(15.0 / fpi);
            ylm_x[l] = c * rhat[1] * (1.0 - 2.0 * rhat[0] * rhat[0]) / rmod;
        }
        else if (l == 5)
        {
            c = sqrt(15.0 / fpi);
            ylm_x[l] = c * rhat[2] * (1.0 - 2.0 * rhat[0] * rhat[0]) / rmod;
        }
        else if (l == 6)
        {
            c = sqrt(5.0 / (4.0 * fpi));
            ylm_x[l] = c * (-6.0 * rhat[0] * rhat[2] * rhat[2]) / rmod;
        }
        else if (l == 7)
        {
            c = sqrt(15.0 / fpi);
            ylm_x[l] = c * (-2.0 * rhat[0] * rhat[1] * rhat[2]) / rmod;
        }
        else if (l == 8)
        {
            c = sqrt(15.0 / (4.0 * fpi));
            ylm_x[l] = c * 2.0 * rhat[0] * (1.0 - rhat[0] * rhat[0] + rhat[1] * rhat[1]) / rmod;
        }
        else if (l == 9)
        {
            c = sqrt(7.0 / fpi) * 5.0 / 2.0;
            ylm_x[l] =
                c * (3.0 * pow(rhat[0], 2.0) * (1.0 - rhat[0] * rhat[0]) -
                     0.6 * (1.0 - rhat[0] * rhat[0])) / rmod;
            if (rmod < 1.0e-8)
                ylm_x[l] = 0.0;
        }
        else if (l == 10)
        {
            c = sqrt(7.0 / fpi) * 5.0 / 2.0;
            ylm_x[l] = c * rhat[0] * rhat[1] * (-3.0 * rhat[1] * rhat[1] + 0.6) / rmod;
        }
        else if (l == 11)
        {
            c = sqrt(7.0 * 15.0 / fpi);
            ylm_x[l] = c * rhat[1] * rhat[2] * (1.0 - 3.0 * pow(rhat[0], 2.0)) / rmod;
        }
        else if (l == 12)
        {
            c = sqrt(7.0 / fpi) * 5.0 / 2.0;
            ylm_x[l] = c * rhat[0] * rhat[2] * (-3.0 * rhat[2] * rhat[2] + 0.6) / rmod;
        }
        else if (l == 13)
        {
            c = sqrt(7.0 * 15.0 / fpi) / 2.0;
            ylm_x[l] =
                c * rhat[0] * rhat[2] * (2.0 -
                                         3.0 * (rhat[0] * rhat[0] - rhat[1] * rhat[1])) / rmod;
        }
        else if (l == 14)
        {
            c = sqrt(7.0 * 15.0 / fpi) / 2.0;
            ylm_x[l] =
                c * rhat[0] * rhat[2] * (-2.0 -
                                         3.0 * (rhat[2] * rhat[2] - rhat[0] * rhat[0])) / rmod;
        }
        else if (l == 15)
        {
            c = sqrt(7.0 * 15.0 / fpi) / 2.0;
            ylm_x[l] =
                c * (rhat[1] * rhat[1] - rhat[2] * rhat[2]) * (1.0 -
                                                               3.0 * rhat[0] * rhat[0]) / rmod;
        }
        else if (l == 16)
        {
            c = sqrt(3.0 * 7.0 / fpi) * 5.0 / 4.0;
            ylm_x[l] =
                c * 4.0 * rhat[0] * (rhat[0] * rhat[0] - pow(rhat[0], 4.0) -
                                     pow(rhat[1], 4.0) - pow(rhat[2], 4.0)) / rmod;
        }
        else if (l == 17)
        {
            c = sqrt(9.0 * 35.0 / fpi) / 2.0;
            ylm_x[l] =
                c * (-4.0 * rhat[0] * rhat[1] * rhat[2]) * (rhat[1] *
                                                            rhat[1] - rhat[2] * rhat[2]) / rmod;
        }
        else if (l == 18)
        {
            c = sqrt(9.0 * 35.0 / fpi) / 2.0;
            ylm_x[l] =
                c * rhat[2] * ((rhat[2] * rhat[2] - 3.0 * rhat[0] * rhat[0]) -
                               4.0 * rhat[0] * rhat[0] * (rhat[2] * rhat[2] -
                                                          rhat[0] * rhat[0])) / rmod;
        }
        else if (l == 19)
        {
            c = sqrt(9.0 * 5.0 / fpi) / 4.0;
            ylm_x[l] =
                c * 4.0 * rhat[0] * (rhat[0] * rhat[0] -
                                     3.0 * rhat[2] * rhat[2] -
                                     (pow(rhat[0], 4.0) -
                                      pow(rhat[1],
                                          4.0) -
                                      6.0 * rhat[2] * rhat[2] * rhat[0] *
                                      rhat[0] +
                                      6.0 * rhat[2] * rhat[2] * rhat[1] * rhat[1])) / rmod;
        }
        else if (l == 20)
        {
            c = sqrt(9.0 * 35.0 / fpi) / 2.0;
            ylm_x[l] =
                c * (3.0 * rhat[0] * rhat[0] * rhat[1] - pow(rhat[1], 3.0) -
                     4.0 * (pow(rhat[0], 4.0) * rhat[1] -
                            rhat[0] * rhat[0] * pow(rhat[1], 3.0))) / rmod;
        }
        else if (l == 21)
        {
            c = sqrt(9.0 * 5.0 / fpi) * 7.0 / 2.0;
            ylm_x[l] =
                c * rhat[1] * (rhat[2] * rhat[2] - 1.0 / 7.0 -
                               4.0 * rhat[0] * rhat[0] * rhat[2] * rhat[2] +
                               2.0 / 7.0 * rhat[0] * rhat[0]) / rmod;
        }
        else if (l == 22)
        {
            c = sqrt(9.0 * 5.0 / fpi) * 7.0 / 2.0;
            ylm_x[l] =
                c * rhat[2] * (rhat[1] * rhat[1] - 1.0 / 7.0 -
                               4.0 * rhat[0] * rhat[0] * rhat[1] * rhat[1] +
                               2.0 / 7.0 * rhat[0] * rhat[0]) / rmod;
        }
        else if (l == 23)
        {
            c = sqrt(9.0 * 5.0 / fpi) * 7.0 / 2.0;
            ylm_x[l] =
                c * 2.0 * rhat[0] * rhat[1] * rhat[2] * (1.0 -
                                                         2.0 * rhat[0] *
                                                         rhat[0] + 1.0 / 7.0) / rmod;
        }
        else if (l == 24)
        {
            c = sqrt(9.0 * 5.0 / (3.0 * fpi)) * 7.0 / 2.0;
            ylm_x[l] =
                c * rhat[0] * (-4.0 * pow(rhat[2], 4.0) -
                               2.0 * (rhat[0] * rhat[0] - pow(rhat[0], 4.0) -
                                      pow(rhat[1],
                                          4.0)) -
                               6.0 / 7.0 * (-2.0 * rhat[2] * rhat[2] - 1.0 +
                                            rhat[0] * rhat[0] + rhat[1] * rhat[1])) / rmod;
        }
        else
            error_handler("higher l not programmed");
    }
}
