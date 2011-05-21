/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <math.h>
#include <stdio.h>
#include "main.h"

void ylmr2 (double *r, double *out)
{
    int l, i;
    double rmod, rhat[3], rr1, rr2, rr3, c;

    rmod = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
    rmod = sqrt (rmod);

    rmod += 1.0e-9;
    for (i = 0; i < 3; i++)
        rhat[i] = r[i] / rmod;

    for (l = 0; l < 25; l++)
    {
        if (l == 0)
        {
            c = sqrt (1.0 / fourPI);
            out[l] = c;
        }
        else if (l == 1)
        {
            c = sqrt (3.0 / fourPI);
            out[l] = c * rhat[0];
        }
        else if (l == 2)
        {
            c = sqrt (3.0 / fourPI);
            out[l] = c * rhat[2];
        }
        else if (l == 3)
        {
            c = sqrt (3.0 / fourPI);
            out[l] = c * rhat[1];
        }
        else if (l == 4)
        {
            c = sqrt (15.0 / fourPI);
            out[l] = c * rhat[0] * rhat[1];
        }
        else if (l == 5)
        {
            c = sqrt (15.0 / fourPI);
            out[l] = c * rhat[0] * rhat[2];
        }
        else if (l == 6)
        {
            c = sqrt (5.0 / (4.0 * fourPI));
            out[l] = c * (3.0 * rhat[2] * rhat[2] - 1.0);
            if (rmod < 1.0e-8)
                out[l] = 0.0;
        }
        else if (l == 7)
        {
            c = sqrt (15.0 / fourPI);
            out[l] = c * rhat[1] * rhat[2];
        }
        else if (l == 8)
        {
            c = sqrt (15.0 / (4 * fourPI));
            out[l] = c * (rhat[0] * rhat[0] - rhat[1] * rhat[1]);
        }
        else if (l == 9)
        {
            c = sqrt (7.0 / fourPI) * 5.0 / 2.0;
            out[l] = c * rhat[0] * (rhat[0] * rhat[0] - 0.6);
        }
        else if (l == 10)
        {
            c = sqrt (7.0 / fourPI) * 5.0 / 2.0;
            out[l] = c * rhat[1] * (rhat[1] * rhat[1] - 0.6);
        }
        else if (l == 11)
        {
            c = sqrt (7.0 * 15.0 / fourPI);
            out[l] = c * rhat[0] * rhat[1] * rhat[2];
        }
        else if (l == 12)
        {
            c = sqrt (7.0 / fourPI) * 5.0 / 2.0;
            out[l] = c * rhat[2] * (rhat[2] * rhat[2] - 0.6);
        }
        else if (l == 13)
        {
            c = sqrt (7.0 * 15.0 / fourPI) / 2.0;
            out[l] = c * rhat[2] * (rhat[0] * rhat[0] - rhat[1] * rhat[1]);
        }
        else if (l == 14)
        {
            c = sqrt (7.0 * 15.0 / fourPI) / 2.0;
            out[l] = c * rhat[1] * (rhat[2] * rhat[2] - rhat[0] * rhat[0]);
        }
        else if (l == 15)
        {
            c = sqrt (7.0 * 15.0 / fourPI) / 2.0;
            out[l] = c * rhat[0] * (rhat[1] * rhat[1] - rhat[2] * rhat[2]);
        }
        else if (l == 16)
        {
            c = sqrt (3.0 * 7.0 / fourPI) * 5.0 / 4.0;
            rr1 = rhat[0] * rhat[0];
            rr2 = rhat[1] * rhat[1];
            rr3 = rhat[2] * rhat[2];
            out[l] = c * (rr1 * rr1 + rr2 * rr2 + rr3 * rr3 - 0.6);
            if (rmod < 1.0e-8)
                out[l] = 0.0;
        }
        else if (l == 17)
        {
            c = sqrt (9.0 * 35.0 / fourPI) / 2.0;
            out[l] = c * rhat[1] * rhat[2] * (rhat[1] * rhat[1] - rhat[2] * rhat[2]);
        }
        else if (l == 18)
        {
            c = sqrt (9.0 * 35.0 / fourPI) / 2.0;
            out[l] = c * rhat[0] * rhat[2] * (rhat[2] * rhat[2] - rhat[0] * rhat[0]);
        }
        else if (l == 19)
        {
            c = sqrt (9.0 * 5.0 / fourPI) / 4.0;
            rr1 = rhat[0] * rhat[0];
            rr2 = rhat[1] * rhat[1];
            rr3 = rhat[2] * rhat[2];
            out[l] = c * (rr1 * rr1 - rr2 * rr2 - 6.0 * rr3 * (rr1 - rr2));
        }
        else if (l == 20)
        {
            c = sqrt (9.0 * 35.0 / fourPI) / 2.0;
            out[l] = c * rhat[0] * rhat[1] * (rhat[0] * rhat[0] - rhat[1] * rhat[1]);
        }
        else if (l == 21)
        {
            c = sqrt (9.0 * 5.0 / fourPI) * 7.0 / 2.0;
            out[l] = c * rhat[0] * rhat[1] * (rhat[2] * rhat[2] - 1.0 / 7.0);
        }
        else if (l == 22)
        {
            c = sqrt (9.0 * 5.0 / fourPI) * 7.0 / 2.0;
            out[l] = c * rhat[0] * rhat[2] * (rhat[1] * rhat[1] - 1.0 / 7.0);
        }
        else if (l == 23)
        {
            c = sqrt (9.0 * 5.0 / fourPI) * 7.0 / 2.0;
            out[l] = c * rhat[1] * rhat[2] * (rhat[0] * rhat[0] - 1.0 / 7.0);
        }
        else if (l == 24)
        {
            c = sqrt (9.0 * 5.0 / (3.0 * fourPI)) * 7.0 / 2.0;
            rr1 = rhat[0] * rhat[0];
            rr2 = rhat[1] * rhat[1];
            rr3 = rhat[2] * rhat[2];
            out[l] =
                c * (rr3 * rr3 - 0.5 * (rr1 * rr1 + rr2 * rr2) -
                     6.0 / 7.0 * (rr3 - 0.5 * (rr1 + rr2)));
        }
        else
            error_handler ("higher l not programmed");
    }
}
