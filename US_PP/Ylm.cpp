#include <math.h>
#include <stdio.h>
#include "const.h"
#include "transition.h"

double Ylm (int l, int m, double * r)
{
    double rmod, rhat[3], rr1, rr2, rr3, c, out;
    int i;

    rmod = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
    rmod = sqrt (rmod);

    int lm = l*l + m;
    rmod += 1.0e-9;
    for (i = 0; i < 3; i++)
        rhat[i] = r[i] / rmod;

    if (lm == 0)
    {
        c = sqrt (1.0 / fourPI);
        out = c;
        return (out);
    }
    if (lm == 1)
    {
        c = sqrt (3.0 / fourPI);
        out = c * rhat[0];
        return (out);
    }
    if (lm == 2)
    {
        c = sqrt (3.0 / fourPI);
        out = c * rhat[1];
        return (out);
    }
    if (lm == 3)
    {
        c = sqrt (3.0 / fourPI);
        out = c * rhat[2];
        return (out);
    }
    if (lm == 4)
    {
        c = sqrt (15.0 / fourPI);
        out = c * rhat[0] * rhat[1];
        return (out);
    }
    if (lm == 5)
    {
        c = sqrt (15.0 / fourPI);
        out = c * rhat[0] * rhat[2];
        return (out);
    }
    if (lm == 6)
    {
        c = sqrt (5.0 / (4.0 * fourPI));
        out = c * (3.0 * rhat[2] * rhat[2] - 1.0);
        if (rmod < 1.0e-8)
            out = 0.0;
        return (out);
    }
    if (lm == 7)
    {
        c = sqrt (15.0 / fourPI);
        out = c * rhat[1] * rhat[2];
        return (out);
    }
    if (lm == 8)
    {
        c = sqrt (15.0 / (4 * fourPI));
        out = c * (rhat[0] * rhat[0] - rhat[1] * rhat[1]);
        return (out);
    }
    if (lm == 9)
    {
        c = sqrt (7.0 / fourPI) * 5.0 / 2.0;
        out = c * rhat[0] * (rhat[0] * rhat[0] - 0.6);
        return (out);
    }
    if (lm == 10)
    {
        c = sqrt (7.0 / fourPI) * 5.0 / 2.0;
        out = c * rhat[1] * (rhat[1] * rhat[1] - 0.6);
        return (out);
    }
    if (lm == 11)
    {
        c = sqrt (7.0 * 15.0 / fourPI);
        out = c * rhat[0] * rhat[1] * rhat[2];
        return (out);
    }
    if (lm == 12)
    {
        c = sqrt (7.0 / fourPI) * 5.0 / 2.0;
        out = c * rhat[2] * (rhat[2] * rhat[2] - 0.6);
        return (out);
    }
    if (lm == 13)
    {
        c = sqrt (7.0 * 15.0 / fourPI) / 2.0;
        out = c * rhat[2] * (rhat[0] * rhat[0] - rhat[1] * rhat[1]);
        return (out);
    }
    if (lm == 14)
    {
        c = sqrt (7.0 * 15.0 / fourPI) / 2.0;
        out = c * rhat[1] * (rhat[2] * rhat[2] - rhat[0] * rhat[0]);
        return (out);
    }
    if (lm == 15)
    {
        c = sqrt (7.0 * 15.0 / fourPI) / 2.0;
        out = c * rhat[0] * (rhat[1] * rhat[1] - rhat[2] * rhat[2]);
        return (out);
    }
    if (lm == 16)
    {
        c = sqrt (3.0 * 7.0 / fourPI) * 5.0 / 4.0;
        rr1 = rhat[0] * rhat[0];
        rr2 = rhat[1] * rhat[1];
        rr3 = rhat[2] * rhat[2];
        out = c * (rr1 * rr1 + rr2 * rr2 + rr3 * rr3 - 0.6);
        if (rmod < 1.0e-8)
            out = 0.0;
        return (out);
    }
    if (lm == 17)
    {
        c = sqrt (9.0 * 35.0 / fourPI) / 2.0;
        out = c * rhat[1] * rhat[2] * (rhat[1] * rhat[1] - rhat[2] * rhat[2]);
        return (out);
    }
    if (lm == 18)
    {
        c = sqrt (9.0 * 35.0 / fourPI) / 2.0;
        out = c * rhat[0] * rhat[2] * (rhat[2] * rhat[2] - rhat[0] * rhat[0]);
        return (out);
    }
    if (lm == 19)
    {
        c = sqrt (9.0 * 5.0 / fourPI) / 4.0;
        rr1 = rhat[0] * rhat[0];
        rr2 = rhat[1] * rhat[1];
        rr3 = rhat[2] * rhat[2];
        out = c * (rr1 * rr1 - rr2 * rr2 - 6.0 * rr3 * (rr1 - rr2));
        return (out);
    }
    if (lm == 20)
    {
        c = sqrt (9.0 * 35.0 / fourPI) / 2.0;
        out = c * rhat[0] * rhat[1] * (rhat[0] * rhat[0] - rhat[1] * rhat[1]);
        return (out);
    }
    if (lm == 21)
    {
        c = sqrt (9.0 * 5.0 / fourPI) * 7.0 / 2.0;
        out = c * rhat[0] * rhat[1] * (rhat[2] * rhat[2] - 1.0 / 7.0);
        return (out);
    }
    if (lm == 22)
    {
        c = sqrt (9.0 * 5.0 / fourPI) * 7.0 / 2.0;
        out = c * rhat[0] * rhat[2] * (rhat[1] * rhat[1] - 1.0 / 7.0);
        return (out);
    }
    if (lm == 23)
    {
        c = sqrt (9.0 * 5.0 / fourPI) * 7.0 / 2.0;
        out = c * rhat[1] * rhat[2] * (rhat[0] * rhat[0] - 1.0 / 7.0);
        return (out);
    }
    if (lm == 24)
    {
        c = sqrt (9.0 * 5.0 / (3.0 * fourPI)) * 7.0 / 2.0;
        rr1 = rhat[0] * rhat[0];
        rr2 = rhat[1] * rhat[1];
        rr3 = rhat[2] * rhat[2];
        out =
            c * (rr3 * rr3 - 0.5 * (rr1 * rr1 + rr2 * rr2) - 6.0 / 7.0 * (rr3 - 0.5 * (rr1 + rr2)));
        return (out);
    }

    return CubicHarmonic(l, m, r);

}
