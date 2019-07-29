#include <math.h>
#include <stdio.h>
#include "main.h"

double get_QnmL (int idx, int ltot, double r, SPECIES * sp)
{
    int i;
    double return_value, r2, power_r;
    double *coeff;

    if(ct.norm_conserving_pp) return 0.0;

    coeff = sp->qfcoef + idx * sp->nlc * sp->nqf + ltot * sp->nqf;

    /*Special treatment for ltot == 0 and small r */
    if ((ltot == 0) && (r < 1e-12))
        return (coeff[0]);


    r2 = r * r;
    return_value = 0.0;
    power_r = ONE;

    for (i = 0; i < sp->nqf; i++)
    {

        return_value += coeff[i] * power_r;

        power_r *= r2;

#if 0
        if (power_r > 1e-10)
            power_r *= r2;
        else
            break;
#endif

    }

    return_value *= pow (r, (double) ltot);

    return (return_value);
}
