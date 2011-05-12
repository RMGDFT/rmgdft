/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <math.h>
#include <stdio.h>
#include "md.h"

REAL get_QnmL(int idx, int ltot, REAL r, SPECIES * sp)
{
    int i;
    REAL return_value, r2, power_r;
    REAL *coeff;


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

    return_value *= pow(r, (REAL) ltot);

    return (return_value);
}
