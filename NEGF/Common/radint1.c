/************************** SVN Revision Information **************************
 **    $Id: radint1.c 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 
#include <math.h>
#include "md.h"
#include <float.h>

REAL radint1 (REAL * func, REAL * r, REAL * rab, int n)
{
    REAL sum, r12, f1, f2, f3;
    int i;

    sum = 0.0;
    r12 = 1.0 / 12.0;
/*	f3=func[0]*rab[0]*r12;*/
    f3 = 0.0;

    for (i = 0; i < n - 1; i = i + 2)
    {
        f1 = f3;
        f2 = func[i] * rab[i] * r[i] * r[i] * r12;
        f3 = func[i + 1] * rab[i + 1] * r[i + 1] * r[i + 1] * r12;
        sum = sum + 4.0 * f1 + 16.0 * f2 + 4.0 * f3;
    }
    return sum;
}                               /*end radint1 */

REAL radint2 (REAL * func, REAL * r, REAL * rab, int n)
{
    REAL sum, r12, f1, f2, f3;
    int i;

    sum = 0.0;
    r12 = 1.0 / 12.0;
/*	f3=func[0]*rab[0]*r12;*/
    f3 = 0.0;

    for (i = 0; i < n - 1; i = i + 2)
    {
        f1 = f3;
        f2 = func[i] * rab[i] * r12;
        f3 = func[i + 1] * rab[i + 1] * r12;
        sum = sum + 4.0 * f1 + 16.0 * f2 + 4.0 * f3;
    }
    return sum;
}                               /*end radint1 */


/******/
