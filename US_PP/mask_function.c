/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include "main.h"
#include <math.h>

REAL mask_function (REAL x)
{
    if (x >= 1.0)
        return 0.0;

    if ((x > 0.0) && (x < 1.0))
        return (0.5*(cos(PI*x)+1));
        //return (15.8231*x*x*x*x*x -47.1834 *x*x*x*x + 50.0273 *x*x*x -20.3746*x*x + 0.7294*x + 0.9872);

    return 1.0;
}

