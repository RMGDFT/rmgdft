
#include "main.h"
#include <float.h>
#include <math.h> 

#define SMALL 1.e-10

void gcxbcp (REAL rho, REAL grad, REAL * enxc, REAL * vxc1, REAL * vxc2)
{
	REAL sx, sc, v1x, v2x, v1c, v2c, arho;


	sx = v1x = v2x = 0.0;
	sc = v1c = v2c = 0.0;
	
	arho = fabs(rho);
	
	/* exchange correction */
	if (arho > SMALL )
		becke88 (arho, grad, &sx, &v1x, &v2x); 

	/* correlation correction */
	if (arho > SMALL )
		perdew86 (arho, grad, &sc, &v1c, &v2c);

	*vxc1 = v1x + v1c;
	*vxc2 = v2x + v2c;
	*enxc = sx + sc;
}









































































































