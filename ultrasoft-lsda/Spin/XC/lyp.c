

#include "main.h"
#include <float.h>
#include <math.h> 


void lyp ( REAL rs, REAL * ec, REAL * vc )
{
	/* C. Lee. W. Yang and R.G. Parr, PRB 37, 785 (1988 ) - LDA part only */
	
	REAL a=0.04918, b=0.132*2.87123400018819108;
	REAL pi43=1.61199195401647, c=0.2533*pi43, d=0.349*pi43;
	/* pi43 = (4 * pi / 3)^(1 / 3) */

	REAL ecrs, ox;

	ecrs = b * exp ( - c * rs );
	ox = 1.0 / ( 1.0 + d * rs );
	*ec = - a * ox * ( 1.0 + ecrs );
	*vc = (*ec) - rs / 3.0 * a * ox * ( d * ox + ecrs * ( d * ox + c ) );
}


void glyp ( REAL rho, REAL grho, REAL * sc, REAL * v1c, REAL * v2c )
{
	/* Lee Yang Parr: gradient correction part */

	REAL a=0.04918, b=0.132, c=0.2533, d=0.349;
	REAL rhom13, rhom43, rhom53, om, xl, ff, dom, dxl;

	rhom13 = pow ( rho, - 1.0 / 3.0 );
	om = exp ( -c * rhom13 ) / ( 1.0 + d * rhom13 );
	xl = 1.0 + ( 7.0 / 3.0 ) * ( c * rhom13 + d * rhom13 / ( 1.0 + d * rhom13 ) );
	ff = a * b * grho * grho / 24.0;
	rhom53 = pow ( rhom13, 5 );
	*sc = ff * rhom53 * om * xl;
	dom = - om * ( c + d + c * d * rhom13 ) / ( 1.0 + d * rhom13 );
	dxl = ( 7.0 / 3.0 ) * ( c + d + 2.0 * c * d * rhom13 + c * d * d * pow ( rhom13, 2 ) ) / pow ( ( 1.0 + d * rhom13 ), 2 );
	rhom43 = pow ( rhom13, 4 );
	*v1c = - ff * rhom43 / 3.0 * ( 5.0 * rhom43 * om * xl + rhom53 * dom * xl + rhom53 * om * dxl );
	*v2c = 2.0 * (*sc) / ( grho * grho );
	*sc = (*sc) / rho;
}
