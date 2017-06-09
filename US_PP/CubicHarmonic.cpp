/************************** SVN Revision Information **************************
 **    $Id: ylm.c 3141 2015-08-07 15:09:38Z luw $    **
******************************************************************************/

#include <math.h>
#include <stdio.h>
#include "const.h"
#include <boost/math/special_functions/spherical_harmonic.hpp>

using boost::math::policies::policy;
using boost::math::policies::promote_double;
typedef policy<promote_double<false> > harmonic_policy;

double CubicHarmonic(int L, int M, double *r)
{

//  L:  momentum
//  M = [0, 2*L+1)
//  spherical m = M-L,  range [-L,L]

    double xymod, rmod, theta, phi;
    double val;
    double eps = 1.0e-10;
    if(L == 0) return sqrt(1.0/fourPI);

    rmod = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
    rmod = sqrt (rmod);

    if(rmod < eps ) return 0.0;
    theta = acos(r[2]/rmod);

    if(M == L)
    {
        phi = 0.0;
        val  =  std::real(boost::math::spherical_harmonic(L, M-L, theta, phi, harmonic_policy()));
        return val;
    }

    xymod = r[0] * r[0] + r[1] * r[1];

    if(xymod < eps *eps)
    {
            return 0.0;
    }

    phi = atan2(r[1], r[0]);
    if (theta < 0.0) theta += PI;
    if(phi < 0.0) phi += twoPI;

    if( M < L)
        val =  std::real(boost::math::spherical_harmonic(L, M-L, theta, phi, harmonic_policy()));
    else
    {
        val =  std::imag(boost::math::spherical_harmonic(L, M-L, theta, phi, harmonic_policy()));
        if( (M-L)%2) val = -val;
    }

    val *= sqrt(2.0);

    return val;

}

