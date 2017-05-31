#include "make_conf.h"
#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "prototypes_on.h"

extern "C" void ZeroBoundary_C(double *a, int ixx, int iyy, int izz)
{
    ZeroBoundary(a, ixx, iyy, izz);
}

void ZeroBoundary(double *a, int ixx, int iyy, int izz)
{

    for(int ix=0;ix<ixx;ix+=ixx-1)
    {
        for(int iy=0;iy<iyy;iy++)
        {
            for(int iz=0;iz<izz;iz++)
            {
                a[ix*iyy*izz+iy*izz+iz] = 0.0;
            }
        }
    }
    for(int ix=0;ix<ixx;ix++)
    {
        for(int iy=0;iy<iyy;iy+=iyy-1)
        {
            for(int iz=0;iz<izz;iz++)
            {
                a[ix*iyy*izz+iy*izz+iz] = 0.0;
            }
        }
    }
    for(int ix=0;ix<ixx;ix++)
    {
        for(int iy=0;iy<iyy;iy++)
        {
            for(int iz=0;iz<izz;iz+=izz-1)
            {
                a[ix*iyy*izz+iy*izz+iz] = 0.0;
            }
        }
    }

}
