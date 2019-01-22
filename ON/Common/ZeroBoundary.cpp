
#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "prototypes_on.h"

extern "C" void ZeroBoundaryC(double *a, int ixx, int iyy, int izz)
{
    ZeroBoundary(a, ixx, iyy, izz);
}

void ZeroBoundary(double *a, int ixx, int iyy, int izz)
{
    ZeroBoundary(a, ixx, iyy, izz, ct.kohn_sham_fd_order / 2);
}

void ZeroBoundary(double *a, int ixx, int iyy, int izz, int width)
{

    for(int ix=0;ix<width;ix++){
        for(int iy=0;iy<iyy;iy++){
            for(int iz=0;iz<izz;iz++){
                a[ix*iyy*izz+iy*izz+iz] = 0.0;
            }
        }
    }

    for(int ix=ixx-width;ix<ixx;ix++){
        for(int iy=0;iy<iyy;iy++){
            for(int iz=0;iz<izz;iz++){
                a[ix*iyy*izz+iy*izz+iz] = 0.0;
            }
        }
    }

    for(int ix=0;ix<ixx;ix++){
        for(int iy=0;iy<width;iy++){
            for(int iz=0;iz<izz;iz++){
                a[ix*iyy*izz+iy*izz+iz] = 0.0;
            }
        }
    }
    for(int ix=0;ix<ixx;ix++){
        for(int iy=iyy-width;iy<iyy;iy++){
            for(int iz=0;iz<izz;iz++){
                a[ix*iyy*izz+iy*izz+iz] = 0.0;
            }
        }
    }
    for(int ix=0;ix<ixx;ix++){
        for(int iy=0;iy<iyy;iy++){
            for(int iz=0;iz<width;iz++){
                a[ix*iyy*izz+iy*izz+iz] = 0.0;
            }
        }
    }
    for(int ix=0;ix<ixx;ix++){
        for(int iy=0;iy<iyy;iy++){
            for(int iz=izz-width;iz<izz;iz++){
                a[ix*iyy*izz+iy*izz+iz] = 0.0;
            }
        }
    }
}
