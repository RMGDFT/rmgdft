
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

for(int ix=0;ix<ixx;ix++){
  for(int iy=0;iy<iyy;iy++){
    for(int iz=0;iz<izz;iz++){
            if((ix==0)||(iy==0)||(iz==0))a[ix*iyy*izz+iy*izz+iz] = 0.0;
            if((ix==(ixx-1))||(iy==(iyy-1))||(iz==(izz-1)))a[ix*iyy*izz+iy*izz+iz] = 0.0;
    }
  }
}
#if 0
for(int ix=0;ix<ixx;ix++){
  for(int iy=0;iy<iyy;iy++){
    for(int iz=0;iz<izz;iz++){
        if(width == 1){
            if((ix==0)||(iy==0)||(iz==0))a[ix*iyy*izz+iy*izz+iz] = 0.0;
            if((ix==(ixx-1))||(iy==(iyy-1))||(iz==(izz-1)))a[ix*iyy*izz+iy*izz+iz] = 0.0;
        }

        if(width == 2){
            if((ix==0)||(iy==0)||(iz==0))a[ix*iyy*izz+iy*izz+iz] = 0.0;
            if((ix==(ixx-1))||(iy==(iyy-1))||(iz==(izz-1)))a[ix*iyy*izz+iy*izz+iz] = 0.0;
            if((ix==1)||(iy==1)||(iz==1))a[ix*iyy*izz+iy*izz+iz] *= 0.5;
            if((ix==(ixx-2))||(iy==(iyy-2))||(iz==(izz-2)))a[ix*iyy*izz+iy*izz+iz] *= 0.5;
        }

        if(width == 3){
            if((ix==0)||(iy==0)||(iz==0))a[ix*iyy*izz+iy*izz+iz] = 0.0;
            if((ix==(ixx-1))||(iy==(iyy-1))||(iz==(izz-1)))a[ix*iyy*izz+iy*izz+iz] = 0.0;
            if((ix==1)||(iy==1)||(iz==1))a[ix*iyy*izz+iy*izz+iz] *= 0.3333;
            if((ix==(ixx-2))||(iy==(iyy-2))||(iz==(izz-2)))a[ix*iyy*izz+iy*izz+iz] *= 0.3333;
            if((ix==2)||(iy==2)||(iz==2))a[ix*iyy*izz+iy*izz+iz] *= 0.6666;
            if((ix==(ixx-3))||(iy==(iyy-3))||(iz==(izz-3)))a[ix*iyy*izz+iy*izz+iz] *= 0.6666;
        }

        if(width == 4){
            if((ix==0)||(iy==0)||(iz==0))a[ix*iyy*izz+iy*izz+iz] = 0.0;
            if((ix==(ixx-1))||(iy==(iyy-1))||(iz==(izz-1)))a[ix*iyy*izz+iy*izz+iz] = 0.0;
            if((ix==1)||(iy==1)||(iz==1))a[ix*iyy*izz+iy*izz+iz] *= 0.0;
            if((ix==(ixx-2))||(iy==(iyy-2))||(iz==(izz-2)))a[ix*iyy*izz+iy*izz+iz] *= 0.0;
            if((ix==2)||(iy==2)||(iz==2))a[ix*iyy*izz+iy*izz+iz] *= 0.0;
            if((ix==(ixx-3))||(iy==(iyy-3))||(iz==(izz-3)))a[ix*iyy*izz+iy*izz+iz] *= 0.0;
            if((ix==3)||(iy==3)||(iz==3))a[ix*iyy*izz+iy*izz+iz] *= 0.0;
            if((ix==(ixx-4))||(iy==(iyy-4))||(iz==(izz-4)))a[ix*iyy*izz+iy*izz+iz] *= 0.0;
            if((ix==4)||(iy==4)||(iz==4))a[ix*iyy*izz+iy*izz+iz] *= 0.0;
            if((ix==(ixx-5))||(iy==(iyy-5))||(iz==(izz-5)))a[ix*iyy*izz+iy*izz+iz] *= 0.0;
        }

        if(width == 5){
            if((ix==0)||(iy==0)||(iz==0))a[ix*iyy*izz+iy*izz+iz] = 0.0;
            if((ix==(ixx-1))||(iy==(iyy-1))||(iz==(izz-1)))a[ix*iyy*izz+iy*izz+iz] = 0.0;
            if((ix==1)||(iy==1)||(iz==1))a[ix*iyy*izz+iy*izz+iz] *= 0.2;
            if((ix==(ixx-2))||(iy==(iyy-2))||(iz==(izz-2)))a[ix*iyy*izz+iy*izz+iz] *= 0.2;
            if((ix==2)||(iy==2)||(iz==2))a[ix*iyy*izz+iy*izz+iz] *= 0.4;
            if((ix==(ixx-3))||(iy==(iyy-3))||(iz==(izz-3)))a[ix*iyy*izz+iy*izz+iz] *= 0.4;
            if((ix==3)||(iy==3)||(iz==3))a[ix*iyy*izz+iy*izz+iz] *= 0.6;
            if((ix==(ixx-4))||(iy==(iyy-4))||(iz==(izz-4)))a[ix*iyy*izz+iy*izz+iz] *= 0.6;
            if((ix==4)||(iy==4)||(iz==4))a[ix*iyy*izz+iy*izz+iz] *= 0.8;
            if((ix==(ixx-5))||(iy==(iyy-5))||(iz==(izz-5)))a[ix*iyy*izz+iy*izz+iz] *= 0.8;
        }


#if 0
        if(width >= 5){
            if((ix==4)||(iy==4)||(iz==4))a[ix*iyy*izz+iy*izz+iz] = 0.0;
            if((ix==(ixx-5))||(iy==(iyy-5))||(iz==(izz-5)))a[ix*iyy*izz+iy*izz+iz] = 0.0;
        }
#endif
    }
  }
}
#endif

}
