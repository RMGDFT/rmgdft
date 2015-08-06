#include "grid.h"
#include "common_prototypes.h"
#include "main.h"
#include <float.h>
#include <math.h>



void get_mgga_exc_vxc (double * rho, double * rho_oppo, double * tau, double * vxc, double * exc)
{



    mgga_libxc (rho, tau, vxc, exc, ct.xctype);



}                               /* end get_mgga_exc_vxc */



/******/
