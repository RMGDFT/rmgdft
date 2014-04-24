#include "grid.h"
#include "common_prototypes.h"
#include "main.h"
#include <float.h>
#include <math.h>



void get_mgga_exc_vxc (rmg_double_t * rho, rmg_double_t * rho_oppo, rmg_double_t * tau, rmg_double_t * vxc, rmg_double_t * exc)
{



    mgga_libxc (rho, tau, vxc, exc, ct.xctype);



}                               /* end get_mgga_exc_vxc */



/******/
