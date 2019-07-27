#include "TradeImages.h"
#include "Lattice.h"
#include "FiniteDiff.h"
#include "common_prototypes.h"
#include "transition.h"

void app_cir_driver (double * a, double * b, int dimx, int dimy, int dimz, int order)
{
    CPP_app_cir_driver<double>(&Rmg_L, Rmg_T, a, b, dimx, dimy, dimz, order);
}
double app_cil_driver (double * a, double * b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz, int order)
{
    return CPP_app_cil_driver<double>(&Rmg_L, Rmg_T, a, b, dimx, dimy, dimz, gridhx, gridhy, gridhz, order);
}

