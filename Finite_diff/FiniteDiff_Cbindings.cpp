#include "TradeImages.h"
#include "Lattice.h"
#include "FiniteDiff.h"
#include "RmgTimer.h"
#include "common_prototypes.h"
#include "rmg_alloc.h"
#include "rmg_error.h"


extern "C" void app_cir_driver (double * a, double * b, int dimx, int dimy, int dimz, int order)
{
    TradeImages T;
    Lattice L;
    CPP_app_cir_driver<double>(&L, &T, a, b, dimx, dimy, dimz, order);
}
extern "C" void app_cir_driver_f (rmg_float_t * a, rmg_float_t * b, int dimx, int dimy, int dimz, int order)
{
    TradeImages T;
    Lattice L;
    CPP_app_cir_driver<float>(&L, &T, a, b, dimx, dimy, dimz, order);
}
extern "C" double app_cil_driver (double * a, double * b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz, int order)
{
    TradeImages T;
    Lattice L;
    return CPP_app_cil_driver<double>(&L, &T, a, b, dimx, dimy, dimz, gridhx, gridhy, gridhz, order);
}
extern "C" double app_cil_driver_f (rmg_float_t * a, rmg_float_t * b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz, int order)
{
    TradeImages T;
    Lattice L;
    return CPP_app_cil_driver<float>(&L, &T, a, b, dimx, dimy, dimz, gridhx, gridhy, gridhz, order);
}

