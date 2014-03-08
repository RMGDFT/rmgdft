#include "TradeImages.h"
#include "FiniteDiff.h"
#include "common_prototypes.h"
#include "rmg_alloc.h"
#include "rmg_error.h"
#include "RmgTimer.h"

using namespace std;

template <typename RmgType>
rmg_double_t CPP_app_cil_driver (RmgType * a, RmgType * b, int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz, int order)
{

    RmgTimer RT("App_cil total time");
    int sbasis;
    rmg_double_t cc;
    void *allocp;
    RmgType *rptr;
    TradeImages T;
    FiniteDiff FD;
    sbasis = (dimx + 4) * (dimy + 4) * (dimz + 4);
    my_malloc (allocp, sbasis + 64, double);
    rptr = (RmgType *)allocp;

    if(order == APP_CI_FOURTH) {
        RmgTimer *RT1 = new RmgTimer("App_cil trade images time");
        T.trade_imagesx (a, rptr, dimx, dimy, dimz, 1, FULL_TRADE);
        delete(RT1);
        cc = FD.app_cil_fourth (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);
    }
    else if(order == APP_CI_SIXTH) {
        RmgTimer *RT1 = new RmgTimer("App_cil trade images time");
        T.trade_imagesx (a, rptr, dimx, dimy, dimz, 2, FULL_TRADE);
        delete(RT1);
        cc = FD.app_cil_sixth (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);
    }
    else {
        rmg_error_handler (__FILE__, __LINE__, "APP_CIL order not programmed yet in app_cil_driver.\n");
    }

    my_free(rptr);
    return cc;

}

extern "C" rmg_double_t app_cil_driver (rmg_double_t * a, rmg_double_t * b, int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz, int order)
{
    return CPP_app_cil_driver<double>(a, b, dimx, dimy, dimz, gridhx, gridhy, gridhz, order);
}
extern "C" rmg_double_t app_cil_driver_f (rmg_float_t * a, rmg_float_t * b, int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz, int order)
{
    return CPP_app_cil_driver<float>(a, b, dimx, dimy, dimz, gridhx, gridhy, gridhz, order);
}
