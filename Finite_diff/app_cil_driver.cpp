#include "TradeImages.h"
#include "FiniteDiff.h"
#include "common_prototypes.h"
#include "rmg_alloc.h"
#include "rmg_error.h"


template <typename RmgType>
rmg_double_t CPP_app_cil_driver (RmgType * a, RmgType * b, int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz, int order)
{

    int P0_BASIS, numgrid, sbasis;
    rmg_double_t cc;
    void *allocp;
    RmgType *rptr;
    TradeImages T;

    P0_BASIS = get_P0_BASIS();
    numgrid = dimx * dimy * dimz;

    sbasis = (dimx + 4) * (dimy + 4) * (dimz + 4);
    my_malloc (allocp, sbasis + 64, double);
    rptr = (RmgType *)allocp;

    if(order == APP_CI_FOURTH) {
        T.CPP_trade_imagesx (a, rptr, dimx, dimy, dimz, 1, FULL_FD);
        if(numgrid == P0_BASIS) {
            cc = FD_app_cil_fourth_global (rptr, b, gridhx, gridhy, gridhz);
        }
        else {
            cc = FD_app_cil_fourth_standard (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);
        }
        my_free(rptr);
        return cc;
    }
    if(order == APP_CI_SIXTH) {
        T.CPP_trade_imagesx (a, rptr, dimx, dimy, dimz, 2, FULL_FD);
        if(numgrid == P0_BASIS) {
            cc = FD_app_cil_sixth_global (rptr, b, gridhx, gridhy, gridhz);
        }
        else {
            cc = FD_app_cil_sixth_standard (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);
        }

        my_free(rptr);
        return cc;
    }

    rmg_error_handler("APP_CIL order not programmed yet in app_cil_driver.\n");

}

extern "C" rmg_double_t app_cil_driver (rmg_double_t * a, rmg_double_t * b, int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz, int order)
{
    CPP_app_cil_driver<double>(a, b, dimx, dimy, dimz, gridhx, gridhy, gridhz, order);
}
extern "C" rmg_double_t app_cil_driver_f (rmg_float_t * a, rmg_float_t * b, int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz, int order)
{
    CPP_app_cil_driver<float>(a, b, dimx, dimy, dimz, gridhx, gridhy, gridhz, order);
}
