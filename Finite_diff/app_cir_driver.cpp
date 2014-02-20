#include "TradeImages.h"
#include "FiniteDiff.h"
#include "common_prototypes.h"
#include "rmg_alloc.h"
#include "rmg_error.h"


template <typename RmgType>
void CPP_app_cir_driver (RmgType * a, RmgType * b, int dimx, int dimy, int dimz, int order)
{

    int P0_BASIS, numgrid, sbasis;
    void *allocp;
    RmgType *rptr;
    TradeImages T;

    P0_BASIS = get_P0_BASIS();
    numgrid = dimx * dimy * dimz;

    sbasis = (dimx + 4) * (dimy + 4) * (dimz + 4);
    my_malloc (allocp, sbasis + 64, double);
    rptr = (RmgType *)allocp;

    if(order == APP_CI_FOURTH) {
        T.CPP_trade_imagesx (a, rptr, dimx, dimy, dimz, 1, CENTRAL_FD);
        if(numgrid == P0_BASIS) {
            FD_app_cir_fourth_global (rptr, b);
        }
        else {
            FD_app_cir_fourth_standard (rptr, b, dimx, dimy, dimz);
        }
        my_free(rptr);
        return;
    }
    if(order == APP_CI_SIXTH) {
        T.CPP_trade_imagesx (a, rptr, dimx, dimy, dimz, 2, FULL_FD);
        if(numgrid == P0_BASIS) {
            FD_app_cir_sixth_global (rptr, b);
        }
        else {
            FD_app_cir_sixth_standard (rptr, b, dimx, dimy, dimz);
        }
        my_free(rptr);
        return;
    }

    rmg_error_handler("APP_CIR order not programmed yet in CPP_app_cir_driver.\n");

}

extern "C" void app_cir_driver (rmg_double_t * a, rmg_double_t * b, int dimx, int dimy, int dimz, int order)
{
    CPP_app_cir_driver<double>(a, b, dimx, dimy, dimz, order);
}
extern "C" void app_cir_driver_f (rmg_float_t * a, rmg_float_t * b, int dimx, int dimy, int dimz, int order)
{
    CPP_app_cir_driver<float>(a, b, dimx, dimy, dimz, order);
}
