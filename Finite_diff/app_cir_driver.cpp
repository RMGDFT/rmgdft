#include "TradeImages.h"
#include "Lattice.h"
#include "FiniteDiff.h"
#include "rmg_error.h"
#include "RmgTimer.h"

template void CPP_app_cir_driver<float>(float *, float *, int, int, int, int);
template void CPP_app_cir_driver<double>(double *, double *, int, int, int, int);

template <typename RmgType>
void CPP_app_cir_driver (RmgType * a, RmgType * b, int dimx, int dimy, int dimz, int order)
{

    RmgTimer RT("App_cir");
    int sbasis;
    TradeImages T;
    Lattice L;
    FiniteDiff FD(&L);;
    sbasis = (dimx + 4) * (dimy + 4) * (dimz + 4);
    RmgType *rptr = new RmgType[sbasis + 64];

    if(order == APP_CI_FOURTH) {
        RmgTimer *RT1 = new RmgTimer("App_cir: trade images");
        T.trade_imagesx (a, rptr, dimx, dimy, dimz, 1, FULL_TRADE);
        delete(RT1);
        FD.app_cir_fourth (rptr, b, dimx, dimy, dimz);
    }
    else if(order == APP_CI_SIXTH) {
        RmgTimer *RT1 = new RmgTimer("App_cir: trade images");
        T.trade_imagesx (a, rptr, dimx, dimy, dimz, 2, FULL_TRADE);
        delete(RT1);
        FD.app_cir_sixth (rptr, b, dimx, dimy, dimz);
    }
    else {
        rmg_error_handler (__FILE__, __LINE__, "APP_CIR order not programmed yet in CPP_app_cir_driver.\n");
    }

    delete [] rptr;
    return;

}

