#include "TradeImages.h"
#include "Lattice.h"
#include "FiniteDiff.h"
#include "rmg_error.h"
#include "RmgTimer.h"

using namespace std;

template double CPP_app_cil_driver<float>(float *, float *, int, int, int, double, double, double, int);
template double CPP_app_cil_driver<double>(double *, double *, int, int, int, double, double, double, int);

template <typename RmgType>
double CPP_app_cil_driver (RmgType * a, RmgType * b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz, int order)
{

    RmgTimer RT("App_cil");
    int sbasis;
    double cc;
    TradeImages T;
    Lattice L;
    FiniteDiff FD(&L);
    sbasis = (dimx + 4) * (dimy + 4) * (dimz + 4);
    RmgType *rptr = new RmgType[sbasis + 64];

    if(order == APP_CI_FOURTH) {
        RmgTimer *RT1 = new RmgTimer("App_cil: trade images");
        T.trade_imagesx (a, rptr, dimx, dimy, dimz, 1, FULL_TRADE);
        delete(RT1);
        cc = FD.app_cil_fourth (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);
    }
    else if(order == APP_CI_SIXTH) {
        RmgTimer *RT1 = new RmgTimer("App_cil: trade images");
        T.trade_imagesx (a, rptr, dimx, dimy, dimz, 2, FULL_TRADE);
        delete(RT1);
        cc = FD.app_cil_sixth (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);
    }
    else {
        rmg_error_handler (__FILE__, __LINE__, "APP_CIL order not programmed yet in app_cil_driver.\n");
    }

    delete [] rptr;
    return cc;

}

