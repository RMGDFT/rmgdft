
#include "const.h"
#include "common_prototypes.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include "transition.h"
#include "rmg_error.h"
#include "FiniteDiff.h"


template void AppCilrDriver<double>(TradeImages *, double *, double *, double *, double *, int, int, int, double, double, double, int);
template void AppCilrDriver<std::complex<double> >(TradeImages *, std::complex<double> *, std::complex<double> *, std::complex<double> *, double *, int, int, int, double, double, double, int);

template  <typename OrbitalType> void AppCilrDriver (TradeImages *T, OrbitalType * psi, OrbitalType * a_psi, OrbitalType *b_psi, double *vtot, 
    int dimx, int dimy, int dimz, double hx, double hy, double hz, int order)
{

    OrbitalType *rptr = new OrbitalType[(dimx + 4) * (dimy + 4) * (dimz + 4)];

    if(order == APP_CI_FOURTH) {

        T->trade_imagesx (psi, rptr, dimx, dimy, dimz, 1, FULL_TRADE);
        AppCilrFourth(rptr, a_psi, b_psi, vtot, dimx, dimy, dimz, hx, hy, hz);

    }
    else if(order == APP_CI_SIXTH) {

        T->trade_imagesx (psi, rptr, dimx, dimy, dimz, 2, FULL_TRADE);
        AppCilrSixth(rptr, a_psi, b_psi, vtot, dimx, dimy, dimz, hx, hy, hz);

    }
    else {
        rmg_error_handler(__FILE__, __LINE__, "Finite difference order not programmed yet in AppCilrDriver.\n");

    }
    
    delete [] rptr;
    

}

