#include <complex>
#include "const.h"
#include "TradeImages.h"
#include "Lattice.h"
#include "FiniteDiff.h"
#include "rmg_error.h"
#include "transition.h"

template double ApplyAOperator<float>(Lattice *, TradeImages *, float *, float *, int, int, int, double, double, double, int);
template double ApplyAOperator<double>(Lattice *, TradeImages *, double *, double *, int, int, int, double, double, double, int);
template double ApplyAOperator<std::complex<float> >(Lattice *, TradeImages *, std::complex<float> *, std::complex<float> *, int, int, int, double, double, double, int);
template double ApplyAOperator<std::complex<double> >(Lattice *, TradeImages *, std::complex<double> *, std::complex<double> *, int, int, int, double, double, double, int);

template <typename DataType>
double ApplyAOperator (Lattice *L, TradeImages *T, DataType *a, DataType *b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz, int order)
{

    if(ct.discretization == MEHRSTELLEN_DISCRETIZATION) {

        return CPP_app_cil_driver (L, T, a, b, dimx, dimy, dimz, gridhx, gridhy, gridhz, order);

    }
    else if(ct.discretization == CENTRAL_DISCRETIZATION) {

        double cc = CPP_app_del2_driver (L, T, a, b, dimx, dimy, dimz, gridhx, gridhy, gridhz, order);
        return cc;

    }
    
    rmg_error_handler(__FILE__, __LINE__, "Unknown discretization method."); 
    return false;
}

