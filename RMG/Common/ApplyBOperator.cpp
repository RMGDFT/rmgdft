#include <complex>
#include "const.h"
#include "TradeImages.h"
#include "Lattice.h"
#include "FiniteDiff.h"
#include "rmg_error.h"
#include "transition.h"

template void ApplyBOperator<float>(Lattice *, TradeImages *, float *, float *, int, int, int, int);
template void ApplyBOperator<double>(Lattice *, TradeImages *, double *, double *, int, int, int, int);
template void ApplyBOperator<std::complex<float> >(Lattice *, TradeImages *, std::complex<float> *, std::complex<float> *, int, int, int, int);
template void ApplyBOperator<std::complex<double> >(Lattice *, TradeImages *, std::complex<double> *, std::complex<double> *, int, int, int, int);

template <typename RmgType>
void ApplyBOperator (Lattice *L, TradeImages *T, RmgType * a, RmgType * b, int dimx, int dimy, int dimz, int order)
{

    if(ct.discretization == MEHRSTELLEN_DISCRETIZATION) {

        CPP_app_cir_driver (L, T, a, b, dimx, dimy, dimz, order);

    }
    else if(ct.discretization == CENTRAL_DISCRETIZATION) {

        for(int ix = 0;ix < dimx*dimy*dimz;ix++) b[ix] = a[ix];

    }
    else {

        rmg_error_handler(__FILE__, __LINE__, "Unknown discretization method.");

    }

}
