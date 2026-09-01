#include <complex>
#include "common_prototypes.h"
template <> double MyConj(double val)
{
    return val;
}
template <> std::complex<double> MyConj(std::complex<double> val)
{
    return std::conj(val);
}
