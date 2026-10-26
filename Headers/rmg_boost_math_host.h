#ifndef RMG_BOOST_MATH_HOST_H
#define RMG_BOOST_MATH_HOST_H
// Force-included in SYCL builds (see CMakeLists.txt). Boost.Math >= 1.86 tags its functions
// SYCL_EXTERNAL when SYCL_LANGUAGE_VERSION is defined, and icpx then fails on the recursive
// and function-pointer code in bessel/roots/minima. RMG uses Boost.Math only on the host, so
// load its config once with the macro hidden; the include guard keeps that choice for every
// later Boost.Math include in the translation unit.
#if defined(SYCL_LANGUAGE_VERSION) && __has_include(<boost/math/tools/config.hpp>)
#pragma push_macro("SYCL_LANGUAGE_VERSION")
#undef SYCL_LANGUAGE_VERSION
#include <boost/math/tools/config.hpp>
#pragma pop_macro("SYCL_LANGUAGE_VERSION")
#endif
#endif
