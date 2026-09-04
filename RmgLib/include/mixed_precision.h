/*
 *
 * Copyright (c) 1995, Emil Briggs
 * Copyright (C) 1998  Emil Briggs, Charles Brabec, Mark Wensell, 
 *                     Dan Sullivan, Chris Rapcewicz, Jerzy Bernholc
 * Copyright (C) 2001  Emil Briggs, Wenchang Lu,
 *                     Marco Buongiorno Nardelli,Charles Brabec, 
 *                     Mark Wensell,Dan Sullivan, Chris Rapcewicz,
 *                     Jerzy Bernholc
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.

 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
*/


/* Helper code used to clean up code that accumulates lower precision mixed real and complex
   data.
*/

#ifndef RMG_mixed_precision_H
#define RMG_mixed_precision_H 1


#include <complex>
#include <type_traits>

template <typename T>
struct is_complex : std::false_type {};

template <typename T>
struct is_complex<std::complex<T>> : std::true_type {};

template <typename T>
inline constexpr bool is_complex_v = is_complex<T>::value;

// 
template <typename T>
using accumulator_t = std::conditional_t<
    is_complex_v<T>,
    std::complex<double>,
    double
>;

// std::complex<float> * std::complex<double>
inline std::complex<double> operator*(const std::complex<float>& lhs, const std::complex<double>& rhs) {
    return std::complex<double>(lhs) * rhs;
}

// std::complex<double> * std::complex<float>
inline std::complex<double> operator*(const std::complex<double>& lhs, const std::complex<float>& rhs) {
    return lhs * std::complex<double>(rhs);
}

// std::complex<float> + std::complex<double>
inline std::complex<double> operator+(const std::complex<float>& lhs, const std::complex<double>& rhs) {
    return std::complex<double>(lhs) + rhs;
}

// std::complex<double> + std::complex<float>
inline std::complex<double> operator+(const std::complex<double>& lhs, const std::complex<float>& rhs) {
    return lhs + std::complex<double>(rhs);
}

#endif
