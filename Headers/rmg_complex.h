#ifndef RMG_rmg_complex_H
#define RMG_rmg_complex_H 1

/* Experimental header that enables mixing of basic arithmetic operations
   for std::complex<float> with double and std::complex<double>. Use with
   care until fully verified.
*/


template <class T, class U>
auto operator*(T const &a, U const &b) -> typename std::common_type<T, U>::type
{
    using rtype = typename std::common_type<T, U>::type;
    return rtype(a) * rtype(b);
}

template <class T, class U>
auto operator/(T const &a, U const &b) -> typename std::common_type<T, U>::type
{
    using rtype = typename std::common_type<T, U>::type;
    return rtype(a) / rtype(b);
}

template <class T, class U>
auto operator+(T const &a, U const &b) -> typename std::common_type<T, U>::type
{
    using rtype = typename std::common_type<T, U>::type;
    return rtype(a) + rtype(b);
}

template <class T, class U>
auto operator-(T const &a, U const &b) -> typename std::common_type<T, U>::type
{
    using rtype = typename std::common_type<T, U>::type;
    return rtype(a) - rtype(b);
}

#endif
