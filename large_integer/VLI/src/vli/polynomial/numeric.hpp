#ifndef VLI_POLYNOMIAL_NUMERIC_HPP
#define VLI_POLYNOMIAL_NUMERIC_HPP

namespace vli {
/*
  The following templates are used in the polynomial (and elsewhere)
  and can be specialized or overloaded for types to provide an optimized
  implementation.
*/

template <typename T>
bool is_zero(T t) {
    return t == 0;
}

template <typename T>
void negate_inplace(T& t) {
    t = -t;
}

template <typename T, typename T2, typename T3>
void multiply_add(T& t, T2 const& t2, T3 const& t3) {
    t += t2 * t3;
}

} // end namespace vli

#endif //VLI_POLYNOMIAL_NUMERIC_HPP
