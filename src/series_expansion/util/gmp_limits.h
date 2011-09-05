#ifndef GMP_LIMITS_H
#define GMP_LIMITS_H

#include <limits>



namespace std
{
  template<>
    struct numeric_limits<mpz_class>
    {
      static const bool is_specialized = true;
      /*
         // TODO
      static int min() throw()
      { return -__INT_MAX__ - 1; }
      static int max() throw()
      { return __INT_MAX__; }

      static const int digits = __glibcxx_digits (int);
      static const int digits10 = __glibcxx_digits10 (int);
      */
      static const bool is_signed = true;
      static const bool is_integer = true;
      static const bool is_exact = true;
      //static const int radix = 2;
      static int epsilon() throw()
      { return 0; }
      static int round_error() throw()
      { return 0; }

      static const int min_exponent = 0;
      static const int min_exponent10 = 0;
      static const int max_exponent = 0;
      static const int max_exponent10 = 0;

      static const bool has_infinity = false;
      static const bool has_quiet_NaN = false;
      static const bool has_signaling_NaN = false;
      static const float_denorm_style has_denorm = denorm_absent;
      static const bool has_denorm_loss = false;

      static int infinity() throw()
      { return static_cast<int>(0); }
      static int quiet_NaN() throw()
      { return static_cast<int>(0); }
      static int signaling_NaN() throw()
      { return static_cast<int>(0); }
      static int denorm_min() throw()
      { return static_cast<int>(0); }

      static const bool is_iec559 = false;
      /*
         //TODO
      static const bool is_bounded = true;
      static const bool is_modulo = true;

      static const bool traps = __glibcxx_integral_traps;
      static const bool tinyness_before = false;
      static const float_round_style round_style = round_toward_zero;
      */
    };
}
#endif //GMP_LIMITS_H
