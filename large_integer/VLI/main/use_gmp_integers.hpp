#ifndef USE_GMP_INTEGERS_HPP
#define USE_GMP_INTEGERS_HPP

#include <gmp.h>
#include <gmpxx.h>


namespace hp2c
{

// A large integer of 128-256 bits (fixed size)
// We use the gmp integer class
typedef mpz_class large_int;


// std out for gmp
std::ostream& operator << (std::ostream& o, mpz_class const &mpz);
//
std::ostream& operator << (std::ostream& o, mpz_class const &mpz)
{
    o<<mpz.get_str();
    return o;
}

}

#endif //USE_GMP_INTEGERS_HPP
