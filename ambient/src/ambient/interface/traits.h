#ifndef AMBIENT_TRAITS_H
#define AMBIENT_TRAITS_H

#include <complex>
namespace ambient {

    struct traits {
    //  typedef std::complex<double> value_type; 
        typedef double value_type; 
        enum {value=128};
    };

}
#endif
