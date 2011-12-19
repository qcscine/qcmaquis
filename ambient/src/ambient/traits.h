#ifndef AMBIENT_TRAITS_H
#define AMBIENT_TRAITS_H

#include <complex>
namespace ambient{
    struct traits
    {
        //typedef std::complex<double> type; 
        typedef double type; 
        enum {value=128}; // C - could be given by the preprocessor -D ? 
    };
}
#endif
