#ifndef __MAQUIS_TYPES_KERNELS_UTILS_HPP__
#define __MAQUIS_TYPES_KERNELS_UTILS_HPP__

#include <limits>
#include "utils/timings.h"

//#define AMBIENT_COMPUTATIONAL_TIMINGS
//#define AMBIENT_CHECK_BOUNDARIES

#ifdef AMBIENT_CHECK_BOUNDARIES
#include <execinfo.h>
#endif

extern "C" {
    double ddot_(const int*, const double*, const int*, const double*, const int*);
}

namespace ambient { namespace numeric { namespace kernels {

    using ambient::numeric::matrix_impl;

    #include "ambient/utils/numeric.h" // BLAS/LAPACK prototypes
    #include "ambient/utils/ceil.h"
   
    #ifdef AMBIENT_COMPUTATIONAL_TIMINGS
        #define __A_TIME_C(name) static __a_timer time(name); time.begin();
        #define __A_TIME_C_STOP time.end();
    #else
        #define __A_TIME_C(name) 
        #define __A_TIME_C_STOP 
    #endif


} } }

#endif
