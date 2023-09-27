/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef PARALLEL_LOOPS_HPP
#define PARALLEL_LOOPS_HPP

#if defined(MAQUIS_OPENMP)
    #define parallel_pragma(a) _Pragma( #a )
    #define threaded_for(...) parallel_pragma(omp parallel for schedule(dynamic, 1)) for(__VA_ARGS__)
    #define parallel_critical parallel_pragma(omp critical)
#else
    #define threaded_for(...) for(__VA_ARGS__)
    #define parallel_critical
#endif

#define omp_for(control_variable, loop_range, ...)  {\
    int dist = loop_range.end() - loop_range.begin(); \
    threaded_for(int loop_control = 0; loop_control < dist; loop_control++){ \
        control_variable = loop_range.begin() + loop_control; \
        __VA_ARGS__ \
    } \
}

#if defined(USE_AMBIENT) && !defined(AMBIENT_SERIAL_COLLECTION)
    #include "ambient/utils/threaded_for_each.hpp"
    #define parallel_for(control_variable, loop_range, ...) ambient::threaded_for_each(loop_range.begin(), loop_range.end(), [&](control_variable) __VA_ARGS__);
#else
    #define parallel_for omp_for
#endif

#endif
