/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Alexandr Kosenkov <alex.kosenkov@gmail.com>
 *                         by Michele Dolfi <dolfim@phys.ethz.ch>
 * 
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 * 
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#ifndef PARALLEL_FOR_HPP
#define PARALLEL_FOR_HPP

#ifdef USE_AMBIENT
    #define select_scope(...) ambient::actor ctxt(__VA_ARGS__)
    #ifndef AMBIENT_SERIAL_FOR
    #define parallel_for(control_variable, loop_range, ...) ambient::threaded_for_each(loop_range.begin(), loop_range.end(), [&](control_variable) __VA_ARGS__);
    #endif
    #define threaded_for(...) for(__VA_ARGS__)
    #define omp_critical
#elif defined(MAQUIS_OPENMP)
    #define select_scope(...) 
    #define parallel_pragma(a) _Pragma( #a )
    #define threaded_for(...) parallel_pragma(omp parallel for schedule(dynamic, 1)) for(__VA_ARGS__)
    #define omp_critical parallel_pragma(omp critical)
#else
    #define select_scope(...) 
    #define threaded_for(...) for(__VA_ARGS__)
    #define omp_critical
#endif

#define omp_for(control_variable, loop_range, ...)  {\
    int dist = loop_range.end() - loop_range.begin(); \
    threaded_for(int loop_control = 0; loop_control < dist; loop_control++){ \
        control_variable = loop_range.begin() + loop_control; \
        __VA_ARGS__ \
    } \
}

#if !defined(USE_AMBIENT) || defined(AMBIENT_SERIAL_FOR)
    #define parallel_for omp_for
#endif

template<typename T>
struct dynamic_range {
    dynamic_range(T first, T second) : first(first), second(second) { }
    T begin() const { return first; }
    T end() const { return second; }
    size_t size() const { return end()-begin(); }
    const T first;
    const T second;
};

template<typename T>
dynamic_range<T> range(T first, T second){
    return dynamic_range<T>(first, second);
}

#endif
