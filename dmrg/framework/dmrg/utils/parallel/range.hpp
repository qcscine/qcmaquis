/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef PARALLEL_RANGE_HPP
#define PARALLEL_RANGE_HPP

#include <cstddef>

namespace parallel {

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

}

#endif
