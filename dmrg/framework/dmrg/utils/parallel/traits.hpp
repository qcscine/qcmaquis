/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef PARALLEL_TRAITS_HPP
#define PARALLEL_TRAITS_HPP

#include <vector>
#include <cassert>

namespace parallel {

    #ifdef USE_AMBIENT
    struct traits {
        typedef int resource_iterator_nop;
        typedef typename ambient::scope::const_iterator resource_iterator;
        static resource_iterator balance(int k, int max){ return ambient::scope::balance(k, max); }
        static resource_iterator permute(int k, const std::vector<int>& s, size_t gran){ if(k >= s.size()) return to_iterator(k % size()); return ambient::scope::permute(k, s, gran); }
        static size_t size(){ return ambient::scope::size(); }
        static size_t distance(const resource_iterator& it){ return (it - ambient::scope::begin()); }
        template<typename T>
        static resource_iterator to_iterator(T offset){ assert(offset < size()); return (ambient::scope::begin() + offset); }
        template<typename T>
        static int placement(const T& obj){ return ambient::get_owner(obj); }
    };
    #else
    struct traits {
        typedef int resource_iterator_nop;
        typedef resource_iterator_nop resource_iterator;
        static resource_iterator balance(int k, int max){ return resource_iterator(); }
        static resource_iterator permute(int k, const std::vector<int>& s, size_t){ return resource_iterator(); }
        static size_t size(){ return 1; }
        template<typename T>
        static resource_iterator to_iterator(T offset){ return resource_iterator(); }
        template<typename T>
        static int placement(const T& obj){ return 0; }
    };
    #endif

}

#endif
