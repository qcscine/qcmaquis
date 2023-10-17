/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */
#ifndef ALIGN_H
#define ALIGN_H

#include <array>

namespace maquis{
    namespace detail {

        typedef std::array<int, 4> index_type;

        // Permutes integral indices to yield the canonical form with i>=j, k>=l
        template<bool P=false>
        inline index_type align(const index_type & idx)
        {
            int i = idx[0], j = idx[1], k = idx[2], l = idx[3];
            if (i<j) std::swap(i,j);
            if (k<l) std::swap(k,l);
            if (i<k) { std::swap(i,k); std::swap(j,l); }
            if (i==k && j<l) { std::swap(j,l); }
            return index_type{i,j,k,l};
        }

        // Disables permutation (e.g. for relativistic calculations)
        template<>
        inline index_type align<true>(const index_type & idx)
        {
            return idx;
        }
    }
}
#endif