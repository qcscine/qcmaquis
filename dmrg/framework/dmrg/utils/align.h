/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
*               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
*               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
*               2014-2014    Sebastian Keller <sebkelle@phys.ethz.ch>
*               2020         Leon Freitag <lefreita@ethz.ch>
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