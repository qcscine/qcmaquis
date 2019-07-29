/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
*               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
*               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
*               2014-2014    Sebastian Keller <sebkelle@phys.ethz.ch>
*               2019         Leon Freitag <lefreita@ethz.ch>
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
#ifndef MEASUREMENT_MAP_H
#define MEASUREMENT_MAP_H

#include <unordered_map>
#include <array>
#include <boost/functional/hash.hpp>

namespace measurements {

    template <std::size_t N>
    using index_array = std::array<int, N>;

    // Hash function for the unordered map
    template <std::size_t N>
    struct index_hash
    {
        public:
            std::size_t operator()(const index_array<N>& id) const
            {
                return boost::hash_range(id.begin(), id.end());
            }
    };

    // A map for measurement results
    template <class V, std::size_t N>
    using measurement_map = std::unordered_map<index_array<N>, V, index_hash<N> >;

    // Partial specializations for 1- and 2-RDM
    template <class V>
    using onerdm_map = measurement_map<V, 2>;

    template <class V>
    using twordm_map = measurement_map<V, 4>;
}

#endif