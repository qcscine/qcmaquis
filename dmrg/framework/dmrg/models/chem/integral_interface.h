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
#ifndef INTEGRAL_INTERFACE_H
#define INTEGRAL_INTERFACE_H

#include <unordered_map>
#include "dmrg/models/chem/util.h"
#include <boost/serialization/serialization.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <boost/serialization/utility.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/functional/hash.hpp>

namespace chem {

    // // Leon: TODO:
    // // distinguish between electronic or vibrational
    // // This way we can merge parse_integrals_vib.h here

    // // The only thing that changes in electronic vs vibrational integrals is the
    // // number of indices in the FCIDUMP file: 4 in the electronic and 6 in the vibrational
    // const int Electronic = 4, Vibrational = 6;
    //
    // template <int N = Electronic>
    // using index_type = std::array<int, N>;

    // Classes to represent integrals

    typedef std::array<int, 4> index_type;

    template <class V>
    using integral_tuple = std::pair<index_type, V>;

    template <class V>
    using integrals = std::vector<integral_tuple<V> >; // TODO: use a map later

    // structs needed for distinguishing whether we have a complex type or not
    // required for proper integral permutation rules in the integral_map
    template<typename T>
    struct is_complex_t : public std::false_type {};
    template<typename T>
    struct is_complex_t<std::complex<T> > : public std::true_type {};

    struct integral_hash
    {
        public:
            std::size_t operator()(const index_type& id) const
            {
                return boost::hash_range(id.begin(), id.end());
            }
    };

    // Map that maps the four indices to an integral value
    // and handles index permutations internally.
    // Indexing is as in the FCIDUMP file, i.e.
    // Orbital indices start from 1 and 2e integrals use all four indices
    // 1e integrals use the first two indices and 0,0 as the 3rd and the 4th index
    // Nuclear repulsion energy uses an index 0,0,0,0

    template <class V>
    class integral_map
    {
        public:
            typedef std::unordered_map<index_type, V, integral_hash> map_t;
            typedef typename map_t::size_type size_type;

            // For complex integrals, use relativistic permutation. Otherwise, use nonrelativistic permutation
            // Maybe these two properties should be decoupled in the future
            typedef typename std::conditional<is_complex_t<V>::value, U1DG, TrivialGroup>::type relativistic_t;

            // Type which returns std::abs(V), for the integral cutoff
            // Not very clean but std::conditional seems not to work here
            typedef typename std::complex<V>::value_type value_type;

            integral_map() = default;

            // Explicit copy using this->operator[]() to avoid potential doubling due to symmetry permutation
            // Not implementing it in the move constructor yet
            integral_map(const map_t & map, value_type cutoff=0.0) : cutoff_(cutoff)
            {
                for (auto&& it: map)
                    (*this)[it->first] = it->second;
            }

            integral_map(map_t && map, value_type cutoff=0.0) : map_(map), cutoff_(cutoff) {};

            // allow initializer lists for construction
            integral_map(std::initializer_list<typename map_t::value_type> l, value_type cutoff=0.0) : integral_map(map_t(l), cutoff) {};

            typedef typename map_t::iterator iterator;
            typedef typename map_t::const_iterator const_iterator;

            iterator begin() { return map_.begin(); };
            const_iterator begin() const { return map_.begin(); };
            iterator end() { return map_.end(); };
            const_iterator end() const { return map_.end(); };

            V& operator[](const index_type & key) { return map_[detail::align<relativistic_t>(key)]; };
            const V& operator[](const index_type & key) const { return map_[detail::align<relativistic_t>(key)]; };
            V& at(const index_type & key) { return map_.at(detail::align<relativistic_t>(key)); };
            const V& at(const index_type & key) const { return map_.at(detail::align<relativistic_t>(key)); };

            size_type size() const { return map_.size(); }
        private:
            friend class boost::serialization::access;

            map_t map_;

            // Integral cutoff
            value_type cutoff_;

            template <typename Archive>
            friend void serialize(Archive& ar, integral_map &i, const unsigned int version)
            {
                ar & i.map_;
            }
    };

    // Serialize the integral into a string

    template <class V>
    std::string serialize(const integral_map<V>& ints)
    {
        std::stringstream ss;
        boost::archive::text_oarchive oa{ss};
        oa << ints;

        return ss.str();
    }

}

#endif