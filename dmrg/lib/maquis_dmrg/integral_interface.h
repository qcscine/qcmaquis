/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
*               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
*               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
*               2014-2014    Sebastian Keller <sebkelle@phys.ethz.ch>
*               2019         Leon Freitag <lefreita@ethz.ch>
*               2020- by Robin Feldmann <robinfe@phys.chem.ethz.ch>
*               2021- by Alberto Baiardi <abaiardi@phys.chem.ethz.ch>
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

#include "dmrg/utils/align.h"
#include <unordered_map>
#include <boost/serialization/serialization.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/functional/hash.hpp>
#include "integral_helper.h"

namespace chem {

/**
 * @brief Class representing the Hamiltonian as a map index_tuple --> factor
 * 
 * Note that permutations are handled internally.
 * Indexing is as in the FCIDUMP file, i.e:
 * 
 * 1. Orbital indices start from 1 and 2e integrals use all four indices
 * 2. 1e integrals use the first two indices and 0,0 as the 3rd and the 4th index
 * 3. Nuclear repulsion energy uses an index 0,0,0,0 
 * 
 * @tparam V Type associated with the scalar factors of the Hamiltonian
 * @tparam HamiltonianType Enum class representing the Hamiltonian type.
 * @tparam Transcorrelation Enum class representing whether the Hamiltonian was transcorrelated.
 */
template<class V, Hamiltonian HamiltonianType=Hamiltonian::Electronic, HamiltonianTransformation Transcorrelation=HamiltonianTransformation::Conventional>
class integral_map
{
public:
    using map_t = std::unordered_map<index_type<HamiltonianType, Transcorrelation>, V, integral_hash<HamiltonianType, Transcorrelation>>;
    using size_type = typename map_t::size_type;
    // Type which returns std::abs(V), for the integral cutoff
    // Not very clean but std::conditional seems not to work here
    using value_type = typename std::complex<V>::value_type;
    using iterator = typename map_t::iterator;
    using const_iterator = typename map_t::const_iterator;

    /** @brief Default constructor */
    integral_map() = default;

    /**
     * @brief Copy constructor
     * Explicit copy using this->operator[]() to avoid potential doubling due to symmetry permutation
     * @param map object that is copied from
     * @param cutoff Threshold for accepting integrals
     */
    explicit integral_map(const map_t & map, value_type cutoff=0.0) : cutoff_(cutoff)
    {
        for (auto&& it: map)
            (*this)[it->first] = it->second;
    }

    /** @brief Move constructor */
    explicit integral_map(map_t && map, value_type cutoff=0.0) : map_(map), cutoff_(cutoff) {};

    /** @brief Initializer list construction */
    integral_map(std::initializer_list<typename map_t::value_type> l, value_type cutoff=0.0) 
        : integral_map(map_t(l), cutoff) {};

    // Iterator classes
    iterator begin() { return map_.begin(); };
    const_iterator begin() const { return map_.begin(); };
    iterator end() { return map_.end(); };
    const_iterator end() const { return map_.end(); };

    // For complex integrals, use relativistic permutation. Otherwise, use nonrelativistic permutation.
    // Maybe these two properties should be decoupled in the future.
    V& operator[](const index_type<HamiltonianType>& key) { return map_[maquis::detail::AlignTraitTypeClass<V>::align(key)]; }
    const V& operator[](const index_type<HamiltonianType>& key) const {  return map_[maquis::detail::AlignTraitTypeClass<V>::align(key)]; }
    V& at(const index_type<HamiltonianType>& key) { return map_.at(maquis::detail::AlignTraitTypeClass<V>::align(key)); }
    const V& at(const index_type<HamiltonianType>& key) const { return map_.at(maquis::detail::AlignTraitTypeClass<V>::align(key)); }

    /** @brief Size getter */
    size_type size() const { return map_.size(); }

private:
    friend class boost::serialization::access;
    // Map storing the data
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
template <class V, Hamiltonian HamiltonianType=Hamiltonian::Electronic>
std::string serialize(const integral_map<V, HamiltonianType>& ints)
{
    std::stringstream ss;
    boost::archive::text_oarchive oa{ss};
    oa << ints;
    return ss.str();
}

} // namespace chem

#endif