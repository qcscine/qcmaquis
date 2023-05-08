/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

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

namespace chem {

/** @brief Enum class distinguishing the possible types of Hamiltonians */
enum class Hamiltonian {Electronic, VibrationalCanonical, VibrationalNMode, PreBO, Vibronic, Excitonic};

/**
 * @brief Constexpr function returning the index of the Hamiltonian map
 *
 * This function is required because different models identify a
 * Hamiltonian term with different formats.
 *
 * - for the electronic Hamiltonian, each term is identified by 4 indices,
 *   because we have a two-body potential, and one index is sufficient
 *   to identify each SQ operator.
 * - for the PreBO Hamiltonian, we still have a two-body Hamiltonian, but
 *   each SQ operator is identified by 2 indices, the particle type and
 *   the orbital number.
 * - the n-mode Hamiltonian has, in principle, an arbitrarily high coupling
 *   degree. We include, here, up to three-body coupling terms, and therefore
 *   we have up to 12 indices (note that in the n-mode Hamiltonian each SQ
 *   operator is identified by 2 indices, as for the PreBO one)
 * - the canonical quantization-based vibration Hamiltonian has a number of indices
 *   equal to the max. order of the Taylor expansion of the PES.
 *   The maximum order of the force constants can be set at compilation.
 */
constexpr int getIndexDim(const Hamiltonian& type) {
    int indexDim=0;
    switch (type) {
        case Hamiltonian::Electronic:
            indexDim = 4;
            break;
        case Hamiltonian::VibrationalCanonical:
            indexDim = 6;
            break;
        // Note that we support so-far only up to 3-body terms
        case Hamiltonian::VibrationalNMode:
            indexDim = 12;
            break;
        // We support up to 2-mode coupling for the vibronic case.
        // The index is, however, 4 because we also include the electronic state index
        case Hamiltonian::Vibronic:
            indexDim = 4;
            break;
        case Hamiltonian::Excitonic:
            indexDim = 2;
            break;
        case Hamiltonian::PreBO:
            indexDim = 8;
            break;
    }
    return indexDim;
}


/** @brief Class associated with the index identifying a single SQ operator */
template <Hamiltonian HamiltonianType=Hamiltonian::Electronic, int N = getIndexDim(HamiltonianType)>
using index_type = std::array<int, N>;

/** @brief Class associated with a single entry of the Hamiltonian */
template <class V, Hamiltonian HamiltonianType=Hamiltonian::Electronic>
using integral_tuple = std::pair<index_type<HamiltonianType>, V>;

/** @brief Class associated with the overall Hamiltonian */
template <class V, Hamiltonian HamiltonianType=Hamiltonian::Electronic>
using integrals = std::vector<integral_tuple<V, HamiltonianType> >; // TODO: use a map later

// Structs needed for distinguishing whether we have a complex type or not
// required for proper integral permutation rules in the integral_map.
template<typename T>
struct is_complex_t : public std::false_type {};
template<typename T>
struct is_complex_t<std::complex<T> > : public std::true_type {};


/** @brief Hasing function for a single Hamiltonian entry */
template <Hamiltonian HamiltonianType=Hamiltonian::Electronic>
struct integral_hash
{
    public:
        std::size_t operator()(const index_type<HamiltonianType>& id) const
        {
            return boost::hash_range(id.begin(), id.end());
        }
};

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
 * @tparam HamiltonianType Enum class representing the
 */
template <class V, Hamiltonian HamiltonianType=Hamiltonian::Electronic>
class integral_map
{
public:
    typedef std::unordered_map<index_type<HamiltonianType>, V, integral_hash<HamiltonianType>> map_t;
    typedef typename map_t::size_type size_type;
    // Type which returns std::abs(V), for the integral cutoff
    // Not very clean but std::conditional seems not to work here
    typedef typename std::complex<V>::value_type value_type;
    typedef typename map_t::iterator iterator;
    typedef typename map_t::const_iterator const_iterator;

    /** @brief Default constructor */
    integral_map() = default;

    /**
     * @brief Copy constructor
     *
     * Explicit copy using this->operator[]() to avoid potential doubling due to symmetry permutation
     *
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

    // For complex integrals, use relativistic permutation. Otherwise, use nonrelativistic permutation
    // Maybe these two properties should be decoupled in the future
    V& operator[](const index_type<HamiltonianType> & key) { return map_[maquis::detail::align<is_complex_t<V>::value>(key)]; };
    const V& operator[](const index_type<HamiltonianType> & key) const { return map_[maquis::detail::align<is_complex_t<V>::value>(key)]; };
    V& at(const index_type<HamiltonianType> & key) { return map_.at(maquis::detail::align<is_complex_t<V>::value>(key)); };
    const V& at(const index_type<HamiltonianType> & key) const { return map_.at(maquis::detail::align<is_complex_t<V>::value>(key)); };

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
