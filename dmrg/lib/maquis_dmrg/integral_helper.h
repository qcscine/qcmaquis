/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2022 Institute for Theoretical Physics, ETH Zurich
*               2022- by Alberto Baiardi <abaiardi@phys.chem.ethz.ch>
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

#ifndef INTEGRAL_HELPER_H
#define INTEGRAL_HELPER_H

#include <complex>

namespace chem {

/** @brief Enum class distinguishing the possible types of Hamiltonians */
enum class Hamiltonian {Electronic, ElectronicTranscorrelated,
                        VibrationalCanonical, VibrationalNMode,
                        PreBO, Vibronic, Excitonic};

/** @brief Constexpr function returning whether an Hamiltonian is Hermitian */
constexpr bool isModelHermitian(const Hamiltonian& type) {
    bool ret = true;
    switch (type) {
        case Hamiltonian::ElectronicTranscorrelated:
            ret = false;
            break;
        default:
            ret = true;
            break;
    }
    return ret;
}

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
 *   equal to the max. order of the Taylor expansion of the PES. Here we include
 *   up to sixth-order force constants.
 */
constexpr int getIndexDim(const Hamiltonian& type) {
    int indexDim=0;
    switch (type) {
        case Hamiltonian::Electronic:
            indexDim = 4;
            break;
        case Hamiltonian::ElectronicTranscorrelated:
            indexDim = 6;
            break;
        case Hamiltonian::VibrationalCanonical:
            indexDim = 6;
            break;
        // Note that we support so-far only up to 3-body terms
        case Hamiltonian::VibrationalNMode:
            indexDim = 12;
            break;
        // We support up to 2-mode coupling for the vibronic case.
        // The index is, however, 4 because we also include the electronic
        // state index (same fore the excitonic case)
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

} // namespace chem

#endif // INTEGRAL_HELPER_H