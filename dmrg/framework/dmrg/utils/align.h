/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
*               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
*               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
*               2014-2014    Sebastian Keller <sebkelle@phys.ethz.ch>
*               2020-        Leon Freitag <lefreita@ethz.ch>
*               2022-        Alberto Baiardi <abaiardi@ethz.ch>
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
#include "dmrg/block_matrix/symmetry.h"

namespace maquis {
namespace detail {

template<class SymmGroup, int NumberOfElements>
auto align(std::initializer_list<int> rhs);

/** @brief Trait class with the index info depending on the Hamiltonian */
template<class SymmGroup>
class AlignTraitClass {
public:
    static constexpr bool doRealign=true;
};

template<>
class AlignTraitClass<TrivialGroup> {
public:
    static constexpr bool doRealign=false;
};

template<>
class AlignTraitClass<U1DG> {
public:
    static constexpr bool doRealign=false;
};

/** @brief Align function for the 4-element array */
inline auto alignArray(const std::array<int, 4>& idx, bool isHermitian) 
{
    int i = idx[0];
    int j = idx[1];
    int k = idx[2];
    int l = idx[3];
    // Same coordinate swap --> symmetry only available if Hermitean
    if (isHermitian) {
        if (i < j)
            std::swap(i, j);
        if (k < l)
            std::swap(k, l);
    }
    // R12 swap
    if (i < k) {
        std::swap(i, k);
        std::swap(j, l);
    }
    if (i == k && j < l) {
        std::swap(j, l);
    }
    return std::array<int, 4>{i, j, k, l};
}

/** @brief Align function for the 6-element array */
inline auto alignArray(const std::array<int, 6>& idx, bool isHermitian) {
    int i = idx[0];
    int j = idx[1];
    int k = idx[2];
    int l = idx[3];
    int m = idx[4];
    int n = idx[5];
    // bra <--> ket symmetry for each coordinate, only valid if Hermitian
    if (isHermitian) {
        if (i < j)
            std::swap(i, j);
        if (k < l)
            std::swap(k, l);
        if (m < n)
            std::swap(m, n);
    }
    // R123 symmetry wrt 1/2
    if (i < k) {
        std::swap(i, k);
        std::swap(j, l);
    }
    // R123 symmetry wrt 2/3
    if (k < m) {
        std::swap(k, m);
        std::swap(l, n);
    }
    // R123 symmetry wrt 1/2 (to be repeated in case the previous swap screwed
    // up stuff)
    if (i < k) {
        std::swap(i, k);
        std::swap(j, l);
    }
    // Same as above, but for the specific case in which a few indices are equal.
    if (i == k && j < l) {
        std::swap(j, l);
    }
    if (k == m && l < n) {
        std::swap(l, n);
    }
    if (i == k && j < l) {
        std::swap(j, l);
    }
    return std::array<int, 6>{i, j, k, l, m, n};
}

/** @brief Same as above, but based on types. Alignes for real, not for complex */
template<class ScalarType>
class AlignTraitTypeClass { };

template<>
class AlignTraitTypeClass<double> {
public:
    static std::array<int, 4> align(const std::array<int, 4>& idx) { return alignArray(idx, true); }
    static std::array<int, 6> align(const std::array<int, 6>& idx) { return alignArray(idx, true); }
};

template<>
class AlignTraitTypeClass<std::complex<double> > {
public:
    static std::array<int, 4> align(const std::array<int, 4>& idx) { return idx; }
    static std::array<int, 6> align(const std::array<int, 6>& idx) { return idx; }
};

} // detail
} // maquis

#endif
