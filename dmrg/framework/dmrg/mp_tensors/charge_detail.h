/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2021 Institute for Theoretical Physics, ETH Zurich
 *               2021- by Alberto Baiardi <abaiardi@ethz.ch>
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

#ifndef CHARGE_DETAIL_H
#define CHARGE_DETAIL_H

#include "dmrg/block_matrix/symmetry.h"

/** 
 * @brief Charge details class 
 * 
 * Class containing symmetry specific functions.
 * For now, the following methods are implemented:
 * 
 *  - `physical`: checks if a given charge is *physical* so it represents a
 *     valid state independently on the specificity of the system.
 *  - `hasLessParticleThan`: compare two charges and checks if the first one 
 *     has a numbero of particles that is <= than the second one. Such a comparison
 *     is implemented (and makes sense) only for U1-like groups.
 */
template<class SymmGroup>
class ChargeDetailClass {
public:
    using ChargeType = typename SymmGroup::charge;
    static bool physical(ChargeType c) { return true; }
    static bool hasLessParticleThan(const ChargeType& inputCharge, const ChargeType& lowerCharge) {
        return true; 
    };
};

/** @brief U1 specialization */
template<>
class ChargeDetailClass<U1> {
public:
    using ChargeType = typename U1::charge;
    static bool physical(ChargeType c) { return c >= 0; }
    static bool hasLessParticleThan(const ChargeType& inputCharge, const ChargeType& refCharge) { return inputCharge <= refCharge; };
};

/** @brief TwoU1PG specialization */
template<>
class ChargeDetailClass<TwoU1PG> {
public:
    using ChargeType = typename TwoU1PG::charge;
    static bool physical(ChargeType c) { return c[0] >= 0 && c[1] >= 0; }
    static bool hasLessParticleThan(const ChargeType& inputCharge, const ChargeType& refCharge) { 
        return inputCharge[0] <= refCharge[0] && inputCharge[1] <= refCharge[1];
    };
};

/** @brief SU2U1 specialization */
template <>
class ChargeDetailClass<SU2U1> {
public:
    using ChargeType = typename SU2U1::charge;
    static bool physical(ChargeType c) { return SU2U1::spin(c) >= 0; }
    static bool hasLessParticleThan(const ChargeType& inputCharge, const ChargeType& lowerCharge) {
        return true;
    }
};

/** @brief SU2U1PG specialization */
template <>
class ChargeDetailClass<SU2U1PG> {
public:
    using ChargeType = SU2U1PG::charge;
    static bool physical(ChargeType c) { return SU2U1PG::spin(c) >= 0; }
    static bool hasLessParticleThan(const ChargeType& inputCharge, const ChargeType& lowerCharge) {
        return true;
    }
};


/** @brief NU1 specialization */
template <int N>
class ChargeDetailClass<NU1_template<N>> {
public:
    using ChargeType = typename NU1_template<N>::charge;
    static bool physical(ChargeType c) {
        return std::all_of(c.begin(), c.end(), [&](int tmp) { return (tmp >= 0); });
    }
    static bool hasLessParticleThan(const ChargeType& inputCharge, const ChargeType& lowerCharge) { 
        return std::mismatch(inputCharge.begin(), inputCharge.end(), lowerCharge.begin(), lowerCharge.end(),
                             std::less_equal<>{}) == std::make_pair(inputCharge.end(), lowerCharge.end());
    }
};

#endif // CHARGE_DETAIL_H