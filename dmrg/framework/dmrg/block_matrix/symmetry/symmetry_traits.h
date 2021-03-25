/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2014-2014 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

#ifndef SYMMETRY_TRAITS_H
#define SYMMETRY_TRAITS_H

#include <boost/type_traits.hpp>

namespace symm_traits {

    // abelian vs. non-abelian

    class AbelianTag {};
    class SU2Tag {};

    template <class SymmGroup>
    struct SymmType { typedef AbelianTag type; };

    template <>
    struct SymmType<SU2U1> { typedef SU2Tag type; };

    template <>
    struct SymmType<SU2U1PG> { typedef SU2Tag type; };



    template <class SymmGroup>
    struct HasU1DG : public std::false_type {};

    template <>
    struct HasU1DG<U1DG> : public std::true_type {};

    template<class SymmGroup, class T=void>
    using enable_if_u1dg_t = typename std::enable_if<symm_traits::HasU1DG<SymmGroup>::value, T>::type;
    template<class SymmGroup, class T=void>
    using disable_if_u1dg_t = typename std::enable_if<!symm_traits::HasU1DG<SymmGroup>::value, T>::type;

    template <class SymmGroup>
    struct Has2U1 : public std::false_type {};

    template <>
    struct Has2U1<TwoU1> : public std::true_type {};
    template <>
    struct Has2U1<TwoU1PG> : public std::true_type {};


    template <class SymmGroup>
    struct HasSU2 : public std::false_type {};

    template <>
    struct HasSU2<SU2U1> : public std::true_type {};
    template <>
    struct HasSU2<SU2U1PG> : public std::true_type {};

    template<class SymmGroup, class T=void>
    using enable_if_su2_t = typename std::enable_if<symm_traits::HasSU2<SymmGroup>::value, T>::type;
    template<class SymmGroup, class T=void>
    using disable_if_su2_t = typename std::enable_if<!symm_traits::HasSU2<SymmGroup>::value, T>::type;

    // point group vs. no point group

    template <class SymmGroup>
    struct HasPG : public std::false_type {};

    template <>
    struct HasPG<TwoU1PG> : public std::true_type {};
    template <>
    struct HasPG<SU2U1PG> : public std::true_type {};
    template <>
    struct HasPG<U1DG> : public std::true_type {};

    template<class SymmGroup, class T=void>
    using enable_if_pg_t = typename std::enable_if<symm_traits::HasPG<SymmGroup>::value, T>::type;
    template<class SymmGroup, class T=void>
    using disable_if_pg_t = typename std::enable_if<!symm_traits::HasPG<SymmGroup>::value, T>::type;


    // chemistry model implemented or not

    template <class SymmGroup>
    struct HasChemModel : public std::false_type {};

    template <>
    struct HasChemModel<TwoU1> : public std::true_type {};
    template <>
    struct HasChemModel<TwoU1PG> : public std::true_type {};
    template <>
    struct HasChemModel<SU2U1> : public std::true_type {};
    template <>
    struct HasChemModel<SU2U1PG> : public std::true_type {};
    template <>
    struct HasChemModel<U1DG> : public std::true_type {};

    template<class SymmGroup, class T=void>
    using enable_if_chemmodel_t = typename std::enable_if<symm_traits::HasChemModel<SymmGroup>::value, T>::type;
    template<class SymmGroup, class T=void>
    using disable_if_chemmodel_t = typename std::enable_if<!symm_traits::HasChemModel<SymmGroup>::value, T>::type;

}

#endif
