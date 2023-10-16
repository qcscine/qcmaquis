/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef SYMMETRY_TRAITS_H
#define SYMMETRY_TRAITS_H

#include <boost/type_traits.hpp>

namespace symm_traits {

// abelian vs. non-abelian

class AbelianTag {};
class SU2Tag {};

template <class SymmGroup>
struct SymmType { typedef AbelianTag type; };
template<>
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

template<>
struct Has2U1<TwoU1> : public std::true_type {};
template<>
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

/** @brief Trait class containing the symmetry name */
template<class SymmGroup>
class SymmetryNameTrait {};

// Various class specialization

template<>
class SymmetryNameTrait<TrivialGroup> {
public:
  static std::string symmName() {
    return "none";
  }
};

template<>
class SymmetryNameTrait<U1> {
public:
  static std::string symmName() {
    return "u1";
  }
};

template<int N>
class SymmetryNameTrait<NU1_template<N>> {
public:
  static const std::string symmName() {
    return "nu1";
  }
};

template<>
class SymmetryNameTrait<TwoU1> {
public:
  static const std::string symmName() {
    return "2u1";
  }
};

template<>
class SymmetryNameTrait<TwoU1PG> {
public:
  static const std::string symmName() {
    return "2u1pg";
  }
};

template<>
class SymmetryNameTrait<SU2U1> {
public:
  static const std::string symmName() {
    return "su2u1";
  }
};

template<>
class SymmetryNameTrait<SU2U1PG> {
public:
  static const std::string symmName() {
    return "su2u1pg";
  }
};

} // namespace symm_traits

#endif
