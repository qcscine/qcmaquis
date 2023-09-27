/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MAQUIS_BLOCK_MATRIX_DETAIL_AMBIENT_HPP
#define MAQUIS_BLOCK_MATRIX_DETAIL_AMBIENT_HPP

#include <ambient/container/numeric/enable_alps_hdf5.hpp>
#include <ambient/container/numeric/matrix.hpp>
#include "dmrg/utils/parallel.hpp"
#include "dmrg/block_matrix/detail/ambient_detail.hpp"
#include "dmrg/block_matrix/detail/alps.hpp"
#include "utils/traits.hpp"

namespace maquis { namespace traits {

    template<typename T, class Allocator> struct scalar_type <ambient::matrix<T,Allocator> > { typedef typename ambient::matrix<T,Allocator>::scalar_type type; };
    template<typename T> struct scalar_type <ambient::diagonal_matrix<T> > { typedef typename ambient::matrix<T>::scalar_type type; };
    template<typename T, class Allocator> struct real_type <ambient::matrix<T,Allocator> > { typedef typename ambient::matrix<T,Allocator>::real_type type; };
    template<typename T> struct real_type <ambient::diagonal_matrix<T> > { typedef typename ambient::matrix<T>::real_type type; };

    template<typename T, class Allocator> struct transpose_view< ambient::matrix<T,Allocator> > { typedef typename ambient::transpose_view<ambient::matrix<T,Allocator> > type; };
    template<typename T> struct transpose_view< ambient::diagonal_matrix<T> > { typedef ambient::diagonal_matrix<T> type; };

    template<class M> struct scalar_type <ambient::tiles<M> > { typedef typename scalar_type<M>::type type; };
    template<class M> struct real_type <ambient::tiles<M> > { typedef typename real_type<M>::type type; };
    template<class M> struct transpose_view< ambient::tiles<M> > { typedef typename ambient::tiles< typename transpose_view<M>::type > type; };

} }

namespace alps { namespace numeric {

    template<typename T, class Allocator> struct associated_vector< ambient::matrix<T,Allocator> > { typedef std::vector<T> type; };
    template<typename T, class Allocator> struct associated_one_matrix< ambient::matrix<T,Allocator> > { typedef ambient::matrix<T,Allocator> type; };
    template<typename T, class Allocator> struct associated_dense_matrix< ambient::matrix<T,Allocator> > { typedef ambient::matrix<T,Allocator> type; };
    template<typename T, class Allocator> struct associated_diagonal_matrix< ambient::matrix<T,Allocator> > { typedef ambient::diagonal_matrix<T> type; };
    template<typename T, class Allocator> struct associated_real_diagonal_matrix< ambient::matrix<T,Allocator> > { typedef ambient::diagonal_matrix<typename maquis::traits::real_type<T>::type> type; };
    template<typename T, class Allocator> struct associated_real_vector< ambient::matrix<T,Allocator> > { typedef std::vector<typename maquis::traits::real_type<T>::type> type; };

    template<class M> struct associated_vector< ambient::tiles<M> > { typedef ambient::tiles< typename associated_vector<M>::type > type; };
    template<class M> struct associated_real_vector< ambient::tiles<M> > { typedef ambient::tiles< typename associated_real_vector<M>::type > type; };
    template<class M> struct associated_one_matrix< ambient::tiles<M> > { typedef ambient::tiles<M> type; };
    template<class M> struct associated_dense_matrix< ambient::tiles<M> > { typedef ambient::tiles<M> type; };
    template<class M> struct associated_diagonal_matrix< ambient::tiles<M> > { typedef ambient::tiles< typename associated_diagonal_matrix<M>::type > type; };
    template<class M> struct associated_real_diagonal_matrix< ambient::tiles<M> > { typedef ambient::tiles< typename associated_real_diagonal_matrix<M>::type > type; };
    
} }

#endif
