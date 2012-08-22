/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2012 by Alexandr Kosenkov <alex.kosenkov@gmail.com>
 *                            Timothee Ewart <timothee.ewart@gmail.com>
 *
 *****************************************************************************/

#ifndef MAQUIS_BLOCK_MATRIX_DETAIL_AMBIENT_HPP
#define MAQUIS_BLOCK_MATRIX_DETAIL_AMBIENT_HPP

#include <ambient/numeric/matrix.hpp>
#include "dmrg/block_matrix/detail/ambient_detail.hpp"
#include "utils/traits.hpp"

namespace maquis { namespace traits {

    template<typename T> struct scalar_type <ambient::numeric::matrix<T> > { typedef typename ambient::numeric::matrix<T>::scalar_type type; };
    template<typename T> struct scalar_type <ambient::numeric::diagonal_matrix<T> > { typedef typename ambient::numeric::matrix<T>::scalar_type type; };
    template<typename T> struct real_type <ambient::numeric::matrix<T> > { typedef typename ambient::numeric::matrix<T>::real_type type; };
    template<typename T> struct real_type <ambient::numeric::diagonal_matrix<T> > { typedef typename ambient::numeric::matrix<T>::real_type type; };

    template<typename T> struct transpose_view< ambient::numeric::matrix<T> > { typedef typename ambient::numeric::transpose_view<ambient::numeric::matrix<T> > type; };
    template<typename T> struct transpose_view< ambient::numeric::diagonal_matrix<T> > { typedef ambient::numeric::diagonal_matrix<T> type; };

    template<class M> struct scalar_type <ambient::numeric::tiles<M> > { typedef typename scalar_type<M>::type type; };
    template<class M> struct real_type <ambient::numeric::tiles<M> > { typedef typename real_type<M>::type type; };
    template<class M> struct transpose_view< ambient::numeric::tiles<M> > { typedef typename ambient::numeric::tiles< typename transpose_view<M>::type > type; };

} }

namespace alps { namespace numeric {

    template<typename T> struct associated_vector< ambient::numeric::matrix<T> > { typedef std::vector<T> type; };
    template<typename T> struct associated_diagonal_matrix< ambient::numeric::matrix<T> > { typedef ambient::numeric::diagonal_matrix<T> type; };
    template<typename T> struct associated_real_diagonal_matrix< ambient::numeric::matrix<T> > { typedef ambient::numeric::diagonal_matrix<typename maquis::traits::real_type<T>::type> type; };
    template<typename T> struct associated_real_vector< ambient::numeric::matrix<T> > { typedef std::vector<typename maquis::traits::real_type<T>::type> type; };

    template<class M> struct associated_vector< ambient::numeric::tiles<M> > { typedef ambient::numeric::tiles< typename associated_vector<M>::type > type; };
    template<class M> struct associated_real_vector< ambient::numeric::tiles<M> > { typedef ambient::numeric::tiles< typename associated_real_vector<M>::type > type; };
    template<class M> struct associated_diagonal_matrix< ambient::numeric::tiles<M> > { typedef ambient::numeric::tiles< typename associated_diagonal_matrix<M>::type > type; };
    template<class M> struct associated_real_diagonal_matrix< ambient::numeric::tiles<M> > { typedef ambient::numeric::tiles< typename associated_real_diagonal_matrix<M>::type > type; };
    
} }

#endif
