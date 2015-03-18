/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2012 by Alexandr Kosenkov <alex.kosenkov@gmail.com>
 *                            Timothee Ewart <timothee.ewart@gmail.com>
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

#ifndef MAQUIS_BLOCK_MATRIX_DETAIL_AMBIENT_HPP
#define MAQUIS_BLOCK_MATRIX_DETAIL_AMBIENT_HPP

#include <ambient/container/numeric/enable_alps_hdf5.hpp>
#include <ambient/container/numeric/matrix.hpp>
#include "dmrg/utils/parallel.hpp"
#include "dmrg/block_matrix/detail/ambient_detail.hpp"
#include "dmrg/block_matrix/detail/alps.hpp"
#include "utils/traits.hpp"

namespace maquis { namespace traits {

    template<typename T, class Allocator> struct scalar_type <ambient::numeric::matrix<T,Allocator> > { typedef typename ambient::numeric::matrix<T,Allocator>::scalar_type type; };
    template<typename T> struct scalar_type <ambient::numeric::diagonal_matrix<T> > { typedef typename ambient::numeric::matrix<T>::scalar_type type; };
    template<typename T, class Allocator> struct real_type <ambient::numeric::matrix<T,Allocator> > { typedef typename ambient::numeric::matrix<T,Allocator>::real_type type; };
    template<typename T> struct real_type <ambient::numeric::diagonal_matrix<T> > { typedef typename ambient::numeric::matrix<T>::real_type type; };

    template<typename T, class Allocator> struct transpose_view< ambient::numeric::matrix<T,Allocator> > { typedef typename ambient::numeric::transpose_view<ambient::numeric::matrix<T,Allocator> > type; };
    template<typename T> struct transpose_view< ambient::numeric::diagonal_matrix<T> > { typedef ambient::numeric::diagonal_matrix<T> type; };

    template<class M> struct scalar_type <ambient::numeric::tiles<M> > { typedef typename scalar_type<M>::type type; };
    template<class M> struct real_type <ambient::numeric::tiles<M> > { typedef typename real_type<M>::type type; };
    template<class M> struct transpose_view< ambient::numeric::tiles<M> > { typedef typename ambient::numeric::tiles< typename transpose_view<M>::type > type; };

} }

namespace alps { namespace numeric {

    template<typename T, class Allocator> struct associated_vector< ambient::numeric::matrix<T,Allocator> > { typedef std::vector<T> type; };
    template<typename T, class Allocator> struct associated_one_matrix< ambient::numeric::matrix<T,Allocator> > { typedef ambient::numeric::matrix<T,Allocator> type; };
    template<typename T, class Allocator> struct associated_dense_matrix< ambient::numeric::matrix<T,Allocator> > { typedef ambient::numeric::matrix<T,Allocator> type; };
    template<typename T, class Allocator> struct associated_diagonal_matrix< ambient::numeric::matrix<T,Allocator> > { typedef ambient::numeric::diagonal_matrix<T> type; };
    template<typename T, class Allocator> struct associated_real_diagonal_matrix< ambient::numeric::matrix<T,Allocator> > { typedef ambient::numeric::diagonal_matrix<typename maquis::traits::real_type<T>::type> type; };
    template<typename T, class Allocator> struct associated_real_vector< ambient::numeric::matrix<T,Allocator> > { typedef std::vector<typename maquis::traits::real_type<T>::type> type; };

    template<class M> struct associated_vector< ambient::numeric::tiles<M> > { typedef ambient::numeric::tiles< typename associated_vector<M>::type > type; };
    template<class M> struct associated_real_vector< ambient::numeric::tiles<M> > { typedef ambient::numeric::tiles< typename associated_real_vector<M>::type > type; };
    template<class M> struct associated_one_matrix< ambient::numeric::tiles<M> > { typedef ambient::numeric::tiles<M> type; };
    template<class M> struct associated_dense_matrix< ambient::numeric::tiles<M> > { typedef ambient::numeric::tiles<M> type; };
    template<class M> struct associated_diagonal_matrix< ambient::numeric::tiles<M> > { typedef ambient::numeric::tiles< typename associated_diagonal_matrix<M>::type > type; };
    template<class M> struct associated_real_diagonal_matrix< ambient::numeric::tiles<M> > { typedef ambient::numeric::tiles< typename associated_real_diagonal_matrix<M>::type > type; };
    
} }

#endif
