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

#ifndef CONTRACTION_ABELIAN_FUNCTORS
#define CONTRACTION_ABELIAN_FUNCTORS


#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"

namespace Abelian {

    namespace detail {

        struct gemm_functor
        {
            template<class Matrix1, class Matrix2, class Matrix3, class SymmGroup>
            void operator()(block_matrix<Matrix1, SymmGroup> const & A,
                            block_matrix<Matrix2, SymmGroup> const & B,
                            block_matrix<Matrix3, SymmGroup> & C)
            {
                gemm(A,B,C);
            }
        };

        struct gemm_trim_left_functor
        {
            template<class Matrix1, class Matrix2, class Matrix3, class SymmGroup>
            void operator()(block_matrix<Matrix1, SymmGroup> const & A,
                            block_matrix<Matrix2, SymmGroup> const & B,
                            block_matrix<Matrix3, SymmGroup> & C)
            {
                gemm_trim_left(A,B,C);
            }
        };

        struct gemm_trim_right_functor
        {
            template<class Matrix1, class Matrix2, class Matrix3, class SymmGroup>
            void operator()(block_matrix<Matrix1, SymmGroup> const & A,
                            block_matrix<Matrix2, SymmGroup> const & B,
                            block_matrix<Matrix3, SymmGroup> & C)
            {
                gemm_trim_right(A,B,C);
            }
        };

    } // namespace detail

    struct Gemms
    {
        typedef detail::gemm_functor gemm;
        typedef detail::gemm_trim_left_functor gemm_trim_left;
        typedef detail::gemm_trim_right_functor gemm_trim_right;
    };

} // namespace abelian

#endif
