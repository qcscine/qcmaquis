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

namespace contraction { // Forward declaration

    template<class Matrix, class OtherMatrix, class SymmGroup>
    void lbtm_kernel(size_t b2,
                     ContractionGrid<Matrix, SymmGroup>& contr_grid,
                     Boundary<OtherMatrix, SymmGroup> const & left,
                     std::vector<block_matrix<Matrix, SymmGroup> > const & left_mult_mps,
                     MPOTensor<Matrix, SymmGroup> const & mpo,
                     DualIndex<SymmGroup> const & ket_basis,
                     Index<SymmGroup> const & right_i,
                     Index<SymmGroup> const & out_left_i,
                     ProductBasis<SymmGroup> const & in_right_pb,
                     ProductBasis<SymmGroup> const & out_left_pb);

    template<class Matrix, class OtherMatrix, class SymmGroup>
    void rbtm_kernel(size_t b1,
                     block_matrix<Matrix, SymmGroup> & ret,
                     Boundary<OtherMatrix, SymmGroup> const & right,
                     std::vector<block_matrix<Matrix, SymmGroup> > const & right_mult_mps,
                     MPOTensor<Matrix, SymmGroup> const & mpo,
                     DualIndex<SymmGroup> const & ket_basis,
                     Index<SymmGroup> const & left_i,
                     Index<SymmGroup> const & out_right_i,
                     ProductBasis<SymmGroup> const & in_left_pb,
                     ProductBasis<SymmGroup> const & out_right_pb);
}

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

struct AbelianGemms
{
    typedef gemm_functor gemm;
    typedef gemm_trim_left_functor gemm_trim_left;
    typedef gemm_trim_right_functor gemm_trim_right;
};

namespace contraction {

    struct lbtm_functor
    {
        template<class Matrix, class OtherMatrix, class SymmGroup>
        void operator()(size_t b2,
                        contraction::ContractionGrid<Matrix, SymmGroup>& contr_grid,
                        Boundary<OtherMatrix, SymmGroup> const & left,
                        std::vector<block_matrix<Matrix, SymmGroup> > const & left_mult_mps,
                        MPOTensor<Matrix, SymmGroup> const & mpo,
                        DualIndex<SymmGroup> const & ket_basis,
                        Index<SymmGroup> const & right_i,
                        Index<SymmGroup> const & out_left_i,
                        ProductBasis<SymmGroup> const & in_right_pb,
                        ProductBasis<SymmGroup> const & out_left_pb)
        {
            contraction::lbtm_kernel(b2, contr_grid, left, left_mult_mps, mpo, ket_basis,
                                     right_i, out_left_i, in_right_pb, out_left_pb);
        }
    };

    struct rbtm_functor
    {
        template<class Matrix, class OtherMatrix, class SymmGroup>
        void operator()(size_t b1,
                        block_matrix<Matrix, SymmGroup> & ret,
                        Boundary<OtherMatrix, SymmGroup> const & right,
                        std::vector<block_matrix<Matrix, SymmGroup> > const & right_mult_mps,
                        MPOTensor<Matrix, SymmGroup> const & mpo,
                        DualIndex<SymmGroup> const & ket_basis,
                        Index<SymmGroup> const & left_i,
                        Index<SymmGroup> const & out_right_i,
                        ProductBasis<SymmGroup> const & in_left_pb,
                        ProductBasis<SymmGroup> const & out_right_pb)
        {
            return contraction::rbtm_kernel(b1, ret, right, right_mult_mps, mpo, ket_basis,
                                            left_i, out_right_i, in_left_pb, out_right_pb);
        }

    };

} // namespace contraction

#endif
