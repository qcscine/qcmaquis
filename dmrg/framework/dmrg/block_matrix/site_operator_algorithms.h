/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2015 Institute for Theoretical Physics, ETH Zurich
 *               2015-2015 by Sebastian Keller <sebkelleb@phys.ethz.ch>
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

#ifndef SITE_OPERATOR_ALGORITHMS_H
#define SITE_OPERATOR_ALGORITHMS_H

//#include <boost/lambda/lambda.hpp>
//#include <boost/function.hpp>

#include "dmrg/utils/logger.h"
#include "dmrg/utils/utils.hpp"
#include "utils/timings.h"
#include "utils/traits.hpp"
#include "utils/bindings.hpp"

#include "dmrg/block_matrix/site_operator.h"

#include "dmrg/utils/parallel.hpp"

template<class Matrix1, class Matrix2, class Matrix3, class SymmGroup, class Scheduler>
void gemm(SiteOperator<Matrix1, SymmGroup> const & A,
          SiteOperator<Matrix2, SymmGroup> const & B,
          SiteOperator<Matrix3, SymmGroup> & C,
          const Scheduler& scheduler = Scheduler())
{
    C.clear();

    typedef typename SymmGroup::charge charge;
    Index<SymmGroup> B_left_basis = B.left_basis();
    for (std::size_t k = 0; k < A.n_blocks(); ++k) {
        std::size_t matched_block = B_left_basis.position(A.basis().right_charge(k));

        if ( matched_block == B.n_blocks() )
            continue;

        std::size_t new_block = C.insert_block(new Matrix3(num_rows(A[k]), num_cols(B[matched_block])),
                                               A.basis().left_charge(k), B.basis().right_charge(matched_block));

        parallel::guard proc(scheduler(k));
        gemm(A[k], B[matched_block], C[new_block]);
    }

    if(scheduler.propagate()){
        C.size_index.resize(C.n_blocks()); // propagating A size_index onto C - otherwise might C.index_sizes();
        for(size_t k = 0; k < A.n_blocks(); ++k){
            size_t matched_block = B_left_basis.position(A.basis().right_charge(k));
            if(matched_block != B.n_blocks())
                C.size_index(C.find_block(A.basis().left_charge(k), B.basis().right_charge(matched_block))) = A.size_index(k);
        }
    }
}

template<class Matrix1, class Matrix2, class Matrix3, class SymmGroup>
void gemm(SiteOperator<Matrix1, SymmGroup> const & A,
          SiteOperator<Matrix2, SymmGroup> const & B,
          SiteOperator<Matrix3, SymmGroup> & C)
{
    gemm(A, B, C, parallel::scheduler_nop());
}

template<class Matrix, class SymmGroup>
SiteOperator<Matrix, SymmGroup> adjoint(SiteOperator<Matrix, SymmGroup> m)
{
    m.adjoint_inplace();
    return m;
}

template<class Matrix, class SymmGroup>
SiteOperator<Matrix, SymmGroup> adjoin(SiteOperator<Matrix, SymmGroup> const & m) // error: it should be adjoin_t_
{
    SiteOperator<Matrix, SymmGroup> ret;
    for (std::size_t k = 0; k < m.n_blocks(); ++k)
        ret.insert_block(m[k],
                         -m.basis().left_charge(k),
                         -m.basis().right_charge(k));
    return ret;
}

template<class Matrix, class SymmGroup>
bool is_hermitian(SiteOperator<Matrix, SymmGroup> const & m)
{
    bool ret = true;
    for (size_t k=0; ret && k < m.n_blocks(); ++k) {
        if (m.basis().left_size(k) != m.basis().right_size(k))
            return false;
        else if (m.basis().left_charge(k) == m.basis().right_charge(k))
            ret = is_hermitian(m[k]);
        else if (! m.has_block(m.basis().right_charge(k), m.basis().left_charge(k)))
            return false;
        else
            ret = ( m[k] == transpose(conj( m(m.basis().right_charge(k), m.basis().left_charge(k)) )) );
    }
    return ret;
}

template <class Matrix, class SymmGroup, class A>
SiteOperator<Matrix, SymmGroup> op_exp_hermitian(Index<SymmGroup> const & phys,
                                                 SiteOperator<Matrix, SymmGroup> M,
                                                 A const & alpha = 1.)
{
    for (typename Index<SymmGroup>::const_iterator it_c = phys.begin(); it_c != phys.end(); it_c++)
        if (M.has_block(it_c->first, it_c->first))
            M(it_c->first, it_c->first) = exp_hermitian(M(it_c->first, it_c->first), alpha);
        else
            M.insert_block(Matrix::identity_matrix(phys.size_of_block(it_c->first)),
                           it_c->first, it_c->first);
    return M;
}

template <class OutOp, class Matrix, class SymmGroup, class A>
OutOp op_exp_hermitian(Index<SymmGroup> const & phys,
                       SiteOperator<Matrix, SymmGroup> const & M,
                       A const & alpha = 1.)
{
    OutOp ret(M.basis());
    for (typename Index<SymmGroup>::const_iterator it_c = phys.begin(); it_c != phys.end(); it_c++)
        if (M.has_block(it_c->first, it_c->first))
            ret(it_c->first, it_c->first) = exp_hermitian(M(it_c->first, it_c->first), alpha);
        else
            ret.insert_block(Matrix::identity_matrix(phys.size_of_block(it_c->first)),
                           it_c->first, it_c->first);
    return ret;
}

namespace detail {
    
    template <class Matrix>
    typename boost::enable_if<boost::is_complex<typename Matrix::value_type>, Matrix>::type
    exp_dispatcher(Matrix const& m, typename Matrix::value_type const& alpha)
    {
        return exp(m, alpha);
    }

    template <class Matrix>
    typename boost::disable_if<boost::is_complex<typename Matrix::value_type>, Matrix>::type
    exp_dispatcher(Matrix const& m, typename Matrix::value_type const& alpha)
    {
        throw std::runtime_error("Exponential of non-hermitian real matrices not implemented!");
        return Matrix();
    }
}

template <class Matrix, class SymmGroup, class A> SiteOperator<Matrix, SymmGroup> op_exp(Index<SymmGroup> const & phys,
                                       SiteOperator<Matrix, SymmGroup> M,
                                       A const & alpha = 1.)
{
    for (typename Index<SymmGroup>::const_iterator it_c = phys.begin(); it_c != phys.end(); it_c++)
        if (M.has_block(it_c->first, it_c->first))
            M(it_c->first, it_c->first) = detail::exp_dispatcher(M(it_c->first, it_c->first), alpha);
        else
            M.insert_block(Matrix::identity_matrix(phys.size_of_block(it_c->first)),
                           it_c->first, it_c->first);
    return M;
}

template<class Matrix1, class Matrix2, class SymmGroup>
void op_kron(Index<SymmGroup> const & phys_A,
             Index<SymmGroup> const & phys_B,
             SiteOperator<Matrix1, SymmGroup> const & A,
             SiteOperator<Matrix1, SymmGroup> const & B,
             SiteOperator<Matrix2, SymmGroup> & C,
             SpinDescriptor<typename symm_traits::SymmType<SymmGroup>::type> lspin
           = SpinDescriptor<typename symm_traits::SymmType<SymmGroup>::type>(),
             SpinDescriptor<typename symm_traits::SymmType<SymmGroup>::type> mspin
           = SpinDescriptor<typename symm_traits::SymmType<SymmGroup>::type>(),
             SpinDescriptor<typename symm_traits::SymmType<SymmGroup>::type> rspin
           = SpinDescriptor<typename symm_traits::SymmType<SymmGroup>::type>(),
             SpinDescriptor<typename symm_traits::SymmType<SymmGroup>::type> tspin
           = SpinDescriptor<typename symm_traits::SymmType<SymmGroup>::type>())
{
    C = SiteOperator<Matrix2, SymmGroup>();

    ProductBasis<SymmGroup> pb_left(phys_A, phys_B);
    ProductBasis<SymmGroup> const& pb_right = pb_left;

    for (int i = 0; i < A.n_blocks(); ++i) {
        for (int j = 0; j < B.n_blocks(); ++j) {
            typename SymmGroup::charge new_right = SymmGroup::fuse(A.basis().right_charge(i), B.basis().right_charge(j));
            typename SymmGroup::charge new_left = SymmGroup::fuse(A.basis().left_charge(i), B.basis().left_charge(j));


            Matrix2 tmp(pb_left.size(A.basis().left_charge(i), B.basis().left_charge(j)),
                       pb_right.size(A.basis().right_charge(i), B.basis().right_charge(j)),
                       0);

            maquis::dmrg::detail::op_kron(tmp, B[j], A[i],
                                          pb_left(A.basis().left_charge(i), B.basis().left_charge(j)),
                                          pb_right(A.basis().right_charge(i), B.basis().right_charge(j)),
                                          A.basis().left_size(i), B.basis().left_size(j),
                                          A.basis().right_size(i), B.basis().right_size(j));

            C.match_and_add_block(tmp, new_left, new_right);
        }
    }
}

namespace ts_ops_detail
{
    template <class Integer>
    std::vector<Integer> allowed_spins(Integer left, Integer right, Integer k1, Integer k2);
}

template<class Matrix1, class Matrix2, class SymmGroup>
void op_kron(Index<SymmGroup> const & phys_A,
             Index<SymmGroup> const & phys_B,
             SiteOperator<Matrix1, SymmGroup> const & Ao,
             SiteOperator<Matrix1, SymmGroup> const & Bo,
             SiteOperator<Matrix2, SymmGroup> & C,
             SpinDescriptor<symm_traits::SU2Tag> lspin,
             SpinDescriptor<symm_traits::SU2Tag> mspin,
             SpinDescriptor<symm_traits::SU2Tag> rspin,
             SpinDescriptor<symm_traits::SU2Tag> target_spin
              = SpinDescriptor<symm_traits::SU2Tag>(-1,0,0))
{
    typedef typename SymmGroup::charge charge;

    ProductBasis<SymmGroup> pb_left(phys_A, phys_B);
    ProductBasis<SymmGroup> const& pb_right = pb_left;

    C = SiteOperator<Matrix2, SymmGroup>();
    SiteOperator<Matrix1, SymmGroup> A = Ao, B = Bo;

    //*************************************
    // expand the small identity to the full one (Hack)

    if (A.spin().get() > 0 && B.spin().get() == 0)
    {
        charge cb = phys_B[1].first, cc = phys_B[2].first;
        if (!B.has_block(cb,cc))
        {
            B.insert_block(Matrix1(1,1,1), cb, cc);
            B.insert_block(Matrix1(1,1,1), cc, cb);
        }
    }
    if (A.spin().get() == 0 && B.spin().get() > 0)
    {
        charge cb = phys_A[1].first, cc = phys_A[2].first;

        if (!A.has_block(cb,cc))
        {
            A.insert_block(Matrix1(1,1,1), cb, cc);
            A.insert_block(Matrix1(1,1,1), cc, cb);
        }
    }

    //*************************************
    // MPO matrix basis spin QN's

    int k1 = A.spin().get(), k2 = B.spin().get(), k, j, jp, jpp;

    j = lspin.get();
    jpp = mspin.get();
    jp = rspin.get();

    std::vector<int> product_spins = ts_ops_detail::allowed_spins(j,jp, k1, k2);
    k = (target_spin.get() > -1) ? target_spin.get() : product_spins[0];

    //*************************************
    // Tensor + Kronecker product

    for (std::size_t i = 0; i < A.n_blocks(); ++i) {
        for (std::size_t j = 0; j < B.n_blocks(); ++j) {
            charge  inA = A.basis().left_charge(i);
            charge outA = A.basis().right_charge(i);
            charge  inB = B.basis().left_charge(j);
            charge outB = B.basis().right_charge(j);

            charge new_left = SymmGroup::fuse(inA, inB);
            charge new_right = SymmGroup::fuse(outA, outB);

            Matrix2 tmp(pb_left.size(inA, inB), pb_right.size(outA, outB), 0);

            std::size_t in_offset = pb_left(inA, inB);
            std::size_t out_offset = pb_right(outA, outB);

            maquis::dmrg::detail::op_kron(tmp, B[j], A[i], in_offset, out_offset,
                                          A.basis().left_size(i), B.basis().left_size(j),
                                          A.basis().right_size(i), B.basis().right_size(j));

            int j1  = std::abs(inA[1]),  j2  = std::abs(inB[1]),  J = productSpin<SymmGroup>(inA, inB);
            int j1p = std::abs(outA[1]), j2p = std::abs(outB[1]), Jp = productSpin<SymmGroup>(outA, outB);

            typename Matrix2::value_type coupling = SU2::mod_coupling(j1,j2,J,k1,k2,k,j1p,j2p,Jp);
            tmp *= coupling;

            C.match_and_add_block(tmp, new_left, new_right);
        }
    }

    //*************************************
    // Matrix basis coupling coefficient, applies uniformly to whole product

    typename Matrix2::value_type coupling = std::sqrt((jpp+1)*(k+1)) * gsl_sf_coupling_6j(j,jp,k,k2,k1,jpp);
    coupling = (((j+jp+k1+k2)/2)%2) ? -coupling : coupling;
    C *= coupling;

    SpinDescriptor<symm_traits::SU2Tag> op_spin(k, j, jp);
    C.spin() = op_spin;
}


// Currently not used

//template<class Matrix, class SymmGroup>
//void op_kron_long(MultiIndex<SymmGroup> const & midx,
//                  typename MultiIndex<SymmGroup>::set_id s,
//                  SiteOperator<Matrix, SymmGroup> const & A,
//                  SiteOperator<Matrix, SymmGroup> const & B,
//                  SiteOperator<Matrix, SymmGroup> const & F,
//                  std::size_t dist,
//                  SiteOperator<Matrix, SymmGroup> & C)
//{
//    assert( midx.size() == 2*(dist+1) );
//    C = SiteOperator<Matrix, SymmGroup>();
//    
//    for (size_t run=0; run<2; ++run) {
//        
//        if (run == 1)
//            C.allocate_blocks();
//        
//        for (index_product_iterator<SymmGroup> it = midx.begin();
//             it != midx.end(); ++it)
//        {
//            bool has_block = A.has_block((*it)[0].first, (*it)[1].first);
//            has_block = has_block && B.has_block((*it)[2*dist].first, (*it)[2*dist+1].first);
//            for (size_t i=1; has_block && i<dist; ++i)
//                has_block = F.has_block((*it)[2*i].first, (*it)[2*i+1].first);
//            
//            if (!has_block)
//                continue;
//            
//            typename Matrix::value_type val = A((*it)[0], (*it)[1]) * B((*it)[2*dist], (*it)[2*dist+1]);
//            for (size_t i=1; i<dist; ++i)
//                val *= F((*it)[2*i], (*it)[2*i+1]);
//            
//            if (val != 0.) {
//                typename MultiIndex<SymmGroup>::coord_t coord_l, coord_r;
//                boost::tie(coord_l, coord_r) = midx.get_coords(s, *it);
//                if (run == 0)
//                    C.reserve(coord_l.first, coord_r.first,
//                              midx.left_size(s, coord_l.first), midx.right_size(s, coord_r.first));
//                else
//                    C(coord_l, coord_r) += val;
//                
//            }
//        }
//        
//    }
//    
//}

#endif
