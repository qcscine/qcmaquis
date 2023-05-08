/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef CONTRACTIONS_SPECIAL_HPP
#define CONTRACTIONS_SPECIAL_HPP

#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/mpotensor.h"

namespace contraction {
   
    template<class Matrix, class OtherMatrix, class SymmGroup>
    MPSTensor<Matrix, SymmGroup>
    site_ortho_boundaries(MPSTensor<Matrix, SymmGroup> const & mps,
                          MPSTensor<Matrix, SymmGroup> const & ortho_mps,
                          block_matrix<OtherMatrix, SymmGroup> const & ortho_left,
                          block_matrix<OtherMatrix, SymmGroup> const & ortho_right)
    {
        ortho_mps.make_right_paired();
        block_matrix<Matrix, SymmGroup> t, t2, t3;
        gemm(ortho_left, ortho_mps.data(), t);
        reshape_right_to_left_new(mps.site_dim(),
                                  ortho_left.left_basis(), ortho_mps.col_dim(),
                                  t, t2);
        gemm(t2, transpose(ortho_right), t3);
        
        mps.make_left_paired();
        t = mps.data();
        reshape_and_pad_left(mps.site_dim(),
                             ortho_left.left_basis(), ortho_right.left_basis(),
                             mps.row_dim(), mps.col_dim(),
                             t3, t);
        
        MPSTensor<Matrix, SymmGroup> t4(mps.site_dim(),
                                        mps.row_dim(), mps.col_dim(),
                                        t, LeftPaired);
        return t4;
    }
    
    template<class Matrix, class SymmGroup>
    block_matrix<Matrix, SymmGroup>
    multiply_with_twosite(block_matrix<Matrix, SymmGroup> const & vec,
                          typename operator_selector<Matrix, SymmGroup>::type const & op,
                          Index<SymmGroup> const & left_i,
                          Index<SymmGroup> const & right_i,
                          Index<SymmGroup> const & phys_i)
    {
        typedef typename SymmGroup::charge charge;
        
        Index<SymmGroup> l_index = phys_i * left_i, r_index = adjoin(phys_i) * right_i;
        common_subset(l_index, r_index);
        block_matrix<Matrix, SymmGroup> ret;
        
        ProductBasis<SymmGroup> left_pb(phys_i, left_i);
        ProductBasis<SymmGroup> right_pb(phys_i, right_i,
                                         boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                             -boost::lambda::_1, boost::lambda::_2));
        ProductBasis<SymmGroup> phys_pb(phys_i, phys_i);
        
        for (size_t ls = 0; ls < left_i.size(); ++ls)
            for (size_t rs = 0; rs < right_i.size(); ++rs)
                for (size_t lps = 0; lps < phys_i.size(); ++lps)
                    for (size_t rps = 0; rps < phys_i.size(); ++rps)
                        for (size_t ilps = 0; ilps < phys_i.size(); ++ilps)
                            for (size_t irps = 0; irps < phys_i.size(); ++irps)
                            {
                                charge lc = left_i[ls].first, rc = right_i[rs].first;
                                charge lpc = phys_i[lps].first, rpc = phys_i[rps].first;
                                charge ilpc = phys_i[ilps].first, irpc = phys_i[irps].first;
                                
                                charge left_vec_charge = SymmGroup::fuse(ilpc, lc);
                                charge right_vec_charge = SymmGroup::fuse(-irpc, rc);
                                
                                charge left_out_charge = SymmGroup::fuse(lpc, lc);
                                charge right_out_charge = SymmGroup::fuse(-rpc, rc);
                                
                                if (left_out_charge != right_out_charge)
                                    continue;
                                if (left_vec_charge != right_vec_charge)
                                    continue;
                                if (! vec.has_block(left_vec_charge, right_vec_charge) )
                                    continue;
                                if (SymmGroup::fuse(lpc, rpc) != SymmGroup::fuse(ilpc, irpc))
                                    continue;

                                if (! ret.has_block(left_out_charge, right_out_charge))
                                    ret.insert_block(Matrix(l_index.size_of_block(left_out_charge), r_index.size_of_block(right_out_charge), 0),
                                                     left_out_charge, right_out_charge);

                                charge both_charge = SymmGroup::fuse(lpc, rpc);
                                
                                assert( op.has_block(both_charge, both_charge) );
                                
                                maquis::dmrg::detail::mwt(ret(left_out_charge, right_out_charge), 
                                                          vec(left_vec_charge, right_vec_charge), 
                                                          op(both_charge, both_charge), 
                                                          left_pb(lpc, lc), right_pb(rpc, rc),
                                                          left_pb(ilpc, lc), right_pb(irpc, rc),
                                                          phys_pb(ilpc, irpc), phys_pb(lpc, rpc),
                                                          left_i[ls].second, right_i[rs].second,
                                                          phys_i[lps].second, phys_i[rps].second,
                                                          phys_i[ilps].second, phys_i[irps].second);
                            }
        
        return ret;
    }
    
    template<class Matrix, class SymmGroup>
    MPSTensor<Matrix, SymmGroup>
    multiply_with_op(MPSTensor<Matrix, SymmGroup> const & mps, typename operator_selector<Matrix, SymmGroup>::type const & op)
    {
        typedef typename SymmGroup::charge charge;
        
        mps.make_left_paired();
        block_matrix<Matrix, SymmGroup> const & vec = mps.data();
        
        Index<SymmGroup> const& phys_i  = mps.site_dim();
        Index<SymmGroup> const& left_i  = mps.row_dim();
        Index<SymmGroup> const& right_i = mps.col_dim();
        
        ProductBasis<SymmGroup> left_pb(phys_i, left_i);
        
        block_matrix<Matrix, SymmGroup> ret_vec;
        Index<SymmGroup> out_left_i;

        for (size_t ls = 0; ls < left_i.size(); ++ls)
            for (size_t rs = 0; rs < right_i.size(); ++rs)
                for (size_t s1 = 0; s1 < phys_i.size(); ++s1)
                    for (size_t s2 = 0; s2 < phys_i.size(); ++s2) {
                        
                        charge lc  = left_i[ls].first, rc  = right_i[rs].first;
                        charge s1c = phys_i[s1].first, s2c = phys_i[s2].first;
                        
                        if (! op.has_block(s1c, s2c) )
                            continue;

                        charge phys_diff = SymmGroup::fuse(s2c, -s1c);
                        
                        charge left_vec_charge  = SymmGroup::fuse(s1c, lc);
                        charge right_vec_charge = rc;

                        if (! vec.has_block(left_vec_charge, right_vec_charge) )
                            continue;
                        
                        if (! out_left_i.has(lc))
                            out_left_i.insert(left_i[ls]);
                        
                        charge left_out_charge  = SymmGroup::fuse(s2c, lc);
                        charge right_out_charge = SymmGroup::fuse(rc, phys_diff);
                        
                        if (! ret_vec.has_block(left_out_charge, right_out_charge) )
                            ret_vec.insert_block(new Matrix(left_pb.size(s2c, lc), right_i[rs].second, 0), left_out_charge, right_out_charge);
                        
                        Matrix & oblock         = ret_vec(left_out_charge, right_out_charge);
                        Matrix const & iblock   = vec(left_vec_charge, right_vec_charge);
                        Matrix const & op_block = op(s1c, s2c);
                        
                        size_t i_l_offset = left_pb(s1c, lc);
                        size_t i_r_offset = 0;
                        size_t l_offset = left_pb(s2c, lc);
                        size_t r_offset = 0;
                        size_t i_op_offset = 0;
                        size_t op_offset = 0;
                        
#ifdef USE_AMBIENT
                        printf("UNOPTIMIZED FUNCTION (MULTIPLY WITH OP!)\n");
#endif
                        for (size_t ss1 = 0; ss1 < phys_i[s1].second; ++ss1) {
                            for (size_t ss2 = 0; ss2 < phys_i[s2].second; ++ss2) {
                                size_t o_l_start = l_offset + ss2*left_i[ls].second;
                                size_t o_r_start = r_offset;
                                
                                size_t i_l_start = i_l_offset + ss1*left_i[ls].second;
                                size_t i_r_start = i_r_offset;
                                
                                typename Matrix::value_type const & op_val = op_block(i_op_offset + ss1,  op_offset + ss2);
                                
                                // TODO: replace by kernel
                                for (size_t rr = 0; rr < right_i[rs].second; ++rr) {
                                    for (size_t ll = 0; ll < left_i[ls].second; ++ll) {
                                        oblock(o_l_start+ll, o_r_start+rr) += iblock(i_l_start+ll, i_r_start+rr) * op_val;
                                    }
                                }
                            }
                        }
                    }
        
        
        MPSTensor<Matrix, SymmGroup> ret(phys_i, out_left_i, ret_vec.right_basis(), ret_vec, LeftPaired);
        assert( ret.reasonable() );
        return ret;
    }
    
    // tested only for diagonal operator
    template<class Matrix, class SymmGroup>
    MPSTensor<Matrix, SymmGroup>
    local_op(MPSTensor<Matrix, SymmGroup> const & mps,
             typename operator_selector<Matrix, SymmGroup>::type const & op)
    {
        typedef typename SymmGroup::charge charge;
        typedef typename operator_selector<Matrix, SymmGroup>::type op_t;
        
        mps.make_left_paired();
        block_matrix<Matrix, SymmGroup> const & vec = mps.data();
        
        Index<SymmGroup> phys_i = mps.site_dim();
        Index<SymmGroup> left_i = mps.row_dim();
        Index<SymmGroup> right_i = mps.col_dim();
        
        ProductBasis<SymmGroup> left_pb(phys_i, left_i);
        
        MPSTensor<Matrix, SymmGroup> ret(phys_i, left_i, right_i, false);
        block_matrix<Matrix, SymmGroup> & ret_vec = ret.data();
        
        for (size_t ls = 0; ls < left_i.size(); ++ls)
            for (size_t rs = 0; rs < right_i.size(); ++rs)
                for (size_t s1 = 0; s1 < phys_i.size(); ++s1)
                    for (size_t s2 = 0; s2 < phys_i.size(); ++s2) {
                        
                        charge lc = left_i[ls].first, rc = right_i[rs].first;
                        charge s1c = phys_i[s1].first, s2c = phys_i[s2].first;
                        
                        charge left_vec_charge = SymmGroup::fuse(s1c, lc);
                        charge right_vec_charge = rc;
                        charge left_out_charge = SymmGroup::fuse(s2c, lc);
                        charge right_out_charge = right_vec_charge;
                        
                        if (! vec.has_block(left_vec_charge, right_vec_charge) )
                            continue;
                        if (! op.has_block(s1c, s2c) )
                            continue;
                        if (! ret_vec.has_block(left_out_charge, right_out_charge) )
                            ret_vec.insert_block(new Matrix(left_pb.size(s2c, lc), right_i[rs].second, 0), left_out_charge, right_out_charge);
                        
                        Matrix & oblock = ret_vec(left_out_charge, right_out_charge);
                        Matrix const & iblock = vec(left_vec_charge, right_vec_charge);
                        Matrix const & op_block = op(s1c, s2c);
                        
                        size_t i_l_offset = left_pb(s1c, lc);
                        size_t i_r_offset = 0;
                        size_t l_offset = left_pb(s2c, lc);
                        size_t r_offset = 0;
                        size_t i_op_offset = 0;
                        size_t op_offset = 0;
                        
                        for (size_t ll = 0; ll < left_i[ls].second; ++ll)
                            for (size_t rr = 0; rr < right_i[rs].second; ++rr)
                                for (size_t ss1 = 0; ss1 < phys_i[s1].second; ++ss1)
                                    for (size_t ss2 = 0; ss2 < phys_i[s2].second; ++ss2) {
#ifndef NDEBUG
                                        oblock(l_offset + ss2*left_i[ls].second + ll,
                                               r_offset + rr);
                                        iblock(i_l_offset + ss1*left_i[ls].second + ll,
                                               i_r_offset + rr);
                                        op_block(i_op_offset + ss1,
                                                 op_offset + ss2);
#endif
                                        oblock(l_offset + ss2*left_i[ls].second + ll,
                                               r_offset + rr) += 
                                        iblock(i_l_offset + ss1*left_i[ls].second + ll,
                                               i_r_offset + rr) *
                                        op_block(i_op_offset + ss1,
                                                 op_offset + ss2);
                                        
                                    }
                    }
        
        
        assert( ret.reasonable() );
        return ret;
    }

    template<class Matrix, class SymmGroup>
    block_matrix<Matrix, SymmGroup>
    density_matrix(MPSTensor<Matrix, SymmGroup> const & bra_tensor,
                   MPSTensor<Matrix, SymmGroup> const & ket_tensor)
    {
        typedef typename SymmGroup::charge charge;
        
        assert( bra_tensor.row_dim() == ket_tensor.row_dim() );
        assert( bra_tensor.col_dim() == ket_tensor.col_dim() );
        assert( bra_tensor.site_dim() == ket_tensor.site_dim() );
        
        bra_tensor.make_left_paired();
        ket_tensor.make_left_paired();
        
        Index<SymmGroup> phys_i = ket_tensor.site_dim();
        Index<SymmGroup> left_i = ket_tensor.row_dim();
        Index<SymmGroup> right_i = ket_tensor.col_dim();
        
        ProductBasis<SymmGroup> left_pb(phys_i, left_i);
        
        block_matrix<Matrix, SymmGroup> const & bra_vec = bra_tensor.data();
        block_matrix<Matrix, SymmGroup> const & ket_vec = ket_tensor.data();
        
        block_matrix<Matrix, SymmGroup> ret(phys_i, phys_i);
        
        for (size_t ls = 0; ls < left_i.size(); ++ls)
            for (size_t rs = 0; rs < right_i.size(); ++rs)
                for (size_t s1 = 0; s1 < phys_i.size(); ++s1)
                    for (size_t s2 = 0; s2 < phys_i.size(); ++s2) {
                        
                        charge lc = left_i[ls].first, rc = right_i[rs].first;
                        charge s1c = phys_i[s1].first, s2c = phys_i[s2].first;
                        
                        charge left_ket_charge = SymmGroup::fuse(s1c, lc);
                        charge left_bra_charge = SymmGroup::fuse(s2c, lc);
                        charge right_charge = rc;
                        
                        if (! ket_vec.has_block(left_ket_charge, right_charge) )
                            continue;
                        if (! bra_vec.has_block(left_bra_charge, right_charge) )
                            continue;
                        if (! ret.has_block(s1c, s2c) )
                            ret.insert_block(new Matrix(phys_i[s1].second, phys_i[s2].second, 0), s1c, s2c);
                        
                        Matrix & oblock = ret(s1c, s2c);
                        Matrix const & ket_block = ket_vec(left_ket_charge, right_charge);
                        Matrix const & bra_block = bra_vec(left_bra_charge, right_charge);
                        
                        size_t l_ket_offset = left_pb(s1c, lc);
                        size_t l_bra_offset = left_pb(s2c, lc);
                        
                        for (size_t ll = 0; ll < left_i[ls].second; ++ll)
                            for (size_t rr = 0; rr < right_i[rs].second; ++rr)
                                for (size_t ss1 = 0; ss1 < phys_i[s1].second; ++ss1)
                                    for (size_t ss2 = 0; ss2 < phys_i[s2].second; ++ss2) {
                                        
                                        oblock(ss1, ss2) +=
                                        utils::conj(bra_block(l_bra_offset + ss2*left_i[ls].second + ll, rr)) *
                                        ket_block(l_ket_offset + ss1*left_i[ls].second + ll, rr);
                                        
                                    }
                    }
        return ret;
    }

    template<class Matrix, class SymmGroup>
    block_matrix<Matrix, SymmGroup>
    density_matrix_2(MPSTensor<Matrix, SymmGroup> const & bra_tensor,
                     MPSTensor<Matrix, SymmGroup> const & ket_tensor)
    {
        typedef typename SymmGroup::charge charge;
        
        assert( bra_tensor.row_dim() == ket_tensor.row_dim() );
        assert( bra_tensor.col_dim() == ket_tensor.col_dim() );
        assert( bra_tensor.site_dim() == ket_tensor.site_dim() );
        
        bra_tensor.make_left_paired();
        ket_tensor.make_left_paired();
        
        Index<SymmGroup> phys_i = ket_tensor.site_dim();
        Index<SymmGroup> left_i = ket_tensor.row_dim();
        Index<SymmGroup> right_i = ket_tensor.col_dim();
        
        block_matrix<Matrix, SymmGroup> ket_mat;
        reshape_left_to_physleft(phys_i, left_i, right_i, ket_tensor.data(), ket_mat);
        block_matrix<Matrix, SymmGroup> bra_mat;
        reshape_left_to_physleft(phys_i, left_i, right_i, bra_tensor.data(), bra_mat);

        block_matrix<Matrix, SymmGroup> dm;
        gemm(ket_mat, transpose(conjugate(bra_mat)), dm);
        
        return dm;
    }

}

#endif
