#ifndef CONTRACTIONS_H
#define CONTRACTIONS_H

#include "mp_tensors/mpstensor.h"
#include "mp_tensors/mpotensor.h"

#include "block_matrix/reshapes.h"
#include "block_matrix/indexing.h"

struct contraction {
    template<class Matrix, class SymmGroup>
    static block_matrix<Matrix, SymmGroup>
    overlap_left_step(MPSTensor<Matrix, SymmGroup> const & bra_tensor,
                      MPSTensor<Matrix, SymmGroup> const & ket_tensor,
                      block_matrix<Matrix, SymmGroup> const & left,
                      block_matrix<Matrix, SymmGroup> * localop = NULL)
    {
        if (localop != NULL)
            throw std::runtime_error("Not implemented!");
        
        assert(left.left_basis() == bra_tensor.left_i);
        assert(left.right_basis() == ket_tensor.left_i);
        
        bra_tensor.make_left_paired();
        ket_tensor.make_right_paired();
        
        block_matrix<Matrix, SymmGroup> t1, t2 = conjugate(bra_tensor.data_), t3;
        gemm(left, ket_tensor.data_, t1);
        reshape_right_to_left(ket_tensor.phys_i, left.left_basis(), ket_tensor.right_i,
                              t1, t3);
        t3 = transpose(t3);
        gemm(t3, t2, t1);
        return transpose(t1);
    }
    
    template<class Matrix, class SymmGroup>
    static block_matrix<Matrix, SymmGroup>
    overlap_right_step(MPSTensor<Matrix, SymmGroup> const & bra_tensor,
                       MPSTensor<Matrix, SymmGroup> const & ket_tensor,
                       block_matrix<Matrix, SymmGroup> const & right,
                       block_matrix<Matrix, SymmGroup> * localop = NULL)
    {
        if (localop != NULL)
            throw std::runtime_error("Not implemented!");
        
        assert(right.left_basis() == bra_tensor.right_i);
        assert(right.right_basis() == ket_tensor.right_i);
        
        bra_tensor.make_right_paired();
        ket_tensor.make_left_paired();
        
        block_matrix<Matrix, SymmGroup> t1, t2 = conjugate(transpose(bra_tensor.data_)), t3 = transpose(right);
        gemm(ket_tensor.data_, t3, t1);
        reshape_left_to_right(ket_tensor.phys_i, ket_tensor.left_i, right.left_basis(),
                              t1, t3);
        
        gemm(t3, t2, t1);
        return transpose(t1);
    }
    
    template<class Matrix, class SymmGroup>
    static Boundary<Matrix, SymmGroup>
    overlap_mpo_left_step(MPSTensor<Matrix, SymmGroup> const & bra_tensor,
                          MPSTensor<Matrix, SymmGroup> const & ket_tensor,
                          Boundary<Matrix, SymmGroup> const & left,
                          MPOTensor<Matrix, SymmGroup> const & mpo)
    {
        assert(left.aux_dim() == mpo.row_dim());
        
        bra_tensor.make_right_paired();
        ket_tensor.make_right_paired();
        
        block_matrix<Matrix, SymmGroup> t1;
        std::vector<block_matrix<Matrix, SymmGroup> > t2(left.aux_dim());
        
        for (std::size_t b = 0; b < left.aux_dim(); ++b)
        {
            gemm(transpose(ket_tensor.data_), left.data_[b], t1);
            gemm(t1, conjugate(bra_tensor.data_), t2[b]);
        }
        
        Boundary<Matrix, SymmGroup> ret;
        ret.data_.resize(mpo.col_dim());
        
        typedef typename SymmGroup::charge charge;
        typedef std::size_t size_t;
        
        for (size_t b1 = 0; b1 < mpo.row_dim(); ++b1)
            for (size_t b2 = 0; b2 < mpo.col_dim(); ++b2)
            {
                block_matrix<Matrix, SymmGroup> const & T = t2[b1];
                block_matrix<Matrix, SymmGroup> const & W = mpo(b1, b2);
                
                if (T.n_blocks() == 0 || W.n_blocks() == 0)
                    continue;
                
                // the boost::bind just turns around the physical charges
                // cf the definition of s1_size, s2_size below
                ProductBasis<SymmGroup> upper_pb(ket_tensor.site_dim(), ket_tensor.col_dim(),
                                                 boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                                     -boost::lambda::_1, boost::lambda::_2)
                                                 );
                ProductBasis<SymmGroup> lower_pb(bra_tensor.site_dim(), bra_tensor.col_dim(),
                                                 boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                                     -boost::lambda::_1, boost::lambda::_2)
                                                 );
                
                for (size_t s1 = 0; s1 < ket_tensor.site_dim().size(); ++s1)
                    for (size_t u = 0; u < ket_tensor.col_dim().size(); ++u)
                    {
                        charge u_charge = ket_tensor.col_dim()[u].first;
                        charge s1_charge = ket_tensor.site_dim()[s1].first;
                        charge tu_charge = SymmGroup::fuse(-s1_charge, u_charge);
                        
                        size_t upper_size = ket_tensor.col_dim()[u].second;
                        size_t s1_size = ket_tensor.site_dim()[s1].second;
                        
                        for (size_t s2 = 0; s2 < bra_tensor.site_dim().size(); ++s2)
                            for (size_t l = 0; l < bra_tensor.col_dim().size(); ++l)
                            {   
                                charge l_charge = bra_tensor.col_dim()[l].first;
                                charge s2_charge = bra_tensor.site_dim()[s2].first;                                
                                charge tl_charge = SymmGroup::fuse(-s2_charge, l_charge);
                                
                                size_t lower_size = bra_tensor.col_dim()[l].second;
                                size_t s2_size = bra_tensor.site_dim()[s2].second;
                                
                                Matrix block(upper_size, lower_size);
                                block *= 0;
                                
                                size_t in_l_offset = upper_pb(s1_charge, u_charge);
                                size_t in_r_offset = lower_pb(s2_charge, l_charge);
                                
                                if (! W.has_block(s1_charge, s2_charge) )
                                    continue;
                                
                                if (! T.has_block(tu_charge, tl_charge) )
                                    continue;
                                
                                for (size_t ss1 = 0; ss1 < s1_size; ++ss1)
                                    for (size_t ss2 = 0; ss2 < s2_size; ++ss2)
                                        for (size_t uu = 0; uu < upper_size; ++uu)
                                            for (size_t ll = 0; ll < lower_size; ++ll)
                                                block(uu, ll) += W(s1_charge, s2_charge)(ss1, ss2) *
                                                T(tu_charge, tl_charge)(in_l_offset+ss1*upper_size+uu,
                                                                        in_r_offset+ss2*lower_size+ll);
                                
                                ret.data_[b2] += block_matrix<Matrix, SymmGroup>(u_charge, l_charge, block);
                            }
                    }
            }
        
        return ret;
    }
    
    template<class Matrix, class SymmGroup>
    static Boundary<Matrix, SymmGroup>
    overlap_mpo_right_step(MPSTensor<Matrix, SymmGroup> const & bra_tensor,
                           MPSTensor<Matrix, SymmGroup> const & ket_tensor,
                           Boundary<Matrix, SymmGroup> const & right,
                           MPOTensor<Matrix, SymmGroup> const & mpo)
    {
        assert(right.aux_dim() == mpo.col_dim());
        
        bra_tensor.make_left_paired();
        ket_tensor.make_left_paired();
        
        block_matrix<Matrix, SymmGroup> t1;
        std::vector<block_matrix<Matrix, SymmGroup> > t2(right.aux_dim());
        
        for (std::size_t b = 0; b < right.aux_dim(); ++b)
        {
            gemm(ket_tensor.data_, right.data_[b], t1);
            gemm(t1, conjugate(transpose(bra_tensor.data_)), t2[b]);
        }
        
        Boundary<Matrix, SymmGroup> ret;
        ret.data_.resize(mpo.row_dim());
        
        typedef typename SymmGroup::charge charge;
        typedef std::size_t size_t;
        
        for (size_t b1 = 0; b1 < mpo.row_dim(); ++b1)
            for (size_t b2 = 0; b2 < mpo.col_dim(); ++b2)
            {
                block_matrix<Matrix, SymmGroup> const & T = t2[b2];
                block_matrix<Matrix, SymmGroup> const & W = mpo(b1, b2);
                
                if (T.n_blocks() == 0 || W.n_blocks() == 0)
                    continue;
                
                ProductBasis<SymmGroup> upper_pb(ket_tensor.site_dim(), ket_tensor.row_dim());
                ProductBasis<SymmGroup> lower_pb(bra_tensor.site_dim(), bra_tensor.row_dim());
                
                for (size_t s1 = 0; s1 < ket_tensor.site_dim().size(); ++s1)
                    for (size_t u = 0; u < ket_tensor.row_dim().size(); ++u)
                    {
                        charge u_charge = ket_tensor.row_dim()[u].first;
                        charge s1_charge = ket_tensor.site_dim()[s1].first;
                        charge tu_charge = SymmGroup::fuse(u_charge, s1_charge);
                        
                        size_t upper_size = ket_tensor.row_dim()[u].second;
                        size_t s1_size = ket_tensor.site_dim()[s1].second;
                        
                        for (size_t s2 = 0; s2 < bra_tensor.site_dim().size(); ++s2)
                            for (size_t l = 0; l < bra_tensor.row_dim().size(); ++l)
                            {   
                                charge l_charge = bra_tensor.row_dim()[l].first;
                                charge s2_charge = bra_tensor.site_dim()[s2].first;
                                charge tl_charge = SymmGroup::fuse(l_charge, s2_charge);
                                
                                size_t lower_size = bra_tensor.row_dim()[l].second;
                                size_t s2_size = bra_tensor.site_dim()[s2].second;
                                
                                Matrix block(upper_size, lower_size);
                                
                                size_t in_l_offset = upper_pb(s1_charge, u_charge);
                                size_t in_r_offset = lower_pb(s2_charge, l_charge);
                                
                                if (! W.has_block(s1_charge, s2_charge) )
                                    continue;
                                
                                if (! T.has_block(tu_charge, tl_charge) )
                                    continue;
                                
                                for (size_t ss1 = 0; ss1 < s1_size; ++ss1)
                                    for (size_t ss2 = 0; ss2 < s2_size; ++ss2)
                                        for (size_t uu = 0; uu < upper_size; ++uu)
                                            for (size_t ll = 0; ll < lower_size; ++ll)
                                                block(uu, ll) += W(s1_charge, s2_charge)(ss1, ss2) *
                                                T(tu_charge, tl_charge)(in_l_offset+ss1*upper_size+uu,
                                                                        in_r_offset+ss2*lower_size+ll);
                                
                                ret.data_[b1] += block_matrix<Matrix, SymmGroup>(u_charge, l_charge, block);
                            }
                    }
            }
        
        return ret;
    }
    
    template<class Matrix, class SymmGroup>
    static MPSTensor<Matrix, SymmGroup>
    site_hamil(MPSTensor<Matrix, SymmGroup> const & ket_tensor,
               Boundary<Matrix, SymmGroup> const & left,
               Boundary<Matrix, SymmGroup> const & right,
               MPOTensor<Matrix, SymmGroup> const & mpo)
    {
        ket_tensor.make_right_paired();
        
        std::vector<block_matrix<Matrix, SymmGroup> > t(left.aux_dim());
        
        for (std::size_t b = 0; b < left.aux_dim(); ++b) {
            gemm(transpose(left.data_[b]), ket_tensor.data_, t[b]);
            block_matrix<Matrix, SymmGroup> tmp;
            reshape_right_to_left<Matrix>(ket_tensor.site_dim(), left.data_[b].right_basis(), ket_tensor.col_dim(),
                                          t[b], tmp);
            swap(t[b], tmp);
        }
        
        Index<SymmGroup> physical_i = ket_tensor.site_dim(), left_i = ket_tensor.row_dim(), right_i = ket_tensor.col_dim();
//        MPSTensor<Matrix, SymmGroup> ret(physical_i, left_i, right_i, false);
//        ret.make_left_paired();
//        ret.data_.clear();
        
        MPSTensor<Matrix, SymmGroup> ret = ket_tensor;
        ret.multiply_by_scalar(0);
        ret.make_left_paired();
        
        typedef typename SymmGroup::charge charge;
        typedef std::size_t size_t;
        
        for (size_t b1 = 0; b1 < left.aux_dim(); ++b1)
            for (size_t b2 = 0; b2 < right.aux_dim(); ++b2)
            {
                block_matrix<Matrix, SymmGroup> T;
                gemm(t[b1], right.data_[b2], T);
                
                block_matrix<Matrix, SymmGroup> const & W = mpo(b1, b2);
                
                if (T.n_blocks() == 0 || W.n_blocks() == 0)
                    continue;
                
                ProductBasis<SymmGroup> out_left_pb(physical_i, left_i);
                ProductBasis<SymmGroup> in_left_pb(physical_i, left.data_[b1].right_basis());
                
                for (size_t s1 = 0; s1 < physical_i.size(); ++s1)
                    for (size_t s2 = 0; s2 < physical_i.size(); ++s2)
                        for (size_t l = 0; l < left_i.size(); ++l)
                            for (size_t r = 0; r < right_i.size(); ++r)
                            {
                                charge T_l_charge = SymmGroup::fuse(physical_i[s1].first, left_i[l].first);
                                charge T_r_charge = right_i[r].first;
                                
                                charge out_l_charge = SymmGroup::fuse(physical_i[s2].first, left_i[l].first);
                                charge out_r_charge = right_i[r].first;
                                
                                if (! T.has_block(T_l_charge, T_r_charge) )
                                    continue;
                                if (! W.has_block(physical_i[s1].first, physical_i[s2].first) )
                                    continue;
                                if (! left.data_[b1].right_basis().has(left_i[l].first) )
                                    continue;
                                if (! right.data_[b2].right_basis().has(right_i[r].first) )
                                    continue;
                                if (out_l_charge != out_r_charge)
                                    continue;
                                
                                size_t in_left_offset = in_left_pb(physical_i[s1].first, left_i[l].first);
                                size_t out_left_offset = out_left_pb(physical_i[s2].first, left_i[l].first);
                                
                                Matrix const & wblock = W(physical_i[s1].first, physical_i[s2].first);
                                Matrix const & iblock = T(T_l_charge, T_r_charge);
                                Matrix oblock(out_left_offset + physical_i[s2].second * left_i[l].second, right_i[r].second);
                                oblock *= 0;

                                /* optimize me */ 
                                for (size_t ss1 = 0; ss1 < physical_i[s1].second; ++ss1)
                                    for (size_t ss2 = 0; ss2 < physical_i[s2].second; ++ss2) {
                                        typename Matrix::value_type wblock_t = wblock(ss1, ss2);
                                        for (size_t rr = 0; rr < right_i[r].second; ++rr) {
                                            typename Matrix::value_type * p1 = &oblock(out_left_offset + ss2*left_i[l].second, rr);
                                            typename Matrix::value_type const * p2 = &iblock(in_left_offset + ss1*left_i[l].second, rr);
                                            for (size_t ll = 0; ll < left_i[l].second; ++ll) {
                                                *(p1++) += *(p2++) * wblock_t;
                                                //oblock(out_left_offset + ss2*left_i[l].second+ll, rr) +=
                                                //iblock(in_left_offset + ss1*left_i[l].second+ll, rr) * wblock(ss1, ss2);
                                            }
                                        }
                                    }
                                
                                ret.data_.match_and_add_block(boost::tuples::make_tuple(oblock, out_l_charge, out_r_charge));
                            }
            }
        
//        ket_tensor.make_left_paired();
//        cout << "ket_tensor: " << ket_tensor << endl;
//        cout << "ret: " << ret << endl;
        
#ifndef NDEBUG
        ket_tensor.make_left_paired();
        assert(ret.data_.left_basis() == ket_tensor.data_.left_basis());
        assert(ret.data_.right_basis() == ket_tensor.data_.right_basis());
#endif
        
//        
//        static int c = 20;
//        if (c-- == 0)
//            exit(0);
        
        return ret;
    }
    
    template<class Matrix, class SymmGroup>
    static Boundary<Matrix, SymmGroup>
    left_boundary_tensor_mpo(MPSTensor<Matrix, SymmGroup> const & mps,
                             Boundary<Matrix, SymmGroup> const & left,
                             MPOTensor<Matrix, SymmGroup> const & mpo)
    {
        mps.make_right_paired();
        
        std::vector<block_matrix<Matrix, SymmGroup> > t(left.aux_dim());
        
        for (std::size_t b = 0; b < left.aux_dim(); ++b) {
            gemm(transpose(left.data_[b]), mps.data_, t[b]);
            block_matrix<Matrix, SymmGroup> tmp;
            reshape_right_to_left<Matrix>(mps.site_dim(), left.data_[b].right_basis(), mps.col_dim(),
                                          t[b], tmp);
            swap(t[b], tmp);
        }
        
        Index<SymmGroup> physical_i = mps.site_dim(), left_i = mps.row_dim(), right_i = mps.col_dim();
        
        Boundary<Matrix, SymmGroup> ret;
        ret.data_.resize(mpo.col_dim());
        
        typedef typename SymmGroup::charge charge;
        typedef std::size_t size_t;
        
        for (size_t b1 = 0; b1 < left.aux_dim(); ++b1)
            for (size_t b2 = 0; b2 < mpo.col_dim(); ++b2)
            {
                block_matrix<Matrix, SymmGroup> const & W = mpo(b1, b2);
                block_matrix<Matrix, SymmGroup> const & T = t[b1];
                
                ProductBasis<SymmGroup> out_left_pb(physical_i, left_i);
                ProductBasis<SymmGroup> in_left_pb(physical_i, left.data_[b1].right_basis());
                
                for (size_t s1 = 0; s1 < physical_i.size(); ++s1)
                    for (size_t s2 = 0; s2 < physical_i.size(); ++s2)
                        for (size_t l = 0; l < left_i.size(); ++l)
                            for (size_t r = 0; r < right_i.size(); ++r)
                            {
                                charge T_l_charge = SymmGroup::fuse(physical_i[s1].first, left_i[l].first);
                                charge T_r_charge = right_i[r].first;
                                
                                charge out_l_charge = SymmGroup::fuse(physical_i[s2].first, left_i[l].first);
                                charge out_r_charge = right_i[r].first;
                                
                                if (! T.has_block(T_l_charge, T_r_charge) )
                                    continue;
                                if (! W.has_block(physical_i[s1].first, physical_i[s2].first) )
                                    continue;
                                if (! left.data_[b1].right_basis().has(left_i[l].first) )
                                    continue;
                                if (! mps.col_dim().has(right_i[r].first) )
                                    continue;
                                
                                size_t in_left_offset = in_left_pb(physical_i[s1].first, left_i[l].first);
                                size_t out_left_offset = out_left_pb(physical_i[s2].first, left_i[l].first);
                                
                                Matrix const & wblock = W(physical_i[s1].first, physical_i[s2].first);
                                Matrix const & iblock = T(T_l_charge, T_r_charge);
                                Matrix oblock(out_left_offset + physical_i[s2].second * left_i[l].second, right_i[r].second);
                                oblock *= 0;
                                
                                /* optimize me */
                                for (size_t ss1 = 0; ss1 < physical_i[s1].second; ++ss1)
                                    for (size_t ss2 = 0; ss2 < physical_i[s2].second; ++ss2) {
                                        typename Matrix::value_type wblock_t = wblock(ss1, ss2);
                                        for (size_t rr = 0; rr < right_i[r].second; ++rr) {
                                            for (size_t ll = 0; ll < left_i[l].second; ++ll) {
                                                oblock(out_left_offset + ss2*left_i[l].second+ll, rr) +=
                                                iblock(in_left_offset + ss1*left_i[l].second+ll, rr) * wblock_t;
                                            }
                                        }
                                    }
                                
                                ret.data_[b2].match_and_add_block(boost::tuples::make_tuple(oblock, out_l_charge, out_r_charge));
                            }
            }
        
        return ret;
    }
    
    template<class Matrix, class SymmGroup>
    static MPSTensor<Matrix, SymmGroup>
    site_hamil2(MPSTensor<Matrix, SymmGroup> const & ket_tensor,
                Boundary<Matrix, SymmGroup> const & left,
                Boundary<Matrix, SymmGroup> const & right,
                MPOTensor<Matrix, SymmGroup> const & mpo)
    {
        Boundary<Matrix, SymmGroup> left_mpo_mps = left_boundary_tensor_mpo(ket_tensor, left, mpo);
        
        MPSTensor<Matrix, SymmGroup> ret = ket_tensor;
        ret.multiply_by_scalar(0);
        ret.make_left_paired();
        
        typedef typename SymmGroup::charge charge;
        typedef std::size_t size_t;
        
        for (size_t b = 0; b < mpo.col_dim(); ++b)
        {
            block_matrix<Matrix, SymmGroup> oblock;
            gemm(left_mpo_mps.data_[b], right.data_[b], oblock);
            for (size_t k = 0; k < oblock.n_blocks(); ++k)
                ret.data_.match_and_add_block(boost::tuples::make_tuple(oblock[k],
                                                                        oblock.left_basis()[k].first,
                                                                        oblock.right_basis()[k].first));
            
        }
        
        return ret;
    }
    
    template<class Matrix, class SymmGroup>
    static MPSTensor<Matrix, SymmGroup>
    predict_new_state_l2r_sweep(MPSTensor<Matrix, SymmGroup> const & mps,
                                MPOTensor<Matrix, SymmGroup> const & mpo,
                                Boundary<Matrix, SymmGroup> const & left,
                                Boundary<Matrix, SymmGroup> const & right,
                                double alpha, double cutoff, std::size_t Mmax)
    {
        mps.make_left_paired();
        block_matrix<Matrix, SymmGroup> dm;
        gemm(mps.data_, transpose(mps.data_), dm);
        
        Boundary<Matrix, SymmGroup> half_dm = left_boundary_tensor_mpo(mps, left, mpo);
        
        mps.make_left_paired();
        for (std::size_t b = 0; b < half_dm.aux_dim(); ++b)
        {
            block_matrix<Matrix, SymmGroup> tdm;
            gemm(half_dm.data_[b], transpose(half_dm.data_[b]), tdm);
            
//            cout << "DM: " << dm.left_basis() << " " << dm.right_basis() << endl;
//            cout << "tDM: " << tdm.left_basis() << " " << tdm.right_basis() << endl;
            
            tdm *= alpha;
            for (std::size_t k = 0; k < tdm.n_blocks(); ++k) {
                if (mps.data_.left_basis().has(tdm.left_basis()[k].first))
                    dm.match_and_add_block(boost::tuples::make_tuple(tdm[k],
                                                                     tdm.left_basis()[k].first,
                                                                     tdm.right_basis()[k].first));
            }
        }
        
        mps.make_left_paired();
        assert(dm.left_basis() == mps.data_.left_basis());
        
        block_matrix<Matrix, SymmGroup> U, V;
        block_matrix<typename blas::associated_diagonal_matrix<Matrix>::type, SymmGroup> S, sqrtS;
        
        svd(dm, U, V, S, cutoff, Mmax);
        
        MPSTensor<Matrix, SymmGroup> ret = mps;
        ret.data_ = U;
        ret.right_i = U.right_basis();
        return ret;
    }
    
    template<class Matrix, class SymmGroup>
    static MPSTensor<Matrix, SymmGroup>
    predict_lanczos_l2r_sweep(MPSTensor<Matrix, SymmGroup> B,
                              MPSTensor<Matrix, SymmGroup> const & psi,
                              MPSTensor<Matrix, SymmGroup> const & A)
    {
        psi.make_left_paired();
        A.make_left_paired();
        
        block_matrix<Matrix, SymmGroup> tmp;
        gemm(transpose(A.data_), psi.data_, tmp);
        
        B.multiply_from_left(tmp);
        return B;
    }
};

#endif
