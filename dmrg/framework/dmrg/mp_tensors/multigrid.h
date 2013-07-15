/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2012 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MULTIGRID_H
#define MULTIGRID_H

#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/optimize.h"

#include <boost/optional.hpp>

struct multigrid {
    
    /******** Isometry T *******
     * Elements stored as:     *
     * T(s, fuse(s1,s2))       *
     *                         *
     *          s1  s2         *
     *          ^   ^          *
     *          |   |          *
     *          ^   ^          *
     *          -----          *
     *          \ T /          *
     *           \ /           *
     *            ^            *
     *            |            *
     *            ^            *
     *            s            *
     ***************************/

    template<class Matrix, class SymmGroup>
    void
    static normalize_restriction_isometry(block_matrix<Matrix, SymmGroup> & T)
    {
        for (size_t k=0; k<T.n_blocks(); k++)
            for (size_t i=0; i<num_rows(T[k]); ++i)
            {
                typename Matrix::value_type nn = 0.;
                for (size_t j=0; j<num_cols(T[k]); ++j) nn += T[k](i,j);
                nn = sqrt(nn);
                for (size_t j=0; j<num_cols(T[k]); ++j) T[k](i,j) /= nn;
            }
    }

    template<class Matrix, class SymmGroup>
    block_matrix<Matrix, SymmGroup>
    static default_restriction_isometry(Index<SymmGroup> const & phys_large,
                                        Index<SymmGroup> const & phys_small)
    {
        ProductBasis<SymmGroup> pb(phys_large, phys_large);
        
        typedef std::size_t size_t;
        typedef typename SymmGroup::charge charge;
        block_matrix<Matrix, SymmGroup> T;
        for (size_t s1=0; s1<phys_large.size(); ++s1)
            for (size_t s2=0; s2<phys_large.size(); ++s2)
            {
                charge s1_charge = phys_large[s1].first;
                charge s2_charge = phys_large[s2].first;
                charge s_charge = SymmGroup::fuse(s1_charge, s2_charge);
                
                size_t s = phys_small.position(s_charge);
                if (s == phys_small.size())
                    continue;
                
                if (! T.has_block(s_charge, s_charge))
                    T.insert_block(new Matrix(phys_small[s].second, pb.size(s1_charge, s2_charge), 0.), s_charge, s_charge);
                
                size_t out_offset = pb(s1_charge, s2_charge);
                Matrix & oblock = T(s_charge, s_charge);
                for(size_t ss1 = 0; ss1 < phys_large[s1].second; ++ss1)
                    for(size_t ss2 = 0; ss2 < phys_large[s2].second; ++ss2)
                        oblock(0, out_offset + ss1*phys_large[s2].second+ss2) = typename Matrix::value_type(1.);
            }
        normalize_restriction_isometry(T);
        return T;
    }
    
    template<class Matrix, class SymmGroup>
    MPSTensor<Matrix, SymmGroup>
    static restriction(MPSTensor<Matrix, SymmGroup> const & M1, MPSTensor<Matrix, SymmGroup> const & M2,
                       block_matrix<Matrix, SymmGroup> const & T)
    {
        M1.make_left_paired();
        M2.make_right_paired();
        
        // contracting M1, M2
        block_matrix<Matrix, SymmGroup> MM;
        gemm(M1.data(), M2.data(), MM);
        
        // init index objects
        typedef typename SymmGroup::charge charge;
        typedef std::size_t size_t;
        
        Index<SymmGroup> left_i = M1.row_dim(), right_i = M2.col_dim();
        Index<SymmGroup> phys_large = M1.site_dim(), phys_small = T.left_basis();
        
        ProductBasis<SymmGroup> in_left(phys_large, left_i);
        ProductBasis<SymmGroup> in_right(phys_large, right_i,
                                         boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                                                 -boost::lambda::_1, boost::lambda::_2));
        ProductBasis<SymmGroup> out_left(phys_small, left_i);
        ProductBasis<SymmGroup> phys_pb(phys_large, phys_large);
        
        // reshape + multiplication
        block_matrix<Matrix, SymmGroup> data;
        for (size_t Mblock = 0; Mblock < MM.n_blocks(); ++Mblock)
            for (size_t Tb = 0; Tb < T.n_blocks(); ++Tb)
                for (size_t s1=0; s1<phys_large.size(); ++s1) {
                    size_t s2 = phys_large.position(SymmGroup::fuse(T.right_basis()[Tb].first,
                                                                    -phys_large[s1].first));
                    if(s2 == phys_large.size()) continue;
                    size_t l = left_i.position(SymmGroup::fuse(MM.left_basis()[Mblock].first,
                                                               -phys_large[s1].first));
                    if(l == left_i.size()) continue;
                    size_t r = right_i.position(SymmGroup::fuse(MM.right_basis()[Mblock].first,
                                                                phys_large[s2].first));
                    if(r == right_i.size()) continue;
                    size_t s = phys_small.position(T.left_basis()[Tb].first);
                    
                    {
                        charge in_l_charge = SymmGroup::fuse(phys_large[s1].first, left_i[l].first);
                        charge in_r_charge = SymmGroup::fuse(-phys_large[s2].first, right_i[r].first);
                        
                        charge T_l_charge = T.left_basis()[Tb].first;
                        charge T_r_charge = SymmGroup::fuse(phys_large[s1].first, phys_large[s2].first);
                        
                        charge out_l_charge = SymmGroup::fuse(phys_small[s].first, left_i[l].first);
                        charge out_r_charge = right_i[r].first;
                        
                        if (! MM.has_block(in_l_charge, in_r_charge) )
                            continue;
                        if (! T.has_block(T_l_charge, T_r_charge) )
                            continue;
                        
                        if (! data.has_block(out_l_charge, out_r_charge))
                            data.insert_block(new Matrix(out_left.size(phys_small[s].first, left_i[l].first), right_i[r].second, 0),
                                              out_l_charge, out_r_charge);
                        
                        size_t in_l_offset = in_left(phys_large[s1].first, left_i[l].first);
                        size_t in_r_offset = in_right(phys_large[s2].first, right_i[r].first);
                        size_t out_l_offset = out_left(phys_small[s].first, left_i[l].first);
                        size_t T_r_offset = phys_pb(phys_large[s1].first, phys_large[s2].first);
                        
                        Matrix const & in_block = MM(in_l_charge, in_r_charge);
                        Matrix const & T_block = T(T_l_charge, T_r_charge);
                        Matrix & out_block = data(out_l_charge, out_r_charge);
                        
                        for (size_t ss=0; ss<phys_small[s].second; ++ss) {
                            for(size_t ss1 = 0; ss1 < phys_large[s1].second; ++ss1)
                                for(size_t ss2 = 0; ss2 < phys_large[s2].second; ++ss2) {
                                    typename Matrix::value_type const & Tval = T_block(ss, T_r_offset + ss1*phys_large[s2].second+ss2);
                                    for(size_t rr = 0; rr < right_i[r].second; ++rr)
                                    /* blas implementation */
                                        maquis::dmrg::detail::iterator_axpy(&in_block(in_l_offset + ss1*left_i[l].second, in_r_offset + ss2*right_i[r].second+rr),
                                                      &in_block(in_l_offset + ss1*left_i[l].second, in_r_offset + ss2*right_i[r].second+rr) + left_i[l].second,
                                                      &out_block(out_l_offset + ss*left_i[l].second, rr),
                                                      Tval);
                                    /* loop implementation:
                                     for(size_t ll = 0; ll < left_i[l].second; ++ll)
                                     out_block(out_l_offset + ss*left_i[l].second+ll, rr)
                                     += in_block(in_l_offset + ss1*left_i[l].second+ll, in_r_offset + ss2*right_i[r].second+rr) * Tval;
                                     */
                                }
                        }
                    }
                    
                }
        
        MPSTensor<Matrix, SymmGroup> M;
        swap(M.data_, data);
        M.cur_storage = LeftPaired;
        M.left_i = left_i;
        M.right_i = right_i;
        M.phys_i = phys_small;
        return M;
    }
    
    template<class Matrix, class SymmGroup>
    MPS<Matrix, SymmGroup>
    static restrict_mps(MPS<Matrix, SymmGroup> const & mps_large,
                        boost::optional<block_matrix<Matrix, SymmGroup> const& > T_ = (boost::optional<block_matrix<Matrix, SymmGroup> const& >()))
    {
        std::size_t LL = mps_large.length();
        assert(LL % 2 == 0);
        std::size_t L = LL/2;
        
        MPS<Matrix, SymmGroup> mps_small; mps_small.resize(L);

        block_matrix<Matrix, SymmGroup> T;
        if (!T_)
            // Assuming same physical index everywhere!
            T = default_restriction_isometry<Matrix>(mps_large.site_dim(0), mps_large.site_dim(0));
        else
            T = T_.get();
        
        for (std::size_t p = 0; p < L; ++p)
        {            
            mps_small[p] = restriction(mps_large[2*p], mps_large[2*p+1], T);
            mps_small.move_normalization_l2r(p, p+1);
        }
        return mps_small;
    }

    template<class Matrix, class SymmGroup>
    MPO<Matrix, SymmGroup>
    static restrict_mpo(MPO<Matrix, SymmGroup> const & mpo_large,
                        boost::optional<block_matrix<Matrix, SymmGroup> const& > T_ = (boost::optional<block_matrix<Matrix, SymmGroup> const& >()))
    {
        std::size_t LL = mpo_large.length();
        assert(LL % 2 == 0);
        std::size_t L = LL/2;
        
        MPO<Matrix, SymmGroup> mpo_small(L);
        
        block_matrix<Matrix, SymmGroup> T;
        if (!T_) {
            // Assuming same physical index everywhere!
            Index<SymmGroup> phys_i;
            for (size_t b1=0; b1<mpo_large[0].col_dim() && phys_i.size()==0; ++b1)
                for (size_t b2=0; b2<mpo_large[0].row_dim() && phys_i.size()==0; ++b2)
                    if (mpo_large[0].has(b1, b2))
                        phys_i = mpo_large[0](b1,b2).left_basis();

            T = default_restriction_isometry<Matrix>(phys_i, phys_i);
        } else {
            T = T_.get();
        }
        
        for (std::size_t p = 0; p < L; ++p)
        {
            Index<SymmGroup> phys_i;
            for (size_t b1=0; b1<mpo_large[2*p].col_dim() && phys_i.size()==0; ++b1)
                for (size_t b2=0; b2<mpo_large[2*p].row_dim() && phys_i.size()==0; ++b2)
                    if (mpo_large[2*p].has(b1, b2))
                        phys_i = mpo_large[2*p](b1,b2).left_basis();
            MPOTensor<Matrix, SymmGroup> mpot = make_twosite_mpo(mpo_large[2*p], mpo_large[2*p+1], phys_i);
            for (size_t b1=0; b1<mpot.col_dim(); ++b1)
                for (size_t b2=0; b2<mpot.row_dim(); ++b2)
                {
                    if (! mpot.has(b1, b2))
                        continue;
                    block_matrix<Matrix, SymmGroup> tmp;
                    gemm(T, mpot(b1, b2), tmp);
                    gemm(tmp, transpose(T), mpot(b1, b2));
                }
            mpo_small[p] = mpot;
        }
        return mpo_small;
    }

    template<class Matrix, class SymmGroup>
    void
    static restriction_old (MPS<Matrix, SymmGroup> const & mps_large,
                            MPS<Matrix, SymmGroup> & mps_small)
    {
        std::size_t L = mps_small.length();
        std::size_t LL = mps_large.length();
        assert(LL == 2*L);
        
        typedef typename SymmGroup::charge charge;
        typedef typename Index<SymmGroup>::basis_iterator bi_t;
        
        
        for (std::size_t p = 0; p < L; ++p)
        {
            //            maquis::cout << std::endl << std::endl << "********************" << std::endl << "starting loop for p = " << p << std::endl;
            
            block_matrix<Matrix, SymmGroup> M;
            
            Index<SymmGroup> alpha_basis = mps_large.row_dim(2*p);
            Index<SymmGroup> b_basis = mps_large.col_dim(2*p);
            Index<SymmGroup> beta_basis = mps_large.col_dim(2*p+1);
            Index<SymmGroup> s_basis = mps_small.site_dim(p);
            Index<SymmGroup> s1_basis = mps_large.site_dim(2*p);
            Index<SymmGroup> s2_basis = mps_large.site_dim(2*p+1);
            
            //            maquis::cout << "alpha_basis:" <<std::endl << alpha_basis << std::endl;
            //            maquis::cout << "beta_basis:" <<std::endl << beta_basis << std::endl;
            //            maquis::cout << "s_basis:" <<std::endl << s_basis << std::endl;
            //            maquis::cout << "s1_basis:" <<std::endl << s1_basis << std::endl;
            //            maquis::cout << "s2_basis:" <<std::endl << s2_basis << std::endl;
            
            
            mps_large[2*p].make_left_paired();
            mps_large[2*p+1].make_right_paired();
            
            ProductBasis<SymmGroup> in1_left(s1_basis, alpha_basis);
            ProductBasis<SymmGroup> in2_right(s2_basis, beta_basis,
                                              boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                                  -boost::lambda::_1, boost::lambda::_2));
            
            ProductBasis<SymmGroup> out_left(s_basis, alpha_basis);
            
            for (bi_t alpha = alpha_basis.basis_begin(); !alpha.end(); ++alpha)
                for (bi_t beta = beta_basis.basis_begin(); !beta.end(); ++beta)
                    for (bi_t b = b_basis.basis_begin(); !b.end(); ++b)
                        for (bi_t s1 = s1_basis.basis_begin(); !s1.end(); ++s1)
                            for (bi_t s2 = s2_basis.basis_begin(); !s2.end(); ++s2)
                            {
                                std::pair<charge, std::size_t> s(SymmGroup::fuse(s1->first, s2->first),
                                                                 s1->second*s2_basis.size_of_block(s2->first)+s2->second);
                                
                                if (!s_basis.has(s.first) || s_basis.size_of_block(s.first) <= s.second)
                                    continue;
                                
                                charge in1_left_c = SymmGroup::fuse(s1->first, alpha->first);
                                charge in1_right_c = b->first;
                                
                                charge in2_left_c = b->first;
                                charge in2_right_c = SymmGroup::fuse(-s2->first, beta->first);
                                
                                if (!mps_large[2*p].data().has_block(in1_left_c, in1_right_c))
                                    continue;
                                if (!mps_large[2*p+1].data().has_block(in2_left_c, in2_right_c))
                                    continue;
                                
                                charge out_left_c = SymmGroup::fuse(s.first, alpha->first);
                                charge out_right_c = beta->first;
                                
                                //                            maquis::cout << "--" << std::endl;
                                //                            maquis::cout << "alpha: " << alpha->first << ", " << alpha->second << std::endl;
                                //                            maquis::cout << "beta: " << beta->first << ", " << beta->second << std::endl;
                                //                            maquis::cout << "s1: " << s1->first << ", " << s1->second << std::endl;
                                //                            maquis::cout << "s2: " << s2->first << ", " << s2->second << std::endl;
                                //                            maquis::cout << "s: " << s.first << ", " << s.second << std::endl;
                                
                                if (!M.has_block(out_left_c, out_right_c))
                                    M.insert_block(new Matrix(out_left.size(s.first, alpha->first),
                                                              beta_basis.size_of_block(beta->first),
                                                              0),
                                                   out_left_c, out_right_c);
                                
                                
                                //                            maquis::cout << "block has size " << out_left.size(s1->first, alpha->first)
                                //                            << "x"
                                //                            << out_right.size(s2->first, beta->first,
                                //                                              boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                //                                                                  -boost::lambda::_1, boost::lambda::_2)) << std::endl;
                                
                                std::size_t in1_left_offset = in1_left(s1->first, alpha->first);
                                std::size_t in1_right_offset = 0;
                                std::size_t in2_left_offset = 0;
                                std::size_t in2_right_offset = in2_right(s2->first, beta->first);
                                
                                std::size_t out_left_offset = out_left(s.first, alpha->first);
                                std::size_t out_right_offset = 0;
                                
                                //                            maquis::cout << "M[" << out_left_c << ", " << out_right_c << "]"
                                //                            << "(" << out_left_offset + s1->second*alpha_basis.size_of_block(alpha->first)+alpha->second << ", " << out_right_offset + s2->second*beta_basis.size_of_block(beta->first)+beta->second << ")" << std::endl;
                                //                            maquis::cout << " = " << M(std::make_pair(out_left_c, out_left_offset + s1->second*alpha_basis.size_of_block(alpha->first)+alpha->second),
                                //                                                    std::make_pair(out_right_c, out_right_offset + s2->second*beta_basis.size_of_block(beta->first)+beta->second)) << std::endl;
                                //                            
                                //                            maquis::cout << "mps[" << in_left_c << ", " << in_right_c << "]"
                                //                            << "(" << in_left_offset + s.second*alpha_basis.size_of_block(alpha->first) + alpha->second << ", " << in_right_offset + beta->second << ")" << std::endl;
                                //                            maquis::cout << " = " << mps_small[p].data()(std::make_pair(in_left_c, in_left_offset + s.second*alpha_basis.size_of_block(alpha->first) + alpha->second),
                                //                                                                      std::make_pair(in_right_c, in_right_offset + beta->second)) << std::endl;
                                
                                M(std::make_pair(out_left_c, out_left_offset + s.second*alpha_basis.size_of_block(alpha->first) + alpha->second),
                                  std::make_pair(out_right_c, out_right_offset + beta->second))
                                += mps_large[2*p].data()(std::make_pair(in1_left_c, in1_left_offset + s1->second*alpha_basis.size_of_block(alpha->first) + alpha->second),
                                                         std::make_pair(in1_right_c, in1_right_offset + b->second))
                                * mps_large[2*p+1].data()(std::make_pair(in2_left_c, in2_left_offset + b->second),
                                                          std::make_pair(in2_right_c, in2_right_offset + s2->second*beta_basis.size_of_block(beta->first) + beta->second))
                                ;
                                
                            }
            
            
            mps_small[p].data() = M;
            mps_small[p].left_i = alpha_basis;
            mps_small[p].right_i = beta_basis;
            mps_small[p].cur_storage = LeftPaired;
            
        }
        
        mps_small.normalize_left();
        
    }
    
    
    template<class Matrix, class SymmGroup>
    void
    static site_extension (MPSTensor<Matrix, SymmGroup> const & mps_small,
                           MPSTensor<Matrix, SymmGroup> & M1,
                           MPSTensor<Matrix, SymmGroup> & M2)
    {        
        typedef typename SymmGroup::charge charge;
        typedef typename Index<SymmGroup>::basis_iterator bi_t;
        
        block_matrix<Matrix, SymmGroup> M;
        
        Index<SymmGroup> alpha_basis = mps_small.row_dim();
        Index<SymmGroup> beta_basis = mps_small.col_dim();
        Index<SymmGroup> s_basis = mps_small.site_dim();
        Index<SymmGroup> s1_basis = M1.site_dim();
        Index<SymmGroup> s2_basis = M2.site_dim();
        
        
        mps_small.make_left_paired();
        
        ProductBasis<SymmGroup> in_left(s_basis, alpha_basis);
        
        ProductBasis<SymmGroup> out_left(s1_basis, alpha_basis);
        ProductBasis<SymmGroup> out_right(s2_basis, beta_basis,
                                          boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                              -boost::lambda::_1, boost::lambda::_2));
        
        for (bi_t alpha = alpha_basis.basis_begin(); !alpha.end(); ++alpha)
            for (bi_t beta = beta_basis.basis_begin(); !beta.end(); ++beta)
                for (bi_t s1 = s1_basis.basis_begin(); !s1.end(); ++s1)
                    for (bi_t s2 = s2_basis.basis_begin(); !s2.end(); ++s2)
                    {
                        std::pair<charge, std::size_t> s(SymmGroup::fuse(s1->first, s2->first),
                                                         s1->second*s2_basis.size_of_block(s2->first)+s2->second);
                        
                        if (!s_basis.has(s.first) || s_basis.size_of_block(s.first) <= s.second)
                            continue;
                        
                        charge in_left_c = SymmGroup::fuse(s.first, alpha->first);
                        charge in_right_c = beta->first;
                        
                        if (!mps_small.data().has_block(in_left_c, in_right_c))
                            continue;
                        
                        charge out_left_c = SymmGroup::fuse(s1->first, alpha->first);
                        charge out_right_c = SymmGroup::fuse(-s2->first, beta->first);
                        
                        if (!M.has_block(out_left_c, out_right_c))
                            M.insert_block(new Matrix(out_left.size(s1->first, alpha->first),
                                                      out_right.size(s2->first, beta->first,
                                                                    boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                                                         -boost::lambda::_1, boost::lambda::_2)),
                                                      0),
                                           out_left_c, out_right_c);
                        
                        std::size_t in_left_offset = in_left(s.first, alpha->first);
                        std::size_t in_right_offset = 0;
                        
                        std::size_t out_left_offset = out_left(s1->first, alpha->first);
                        std::size_t out_right_offset = out_right(s2->first, beta->first);
                        
                        M(std::make_pair(out_left_c, out_left_offset + s1->second*alpha_basis.size_of_block(alpha->first)+alpha->second),
                          std::make_pair(out_right_c, out_right_offset + s2->second*beta_basis.size_of_block(beta->first)+beta->second))
                        = mps_small.data()(std::make_pair(in_left_c, in_left_offset + s.second*alpha_basis.size_of_block(alpha->first) + alpha->second),
                                           std::make_pair(in_right_c, in_right_offset + beta->second));
                        
                    }
        
        
        block_matrix<Matrix, SymmGroup> V, left, right;
        block_matrix<typename alps::numeric::associated_diagonal_matrix<Matrix>::type, SymmGroup> S;
        
        svd(M, left, V, S);
        gemm(S, V, right);
        
        M1.data() = left; // replace_left_paired cannot be used here
        M1.left_i = alpha_basis;
        M1.right_i = left.right_basis();
        M1.cur_storage = LeftPaired;
        
        M2.data() = right; // replace_right_paired cannot be used here
        M2.right_i = beta_basis;
        M2.left_i = right.left_basis();
        M2.cur_storage = RightPaired;
        
    }

    
    template<class Matrix, class SymmGroup>
    void
    static extension (MPS<Matrix, SymmGroup> const & mps_small,
                      MPS<Matrix, SymmGroup> & mps_large)
    {
        std::size_t L1 = mps_small.length();
        std::size_t L2 = mps_large.length();
        assert(L2 == 2*L1);
                
        block_matrix<Matrix, SymmGroup> t_norm;
        
        for (std::size_t p = 0; p < L1; ++p)
        {            
            site_extension(mps_small[p], mps_large[2*p], mps_large[2*p+1]);
                        
            // Normalization
            if (p != 0)
                mps_large[2*p].multiply_from_left(t_norm);
            t_norm = mps_large[2*p].normalize_left(DefaultSolver());
            mps_large[2*p+1].multiply_from_left(t_norm);
            t_norm = mps_large[2*p+1].normalize_left(DefaultSolver());
        }
        
    }
    
    
    template<class Matrix, class SymmGroup>
    void
    static extension_optim (BaseParameters & parms,
                            MPS<Matrix, SymmGroup> mps_small,
                            MPS<Matrix, SymmGroup> & mps_large,
                            std::vector<MPO<Matrix, SymmGroup> > const & mpos_mix)
    {
                
        std::size_t L = mps_small.length();
        std::size_t LL = mps_large.length();
        assert(LL == 2*L);
        
        int sweep = 0;
        
        /*
         * INIT SWEEP PARAMETERS
         */
        double alpha;
        int ngs = parms.get<int>("ngrowsweeps"), nms = parms.get<int>("nmainsweeps");
        if (sweep < ngs)
            alpha = parms.get<double>("alpha_initial");
        else if (sweep < ngs + nms)
            alpha = parms.get<double>("alpha_main");
        else
            alpha = parms.get<double>("alpha_final");
        
        
        double cutoff;
        if (sweep >= parms.get<int>("ngrowsweeps"))
            cutoff = parms.get<double>("truncation_final");
        else
            cutoff = log_interpolate(parms.get<double>("truncation_initial"), parms.get<double>("truncation_final"), parms.get<int>("ngrowsweeps"), sweep);
        
        std::size_t Mmax;
        if (parms.is_set("sweep_bond_dimensions")) {
            std::vector<std::size_t> ssizes = parms.get<std::vector<std::size_t> >("sweep_bond_dimensions");
            if (sweep >= ssizes.size())
                Mmax = *ssizes.rbegin();
            else
                Mmax = ssizes[sweep];
        } else
            Mmax = parms.get<std::size_t>("max_bond_dimension");
        
        
        mps_small.normalize_right();
        /*
         * INIT BOUNDARIES
         */
        
        std::vector<Boundary<Matrix, SymmGroup> > left_(LL+1), right_(L+1);        
        {
            left_[0] = mps_small.left_boundary();
            right_[L] = mps_small.right_boundary();
            
            for(int i = L-1; i >= 0; --i) {
                right_[i] = contraction::overlap_mpo_right_step(mps_small[i], mps_small[i], right_[i+1], mpos_mix[0][i]);
            }
        }
        
        block_matrix<Matrix, SymmGroup> t_norm;
        
        MPSTensor<Matrix, SymmGroup> Msmall = mps_small[0];
        
        typedef typename SymmGroup::charge charge;
        typedef typename Index<SymmGroup>::basis_iterator bi_t;
        
        for (std::size_t p = 0; p < L; ++p)
        {
            
            site_extension(Msmall, mps_large[2*p], mps_large[2*p+1]);
            
            assert( mps_large[2*p].num_check() );
            assert( mps_large[2*p+1].num_check() );
            
            // Normalizing the MPS
            mps_large[2*p+1].divide_by_scalar(mps_large[2*p+1].scalar_norm()); // that may playout
            assert( mps_large[2*p].num_check() );
            assert( mps_large[2*p+1].num_check() );
            t_norm = mps_large[2*p+1].normalize_right(DefaultSolver());
            mps_large[2*p].multiply_from_right(t_norm);
            mps_large[2*p].divide_by_scalar(mps_large[2*p].scalar_norm());
            assert( mps_large[2*p].num_check() );
            assert( mps_large[2*p+1].num_check() );

            
            Boundary<Matrix, SymmGroup> right;
            
            
            Boundary<Matrix, SymmGroup> right_mixed;
            if (p<L-1) {
                right_mixed = right_[p+2];
                right_mixed = contraction::overlap_mpo_right_step(mps_small[p+1], mps_small[p+1], right_mixed, mpos_mix[p+1][2*(p+1)]);
            }

            // Testing energy calculations
            if (false) {
                
                MPS<Matrix, SymmGroup> mps_mixed = mps_small;
                mps_mixed.resize(2*(p+1) + (L-1-p));
                
                std::copy(mps_large.const_begin(), mps_large.const_begin()+2*(p+1), mps_mixed.begin());
                std::copy(mps_small.const_begin()+(p+1), mps_small.const_end(), mps_mixed.begin()+2*(p+1));
                
                
                if (false && p == 0) {
                    for (int shift=0; shift<4; shift++) {
                        maquis::cout << "MPO(" << 2*p+shift << ")::" << std::endl;
                        for (int i=0; i<mpos_mix[p+1][2*p+shift].row_dim(); ++i)
                            for(int j=0; j<mpos_mix[p+1][2*p+shift].col_dim(); ++j)
                                maquis::cout << "(" << i << " --> " << j << "):" << std::endl << mpos_mix[p+1][2*p+shift](i, j);
                    }
                    
                    maquis::cout << "MPO_large(last)::" << std::endl;
                    for (int i=0; i<mpos_mix[L][LL-1].row_dim(); ++i)
                        for(int j=0; j<mpos_mix[L][LL-1].col_dim(); ++j)
                            maquis::cout << "(" << i << " --> " << j << "):" << std::endl << mpos_mix[L][LL-1](i, j);
                    exit(-1);
                    
                }
                
                {
                    maquis::cout << "Norm " << norm(mps_mixed) << std::endl;
                    maquis::cout << "Energy " << "finegraining_fullmpomix " << maquis::real(expval(mps_mixed, mpos_mix[p+1])) << std::endl;
                }
                
            }
            
            /*
             * OPTIMIZATION 2*p
             */
            {
                right = (p<L-1) ? right_mixed : right_[p+1];
                right = contraction::overlap_mpo_right_step(mps_large[2*p+1], mps_large[2*p+1], right, mpos_mix[p+1][2*p+1]);
                if (p == 0)
                    left_[0] = mps_large.left_boundary();
                
                SiteProblem<Matrix, SymmGroup> sp(left_[2*p], right, mpos_mix[p+1][2*p]);

                
                // solver
                if (parms.get<bool>("finegrain_optim"))
                {
                    
                    std::pair<double, MPSTensor<Matrix, SymmGroup> > res;
                    timeval now, then;
                    if (parms["eigensolver"] == std::string("IETL")) {
                        BEGIN_TIMING("IETL")
                        res = solve_ietl_lanczos(sp, mps_large[2*p], parms);
                        END_TIMING("IETL")
                    } else if (parms["eigensolver"] == std::string("IETL_JCD")) {
                        BEGIN_TIMING("JCD")
                        res = solve_ietl_jcd(sp, mps_large[2*p], parms);
                        END_TIMING("JCD")
                    } else {
                        throw std::runtime_error("I don't know this eigensolver.");
                    }
                    
                    mps_large[2*p] = res.second;
                    
                    maquis::cout << "Energy " << "finegraining_1 " << res.first << std::endl;
                    storage::log << std::make_pair("Energy", res.first);

                } else if (true) {
                    // Compute Energy
                    MPSTensor<Matrix, SymmGroup> vec2 =
                    contraction::site_hamil2(mps_large[2*p], sp.left, sp.right, sp.mpo);
                    double energy = mps_large[2*p].scalar_overlap(vec2);
                    maquis::cout << "Energy " << "finegraining_00 " << energy << std::endl;
                    storage::log << std::make_pair("Energy", energy);

                }
                
                // growing
                /*
                
                maquis::cout << "Growing, alpha = " << alpha << std::endl;
                mps_large.grow_l2r_sweep(mpos_mix[L][2*p], left_[2*p], right,
                                         2*p, alpha, cutoff, Mmax);
                 */
                 
                t_norm = mps_large[2*p].normalize_left(DefaultSolver());
                mps_large[2*p+1].multiply_from_left(t_norm);

            }
            
            /*
             * OPTIMIZATION 2*p+1
             */
            {
                right = (p<L-1) ? right_mixed : right_[p+1];
                left_[2*p+1] = contraction::overlap_mpo_left_step(mps_large[2*p], mps_large[2*p],
                                                                  left_[2*p], mpos_mix[L][2*p]);
                
                SiteProblem<Matrix, SymmGroup> sp(left_[2*p+1], right, mpos_mix[p+1][2*p+1]);

                // solver
                if (parms.get<bool>("finegrain_optim"))
                {
                    
                    std::pair<double, MPSTensor<Matrix, SymmGroup> > res;
                    timeval now, then;
                    if (parms["eigensolver"] == std::string("IETL")) {
                        BEGIN_TIMING("IETL")
                        res = solve_ietl_lanczos(sp, mps_large[2*p+1], parms);
                        END_TIMING("IETL")
                    } else if (parms["eigensolver"] == std::string("IETL_JCD")) {
                        BEGIN_TIMING("JCD")
                        res = solve_ietl_jcd(sp, mps_large[2*p+1], parms);
                        END_TIMING("JCD")
                    } else {
                        throw std::runtime_error("I don't know this eigensolver.");
                    }
                    
                    mps_large[2*p+1] = res.second;
                    
                    maquis::cout << "Energy " << "finegraining_2 " << res.first << std::endl;
                    storage::log << std::make_pair("Energy", res.first);

                } else if (true) {
                    // Compute Energy
                    MPSTensor<Matrix, SymmGroup> vec2 =
                    contraction::site_hamil2(mps_large[2*p+1], sp.left, sp.right, sp.mpo);
                    double energy = mps_large[2*p+1].scalar_overlap(vec2);
                    maquis::cout << "Energy " << "finegraining_01 " << energy << std::endl;
                    storage::log << std::make_pair("Energy", energy);
                }
                
                // growing
                /*
                
                if (p < L-1) {
                    maquis::cout << "Growing, alpha = " << alpha << std::endl;
                    MPSTensor<Matrix, SymmGroup> new_mps =
                    contraction::predict_new_state_l2r_sweep(mps_large[2*p+1], mpos_mix[L][2*p+1], left_[2*p+1], right, alpha, cutoff, Mmax);
                    // New tensor for next iteration
                    Msmall = contraction::predict_lanczos_l2r_sweep(mps_small[p+1],
                                                                    mps_large[2*p+1], new_mps);
                    mps_large[2*p+1] = new_mps;
                } else {
                    block_matrix<Matrix, SymmGroup> t = mps_large[2*p+1].normalize_left(DefaultSolver());
                }

                */

                
                t_norm = mps_large[2*p+1].normalize_left(DefaultSolver());

                if (p < L-1) {
                    Msmall = mps_small[p+1];
                    Msmall.multiply_from_left(t_norm);
                }
            }
            
            // Preparing left boundary
            if (p < L-1) {
                left_[2*p+2] = contraction::overlap_mpo_left_step(mps_large[2*p+1], mps_large[2*p+1],
                                                                  left_[2*p+1], mpos_mix[L][2*p+1]);
            }
            
        }
        
#ifndef NDEBUG
        double energy = maquis::real(expval(mps_large, mpos_mix[L]));
        maquis::cout << "Energy " << "finegraining_final " << energy << std::endl;
#endif
        
    }
    
};

#endif
