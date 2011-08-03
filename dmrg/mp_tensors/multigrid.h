/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MULTIGRID_H
#define MULTIGRID_H

#include "utils/zout.hpp"

#include "mp_tensors/mpstensor.h"
#include "block_matrix/indexing.h"

struct multigrid {
    
    template<class Matrix, class SymmGroup>
    MPS<Matrix, SymmGroup>
    static restriction (MPS<Matrix, SymmGroup> const & mps_large,
                        MPS<Matrix, SymmGroup> & mps_small)
    {
        std::size_t L = mps_small.length();
        std::size_t LL = mps_large.length();
        assert(LL == 2*L);
        
        typedef typename SymmGroup::charge charge;
        typedef typename Index<SymmGroup>::basis_iterator bi_t;
        
        
        for (std::size_t p = 0; p < L; ++p)
        {
            //            std::cout << std::endl << std::endl << "********************" << std::endl << "starting loop for p = " << p << std::endl;
            
            block_matrix<Matrix, SymmGroup> M;
            
            Index<SymmGroup> alpha_basis = mps_large.row_dim(2*p);
            Index<SymmGroup> b_basis = mps_large.col_dim(2*p);
            Index<SymmGroup> beta_basis = mps_large.col_dim(2*p+1);
            Index<SymmGroup> s_basis = mps_small.site_dim(p);
            Index<SymmGroup> s1_basis = mps_large.site_dim(2*p);
            Index<SymmGroup> s2_basis = mps_large.site_dim(2*p+1);
            
            //            std::cout << "alpha_basis:" <<std::endl << alpha_basis << std::endl;
            //            std::cout << "beta_basis:" <<std::endl << beta_basis << std::endl;
            //            std::cout << "s_basis:" <<std::endl << s_basis << std::endl;
            //            std::cout << "s1_basis:" <<std::endl << s1_basis << std::endl;
            //            std::cout << "s2_basis:" <<std::endl << s2_basis << std::endl;
            
            
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
                                
                                //                            std::cout << "--" << std::endl;
                                //                            std::cout << "alpha: " << alpha->first << ", " << alpha->second << std::endl;
                                //                            std::cout << "beta: " << beta->first << ", " << beta->second << std::endl;
                                //                            std::cout << "s1: " << s1->first << ", " << s1->second << std::endl;
                                //                            std::cout << "s2: " << s2->first << ", " << s2->second << std::endl;
                                //                            std::cout << "s: " << s.first << ", " << s.second << std::endl;
                                
                                if (!M.has_block(out_left_c, out_right_c))
                                    M.insert_block(Matrix(out_left.size(s.first, alpha->first),
                                                          beta_basis.size_of_block(beta->first),
                                                          0),
                                                   out_left_c, out_right_c);
                                
                                
                                //                            std::cout << "block has size " << out_left.size(s1->first, alpha->first)
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
                                
                                //                            std::cout << "M[" << out_left_c << ", " << out_right_c << "]"
                                //                            << "(" << out_left_offset + s1->second*alpha_basis.size_of_block(alpha->first)+alpha->second << ", " << out_right_offset + s2->second*beta_basis.size_of_block(beta->first)+beta->second << ")" << std::endl;
                                //                            std::cout << " = " << M(std::make_pair(out_left_c, out_left_offset + s1->second*alpha_basis.size_of_block(alpha->first)+alpha->second),
                                //                                                    std::make_pair(out_right_c, out_right_offset + s2->second*beta_basis.size_of_block(beta->first)+beta->second)) << std::endl;
                                //                            
                                //                            std::cout << "mps[" << in_left_c << ", " << in_right_c << "]"
                                //                            << "(" << in_left_offset + s.second*alpha_basis.size_of_block(alpha->first) + alpha->second << ", " << in_right_offset + beta->second << ")" << std::endl;
                                //                            std::cout << " = " << mps_small[p].data()(std::make_pair(in_left_c, in_left_offset + s.second*alpha_basis.size_of_block(alpha->first) + alpha->second),
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
    static extension (MPS<Matrix, SymmGroup> const & mps_small,
                      MPS<Matrix, SymmGroup> & mps_large)
    {
        std::size_t L1 = mps_small.length();
        std::size_t L2 = mps_large.length();
        assert(L2 == 2*L1);
                        
        typedef typename SymmGroup::charge charge;
        typedef typename Index<SymmGroup>::basis_iterator bi_t;

        for (std::size_t p = 0; p < L1; ++p)
        {
//            std::cout << std::endl << std::endl << "********************" << std::endl << "starting loop for p = " << p << std::endl;
            
            block_matrix<Matrix, SymmGroup> M;
            
            Index<SymmGroup> alpha_basis = mps_small.row_dim(p);
            Index<SymmGroup> beta_basis = mps_small.col_dim(p);
            Index<SymmGroup> s_basis = mps_small.site_dim(p);
            Index<SymmGroup> s1_basis = mps_large.site_dim(2*p);
            Index<SymmGroup> s2_basis = mps_large.site_dim(2*p+1);
            
//            std::cout << "alpha_basis:" <<std::endl << alpha_basis << std::endl;
//            std::cout << "beta_basis:" <<std::endl << beta_basis << std::endl;
//            std::cout << "s_basis:" <<std::endl << s_basis << std::endl;
//            std::cout << "s1_basis:" <<std::endl << s1_basis << std::endl;
//            std::cout << "s2_basis:" <<std::endl << s2_basis << std::endl;
            
            
            mps_small[p].make_left_paired();

            ProductBasis<SymmGroup> in_left(s_basis, alpha_basis);
            
            ProductBasis<SymmGroup> out_left(s1_basis, alpha_basis);
            ProductBasis<SymmGroup> out_right(s2_basis, beta_basis,
                                              boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                                  -boost::lambda::_1, boost::lambda::_2));
            
//            std::cout << mps_small[p].data() << std::endl;
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
                                                        
                            if (!mps_small[p].data().has_block(in_left_c, in_right_c))
                                continue;

                            charge out_left_c = SymmGroup::fuse(s1->first, alpha->first);
                            charge out_right_c = SymmGroup::fuse(-s2->first, beta->first);
                            
//                            std::cout << "--" << std::endl;
//                            std::cout << "alpha: " << alpha->first << ", " << alpha->second << std::endl;
//                            std::cout << "beta: " << beta->first << ", " << beta->second << std::endl;
//                            std::cout << "s1: " << s1->first << ", " << s1->second << std::endl;
//                            std::cout << "s2: " << s2->first << ", " << s2->second << std::endl;
//                            std::cout << "s: " << s.first << ", " << s.second << std::endl;

                            if (!M.has_block(out_left_c, out_right_c))
                                M.insert_block(Matrix(out_left.size(s1->first, alpha->first),
                                                      out_right.size(s2->first, beta->first,
                                                                     boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                                                         -boost::lambda::_1, boost::lambda::_2)),
                                                      0),
                                               out_left_c, out_right_c);

                            
//                            std::cout << "block has size " << out_left.size(s1->first, alpha->first)
//                            << "x"
//                            << out_right.size(s2->first, beta->first,
//                                              boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
//                                                                  -boost::lambda::_1, boost::lambda::_2)) << std::endl;
                            
                            std::size_t in_left_offset = in_left(s.first, alpha->first);
                            std::size_t in_right_offset = 0;
                            
                            std::size_t out_left_offset = out_left(s1->first, alpha->first);
                            std::size_t out_right_offset = out_right(s2->first, beta->first);
                            
//                            std::cout << "M[" << out_left_c << ", " << out_right_c << "]"
//                            << "(" << out_left_offset + s1->second*alpha_basis.size_of_block(alpha->first)+alpha->second << ", " << out_right_offset + s2->second*beta_basis.size_of_block(beta->first)+beta->second << ")" << std::endl;
//                            std::cout << " = " << M(std::make_pair(out_left_c, out_left_offset + s1->second*alpha_basis.size_of_block(alpha->first)+alpha->second),
//                                                    std::make_pair(out_right_c, out_right_offset + s2->second*beta_basis.size_of_block(beta->first)+beta->second)) << std::endl;
//                            
//                            std::cout << "mps[" << in_left_c << ", " << in_right_c << "]"
//                            << "(" << in_left_offset + s.second*alpha_basis.size_of_block(alpha->first) + alpha->second << ", " << in_right_offset + beta->second << ")" << std::endl;
//                            std::cout << " = " << mps_small[p].data()(std::make_pair(in_left_c, in_left_offset + s.second*alpha_basis.size_of_block(alpha->first) + alpha->second),
//                                                                      std::make_pair(in_right_c, in_right_offset + beta->second)) << std::endl;
                            
                            M(std::make_pair(out_left_c, out_left_offset + s1->second*alpha_basis.size_of_block(alpha->first)+alpha->second),
                              std::make_pair(out_right_c, out_right_offset + s2->second*beta_basis.size_of_block(beta->first)+beta->second))
                            = mps_small[p].data()(std::make_pair(in_left_c, in_left_offset + s.second*alpha_basis.size_of_block(alpha->first) + alpha->second),
                                                  std::make_pair(in_right_c, in_right_offset + beta->second));
                            
                        }
            
            
            block_matrix<Matrix, SymmGroup> U, V, left, right;
            block_matrix<typename blas::associated_diagonal_matrix<Matrix>::type, SymmGroup> S, Ssqrt;
            
            svd(M, U, V, S);
            Ssqrt = sqrt(S);
            gemm(U, Ssqrt, left);
            gemm(Ssqrt, V, right);

            mps_large[2*p].data() = left;
            mps_large[2*p].left_i = alpha_basis;
            mps_large[2*p].right_i = left.right_basis();
            mps_large[2*p].cur_storage = LeftPaired;
            
            mps_large[2*p+1].data() = right;
            mps_large[2*p+1].right_i = beta_basis;
            mps_large[2*p+1].left_i = right.left_basis();
            mps_large[2*p+1].cur_storage = RightPaired;
            
        }
                
        mps_large.normalize_left();
    }
    
};

#endif
