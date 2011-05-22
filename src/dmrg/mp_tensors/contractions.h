/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef CONTRACTIONS_H
#define CONTRACTIONS_H

#include "utils/zout.hpp"

#include "mp_tensors/mpstensor.h"
#include "mp_tensors/mpotensor.h"

#include "mp_tensors/reshapes.h"
#include "block_matrix/indexing.h"

#include "utils/iterator_blas1.h"
#include "utils/sizeof.h"

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
                if (!mpo.has(b1, b2))
                    continue;
                
                block_matrix<Matrix, SymmGroup> const & W = mpo(b1, b2);
                
                if (W.n_blocks() == 0)
                    continue;
                
                block_matrix<Matrix, SymmGroup> T;
                gemm(t[b1], right.data_[b2], T);
                
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
                                
                                ret.data_.match_and_add_block(oblock, out_l_charge, out_r_charge);
                            }
            }
        
#ifndef NDEBUG
        ket_tensor.make_left_paired();
        assert(ret.data_.left_basis() == ket_tensor.data_.left_basis());
        assert(ret.data_.right_basis() == ket_tensor.data_.right_basis());
#endif
        
        return ret;
    }
    
    template<class Matrix, class SymmGroup>
    static Boundary<Matrix, SymmGroup>
    left_boundary_tensor_mpo(MPSTensor<Matrix, SymmGroup> const & mps,
                             Boundary<Matrix, SymmGroup> const & left,
                             MPOTensor<Matrix, SymmGroup> const & mpo)
    {
        static Timer
        reshape_timer("Reshape in lbtm"),
        loop_timer("Prod with MPO in lbtm"),
        loop1_timer("Prod of left and MPS in lbtm");
        
        mps.make_right_paired();
        
        std::vector<block_matrix<Matrix, SymmGroup> > t(left.aux_dim());
        
        loop1_timer.begin();
        size_t loop_max = left.aux_dim();
#ifndef MPI_PARALLEL
#pragma omp parallel for schedule(guided)
#endif
        for (std::size_t b = 0; b < loop_max; ++b) {
            block_matrix<Matrix, SymmGroup> tmp;
            gemm(transpose(left.data_[b]), mps.data_, tmp);
            reshape_timer.begin();
            reshape_right_to_left<Matrix>(mps.site_dim(), left.data_[b].right_basis(), mps.col_dim(),
                                          tmp, t[b]);
            reshape_timer.end();
        }
        loop1_timer.end();
        
        Index<SymmGroup> physical_i = mps.site_dim(), left_i = mps.row_dim(), right_i = mps.col_dim();
        
        Boundary<Matrix, SymmGroup> ret;
        ret.data_.resize(mpo.col_dim());
        
        typedef typename SymmGroup::charge charge;
        typedef std::size_t size_t;
        
        loop_timer.begin();
        
        mps.make_left_paired();
        
        loop_max = mpo.col_dim();
#ifndef MPI_PARALLEL
#pragma omp parallel for schedule(guided)
#endif
        for (size_t b2 = 0; b2 < loop_max; ++b2) {
            for (int run = 0; run < 2; ++run) {
                if (run == 1)
                    ret.data_[b2].allocate_blocks();
                bool pretend = (run == 0);
                
                for (size_t b1 = 0; b1 < left.aux_dim(); ++b1) {
                    if (!mpo.has(b1, b2))
                        continue;
                    
                    block_matrix<Matrix, SymmGroup> const & W = mpo(b1, b2);
                    if (W.n_blocks() == 0)
                        continue;
                    
                    block_matrix<Matrix, SymmGroup> const & T = t[b1];
                    
                    // note to self: I think these are always the same
                    ProductBasis<SymmGroup> out_left_pb(physical_i, left_i);
                    ProductBasis<SymmGroup> in_left_pb(physical_i, left.data_[b1].right_basis());
                    
                    for (size_t w_block = 0; w_block < W.n_blocks(); ++w_block)
                    {
                        assert( physical_i.has(W.left_basis()[w_block].first) );
                        assert( physical_i.has(W.right_basis()[w_block].first) );
                        
                        size_t s1 = physical_i.position(W.left_basis()[w_block].first);
                        size_t s2 = physical_i.position(W.right_basis()[w_block].first);
                        
                        for (size_t t_block = 0; t_block < T.n_blocks(); ++t_block)
                        {
                            size_t r = right_i.position(T.right_basis()[t_block].first);
                            size_t l = left_i.position(SymmGroup::fuse(T.left_basis()[t_block].first,
                                                                       -physical_i[s1].first));
                            
                            if (l >= left_i.size())
                                continue;
                            if (r >= right_i.size())
                                continue;
                            
                            {
                                charge T_l_charge = SymmGroup::fuse(physical_i[s1].first, left_i[l].first);
                                charge T_r_charge = right_i[r].first;
                                
                                if (! T.has_block(T_l_charge, T_r_charge) )
                                    continue;
                                
                                charge out_l_charge = SymmGroup::fuse(physical_i[s2].first, left_i[l].first);
                                charge out_r_charge = right_i[r].first;
                                
                                if (! left.data_[b1].right_basis().has(left_i[l].first) )
                                    continue;
                                if (! mps.col_dim().has(right_i[r].first) )
                                    continue;
                                if (! mps.data().left_basis().has(out_l_charge) )
                                    continue;
                                
                                size_t in_left_offset = in_left_pb(physical_i[s1].first, left_i[l].first);
                                size_t out_left_offset = out_left_pb(physical_i[s2].first, left_i[l].first);
                                
                                if (!pretend) {
                                    Matrix const & wblock = W(physical_i[s1].first, physical_i[s2].first);
                                    Matrix const & iblock = T(T_l_charge, T_r_charge);
                                    Matrix & oblock = ret.data_[b2](out_l_charge, out_r_charge);
                                    
                                    #ifdef MPI_PARALLEL
                                    ambient::push(ambient::lb_tensor_mpo_l<Matrix::value_type>, ambient::lb_tensor_mpo_c<Matrix::value_type>,
                                                  oblock, iblock, wblock, out_left_offset, in_left_offset, 
                                                  physical_i[s1].second, physical_i[s2].second, left_i[l].second, right_i[r].second);
                                    #else
                                    for(size_t ss1 = 0; ss1 < physical_i[s1].second; ++ss1)
                                        for(size_t ss2 = 0; ss2 < physical_i[s2].second; ++ss2) {
                                            typename Matrix::value_type wblock_t = wblock(ss1, ss2);
                                            for(size_t rr = 0; rr < right_i[r].second; ++rr) {
                                                iterator_axpy(&iblock(in_left_offset + ss1*left_i[l].second, rr),
                                                              &iblock(in_left_offset + ss1*left_i[l].second, rr) + left_i[l].second, // bugbug
                                                              &oblock(out_left_offset + ss2*left_i[l].second, rr),
                                                              wblock_t);
                                            }
                                        }
                                    #endif
                                }
                                
                                if (pretend)
                                    ret.data_[b2].reserve(out_l_charge, out_r_charge,
                                                          mps.data().left_basis().size_of_block(out_l_charge),
                                                          right_i[r].second);
                            }
                        }
                    }
                }
            }
        }
        loop_timer.end();
        
        return ret;
    }
    
    template<class Matrix, class SymmGroup>
    static Boundary<Matrix, SymmGroup>
    right_boundary_tensor_mpo(MPSTensor<Matrix, SymmGroup> const & mps,
                              Boundary<Matrix, SymmGroup> const & right,
                              MPOTensor<Matrix, SymmGroup> const & mpo)
    {
        static Timer timer("right_boundary_tensor_mpo");
        timer.begin();
        
        mps.make_left_paired();
        
        std::vector<block_matrix<Matrix, SymmGroup> > t(right.aux_dim());
        
        size_t loop_max = right.aux_dim();
#ifndef MPI_PARALLEL
#pragma omp parallel for schedule(guided)
#endif
        for(std::size_t b = 0; b < loop_max; ++b){
            gemm(mps.data_, right.data_[b], t[b]);
            block_matrix<Matrix, SymmGroup> tmp;
            reshape_left_to_right<Matrix>(mps.site_dim(), mps.row_dim(), right.data_[b].right_basis(),
                                          t[b], tmp);
            swap(t[b], tmp);
        }
        
        Index<SymmGroup> physical_i = mps.site_dim(), left_i = mps.row_dim(), right_i = mps.col_dim();
        
        Boundary<Matrix, SymmGroup> ret;
        ret.data_.resize(mpo.row_dim());
        
        typedef typename SymmGroup::charge charge;
        typedef std::size_t size_t;
        
        mps.make_right_paired();
        
        loop_max = mpo.row_dim();
#ifndef MPI_PARALLEL
#pragma omp parallel for schedule(guided)
#endif
        for(size_t b1 = 0; b1 < loop_max; ++b1) {
            for(int run = 0; run < 2; ++run) {
                if(run == 1)
                    ret.data_[b1].allocate_blocks();
                bool pretend = (run == 0);
                for(size_t b2 = 0; b2 < mpo.col_dim(); ++b2)
                {
                    if(!mpo.has(b1, b2))
                        continue;
                    
                    block_matrix<Matrix, SymmGroup> const & W = mpo(b1, b2);
                    if (W.n_blocks() == 0)
                        continue;
                    
                    block_matrix<Matrix, SymmGroup> const & T = t[b2];
                    
                    ProductBasis<SymmGroup> out_right_pb(physical_i, right_i,
                                                         boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                                             -boost::lambda::_1, boost::lambda::_2));
                    ProductBasis<SymmGroup> in_right_pb(physical_i, right.data_[b2].right_basis(),
                                                        boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                                            -boost::lambda::_1, boost::lambda::_2));
                    
                    for (size_t w_block = 0; w_block < W.n_blocks(); ++w_block)
                    {
                        size_t s1 = physical_i.position(W.left_basis()[w_block].first);
                        size_t s2 = physical_i.position(W.right_basis()[w_block].first);
                        
                        for (size_t t_block = 0; t_block < T.n_blocks(); ++t_block)
                        {
                            size_t l = left_i.position(T.left_basis()[t_block].first);
                            size_t r = right_i.position(SymmGroup::fuse(physical_i[s1].first,
                                                                        T.right_basis()[t_block].first));
                            
                            if (l >= left_i.size())
                                continue;
                            if (r >= right_i.size())
                                continue;
                            
                            {   
                                charge T_l_charge = T.left_basis()[t_block].first;
                                charge T_r_charge = T.right_basis()[t_block].first;
                                
                                charge out_l_charge = left_i[l].first;
                                charge out_r_charge = SymmGroup::fuse(-physical_i[s2].first,
                                                                      right_i[r].first);
                                
                                assert( T_r_charge == SymmGroup::fuse(-physical_i[s1].first,
                                                                      right_i[r].first) );
                                
                                if (! right.data_[b2].right_basis().has(right_i[r].first) )
                                    continue;
                                if (! mps.row_dim().has(left_i[l].first) )
                                    continue;
                                if (! mps.data().right_basis().has(out_r_charge) )
                                    continue;
                                
                                size_t in_right_offset = in_right_pb(physical_i[s1].first, right_i[r].first);
                                size_t out_right_offset = out_right_pb(physical_i[s2].first, right_i[r].first);
                                
                                if (!pretend) {
                                    const Matrix & wblock = W(physical_i[s1].first, physical_i[s2].first);
                                    const Matrix & iblock = T(T_l_charge, T_r_charge);
                                    Matrix & oblock = ret.data_[b1](out_l_charge, out_r_charge);

                                    //printf("contraction: %d %d , %d %d\n", oblock.num_rows(), oblock.num_cols(), iblock.num_rows(), iblock.num_cols());
                                    #ifdef MPI_PARALLEL
                                    ambient::push(ambient::rb_tensor_mpo_l<Matrix::value_type>, ambient::rb_tensor_mpo_c<Matrix::value_type>,
                                                  oblock, iblock, wblock, out_right_offset, in_right_offset, 
                                                  physical_i[s1].second, physical_i[s2].second, left_i[l].second, right_i[r].second);
                                    #else
                                    for (size_t ss1 = 0; ss1 < physical_i[s1].second; ++ss1)
                                        for (size_t ss2 = 0; ss2 < physical_i[s2].second; ++ss2) {
                                            typename Matrix::value_type wblock_t = wblock(ss1, ss2);
                                            for (size_t rr = 0; rr < right_i[r].second; ++rr)
                                                for (size_t ll = 0; ll < left_i[l].second; ++ll) {
                                                    oblock(ll, out_right_offset + ss2*right_i[r].second+rr) +=
                                                    iblock(ll, in_right_offset + ss1*right_i[r].second+rr) * wblock_t;
                                                }
                                        }
                                    #endif
                                }
                                
                                if (pretend)
                                    ret.data_[b1].reserve(out_l_charge, out_r_charge,
                                                          left_i[l].second,
                                                          mps.data().right_basis().size_of_block(out_r_charge));
                            }
                        }
                    }
                }
            }
        }
        
        timer.end();
        return ret;
    }
    
    
    
    template<class Matrix, class SymmGroup>
    static Boundary<Matrix, SymmGroup>
    overlap_mpo_left_step(MPSTensor<Matrix, SymmGroup> const & bra_tensor,
                          MPSTensor<Matrix, SymmGroup> const & ket_tensor,
                          Boundary<Matrix, SymmGroup> const & left,
                          MPOTensor<Matrix, SymmGroup> const & mpo)
    {
        static Timer timer("overlap_mpo_left_step");
        timer.begin();
        
        Boundary<Matrix, SymmGroup> lbtm = left_boundary_tensor_mpo(ket_tensor, left, mpo);
        
        bra_tensor.make_left_paired();
        Boundary<Matrix, SymmGroup> ret;
        ret.data_.resize(mpo.col_dim());
        
        std::size_t loop_max = mpo.col_dim();
#ifndef MPI_PARALLEL
#pragma omp parallel for schedule(guided)
#endif
        for (std::size_t b = 0; b < loop_max; ++b)
            gemm(transpose(lbtm.data_[b]), bra_tensor.data(), ret.data_[b]);
        
        timer.end();
        
        return ret;
    }
    
    template<class Matrix, class SymmGroup>
    static Boundary<Matrix, SymmGroup>
    overlap_mpo_right_step(MPSTensor<Matrix, SymmGroup> const & bra_tensor,
                           MPSTensor<Matrix, SymmGroup> const & ket_tensor,
                           Boundary<Matrix, SymmGroup> const & right,
                           MPOTensor<Matrix, SymmGroup> const & mpo)
    {
        static Timer timer("overlap_mpo_right_step");
        timer.begin();
        
        Boundary<Matrix, SymmGroup> rbtm = right_boundary_tensor_mpo(ket_tensor, right, mpo);
        
        bra_tensor.make_right_paired();
        Boundary<Matrix, SymmGroup> ret;
        ret.data_.resize(mpo.row_dim());
        
        std::size_t loop_max = mpo.row_dim();
        block_matrix<Matrix, SymmGroup> tmp = transpose(bra_tensor.data());
#ifndef MPI_PARALLEL
#pragma omp parallel for schedule(guided)
#endif
        for (std::size_t b = 0; b < loop_max; ++b)
            gemm(rbtm.data_[b], tmp, ret.data_[b]);
        
        timer.end();
        
        return ret;
    }
    
    template<class Matrix, class SymmGroup>
    static MPSTensor<Matrix, SymmGroup>
    site_hamil2(MPSTensor<Matrix, SymmGroup> const & ket_tensor,
                Boundary<Matrix, SymmGroup> const & left,
                Boundary<Matrix, SymmGroup> const & right,
                MPOTensor<Matrix, SymmGroup> const & mpo)
    {
        static Timer lbtm("lbtm in site_hamil2"), loop("loop in site_hamil2"), all("site_hamil all");
        all.begin();
        
        lbtm.begin();
        Boundary<Matrix, SymmGroup> left_mpo_mps = left_boundary_tensor_mpo(ket_tensor, left, mpo);
        lbtm.end();
        MPSTensor<Matrix, SymmGroup> ret = ket_tensor;
        ret.multiply_by_scalar(0);
        ret.make_left_paired();
        
        typedef typename SymmGroup::charge charge;
        typedef std::size_t size_t;
        
        size_t loop_max = mpo.col_dim();
        
        loop.begin();
#ifndef MPI_PARALLEL
#pragma omp parallel for schedule(guided)
#endif
        for (size_t b = 0; b < loop_max; ++b)
        {
            block_matrix<Matrix, SymmGroup> oblock;
            gemm(left_mpo_mps.data_[b], right.data_[b], oblock);
            for (size_t k = 0; k < oblock.n_blocks(); ++k)
#ifndef MPI_PARALLEL
#pragma omp critical
#endif
                ret.data_.match_and_add_block(oblock[k],
                                              oblock.left_basis()[k].first,
                                              oblock.right_basis()[k].first);
            
        }
        loop.end();
        
        all.end();
        
        return ret;
    }
    
    template<class Matrix, class SymmGroup>
    static MPSTensor<Matrix, SymmGroup>
    predict_new_state_l2r_sweep(MPSTensor<Matrix, SymmGroup> const & mps,
                                MPOTensor<Matrix, SymmGroup> const & mpo,
                                Boundary<Matrix, SymmGroup> const & left,
                                Boundary<Matrix, SymmGroup> const & right,
                                double alpha, double cutoff, std::size_t Mmax,
                                Logger & logger)
    {
        static Timer timer("predict_new_state_l2r_sweep");
        timer.begin();
        
        mps.make_left_paired();
        block_matrix<Matrix, SymmGroup> dm;
        gemm(mps.data_, transpose(mps.data_), dm);
        
        Boundary<Matrix, SymmGroup> half_dm = left_boundary_tensor_mpo(mps, left, mpo);
        
        mps.make_left_paired();
        for (std::size_t b = 0; b < half_dm.aux_dim(); ++b)
        {
            block_matrix<Matrix, SymmGroup> tdm;
            gemm(half_dm.data_[b], transpose(half_dm.data_[b]), tdm);
            
            tdm *= alpha;
            for (std::size_t k = 0; k < tdm.n_blocks(); ++k) {
                if (mps.data_.left_basis().has(tdm.left_basis()[k].first))
                    dm.match_and_add_block(tdm[k],
                                           tdm.left_basis()[k].first,
                                           tdm.right_basis()[k].first);
            }
        }
        
        mps.make_left_paired();
        assert(dm.left_basis() == mps.data_.left_basis());
        
        block_matrix<Matrix, SymmGroup> U, V;
        block_matrix<typename blas::associated_diagonal_matrix<Matrix>::type, SymmGroup> S;
        
        syev_truncate(dm, U, S, cutoff, Mmax, logger);
        
        MPSTensor<Matrix, SymmGroup> ret = mps;
        ret.data_ = U;
        ret.right_i = U.right_basis();
        
        timer.end();
        return ret;
    }
    
    template<class Matrix, class SymmGroup>
    static MPSTensor<Matrix, SymmGroup>
    predict_lanczos_l2r_sweep(MPSTensor<Matrix, SymmGroup> B,
                              MPSTensor<Matrix, SymmGroup> const & psi,
                              MPSTensor<Matrix, SymmGroup> const & A)
    {
        static Timer timer("predict_lanczos_l2r_sweep");
        timer.begin();
        
        psi.make_left_paired();
        A.make_left_paired();
        
        block_matrix<Matrix, SymmGroup> tmp;
        gemm(transpose(A.data_), psi.data_, tmp);
        
        B.multiply_from_left(tmp);
        
        timer.end();
        return B;
    }
    
    template<class Matrix, class SymmGroup>
    static MPSTensor<Matrix, SymmGroup>
    predict_new_state_r2l_sweep(MPSTensor<Matrix, SymmGroup> const & mps,
                                MPOTensor<Matrix, SymmGroup> const & mpo,
                                Boundary<Matrix, SymmGroup> const & left,
                                Boundary<Matrix, SymmGroup> const & right,
                                double alpha, double cutoff, std::size_t Mmax,
                                Logger & logger)
    {
        static Timer timer("predict_new_state_r2l_sweep");
        timer.begin();
        
        mps.make_right_paired();
        block_matrix<Matrix, SymmGroup> dm;
        gemm(transpose(mps.data_), mps.data_, dm);
            
        Boundary<Matrix, SymmGroup> half_dm = right_boundary_tensor_mpo(mps, right, mpo);
        
        mps.make_right_paired();
        for (std::size_t b = 0; b < half_dm.aux_dim(); ++b)
        {
            block_matrix<Matrix, SymmGroup> tdm;
            gemm(transpose(half_dm.data_[b]), half_dm.data_[b], tdm);
            
            tdm *= alpha;
            for (std::size_t k = 0; k < tdm.n_blocks(); ++k) {
                if (mps.data_.right_basis().has(tdm.right_basis()[k].first))
                    dm.match_and_add_block(tdm[k],
                                           tdm.left_basis()[k].first,
                                           tdm.right_basis()[k].first);
            }
        }
        
        mps.make_right_paired();
        assert(dm.right_basis() == mps.data_.right_basis());
        
        block_matrix<Matrix, SymmGroup> U, V;
        block_matrix<typename blas::associated_diagonal_matrix<Matrix>::type, SymmGroup> S;
        
        syev_truncate(dm, U, S, cutoff, Mmax, logger);
        V = transpose(U);
        
        MPSTensor<Matrix, SymmGroup> ret = mps;
        ret.data_ = V;
        ret.left_i = V.left_basis();
        
        timer.end();
        return ret; 
    }
    
    template<class Matrix, class SymmGroup>
    static MPSTensor<Matrix, SymmGroup>
    predict_lanczos_r2l_sweep(MPSTensor<Matrix, SymmGroup> B,
                              MPSTensor<Matrix, SymmGroup> const & psi,
                              MPSTensor<Matrix, SymmGroup> const & A)
    {
        static Timer timer("predict_lanczos_r2l_sweep");
        timer.begin();
        
        psi.make_right_paired();
        A.make_right_paired();
        
        block_matrix<Matrix, SymmGroup> tmp;
        gemm(psi.data_, transpose(A.data_), tmp);
        
        B.multiply_from_right(tmp);
        
        timer.end();
        return B;
    }
    
    template<class Matrix, class SymmGroup>
    static block_matrix<Matrix, SymmGroup>
    multiply_with_twosite(block_matrix<Matrix, SymmGroup> const & vec,
                          block_matrix<Matrix, SymmGroup> const & op,
                          Index<SymmGroup> const & left_i,
                          Index<SymmGroup> const & right_i,
                          Index<SymmGroup> const & phys_i)
    {
        using std::size_t;
        typedef typename SymmGroup::charge charge;
        
        // they should have the same structure
        block_matrix<Matrix, SymmGroup> ret(vec.left_basis(), vec.right_basis());
        
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
                                charge lc = left_i[ls].first, rc = right_i[rs].first, lpc = phys_i[lps].first, rpc = phys_i[rps].first;
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
                                
        //                        cout << "Yes!" << endl;
                                
                                charge both_charge = SymmGroup::fuse(lpc, rpc);
                                
                                assert( op.has_block(both_charge, both_charge) );
                                
                                Matrix & oblock = ret(left_out_charge, right_out_charge);
                                Matrix const & iblock = vec(left_vec_charge, right_vec_charge);
                                Matrix const & op_block = op(both_charge, both_charge);
                                
                                size_t i_l_offset = left_pb(ilpc, lc);
                                size_t i_r_offset = right_pb(irpc, rc);
                                size_t l_offset = left_pb(lpc, lc);
                                size_t r_offset = right_pb(rpc, rc);
                                size_t i_op_offset = phys_pb(ilpc, irpc);
                                size_t op_offset = phys_pb(lpc, rpc);
                                
//                                cout << ilpc << " " << irpc << " " << lpc << " " << rpc << endl;
                                
                                for (size_t ll = 0; ll < left_i[ls].second; ++ll)
                                    for (size_t rr = 0; rr < right_i[rs].second; ++rr)
                                        for (size_t lp = 0; lp < phys_i[lps].second; ++lp)
                                            for (size_t rp = 0; rp < phys_i[rps].second; ++rp)
                                                for (size_t ilp = 0; ilp < phys_i[lps].second; ++ilp)
                                                    for (size_t irp = 0; irp < phys_i[rps].second; ++irp) {
        //                                                cout << "Yo" << endl;
//                                                        cout << i_op_offset + ilp*phys_i[irps].second + irp << " " << op_offset + lp*phys_i[rps].second + rp << endl;
//                                                        oblock(l_offset + lp*left_i[ls].second + ll,
//                                                               r_offset + rp*right_i[rs].second + rr);
//                                                        iblock(i_l_offset + ilp*left_i[ls].second + ll,
//                                                               i_r_offset + irp*right_i[rs].second + rr);
//                                                        op_block(i_op_offset + ilp*phys_i[irps].second + irp,
//                                                                 op_offset + lp*phys_i[rps].second + rp);
                                                        oblock(l_offset + lp*left_i[ls].second + ll,
                                                               r_offset + rp*right_i[rs].second + rr) += 
                                                        iblock(i_l_offset + ilp*left_i[ls].second + ll,
                                                               i_r_offset + irp*right_i[rs].second + rr) *
                                                        op_block(i_op_offset + ilp*phys_i[irps].second + irp,
                                                                 op_offset + lp*phys_i[rps].second + rp);
                                                    }
                            }
        
        return ret;
    }
};

#endif
