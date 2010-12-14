#ifndef CONTRACTIONS_H
#define CONTRACTIONS_H

#include "mpstensor.h"
#include "mpotensor.h"

#include "reshapes.h"
#include "indexing.h"

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
        bra_tensor.make_right_paired();
        ket_tensor.make_right_paired();
        
        block_matrix<Matrix, SymmGroup> t1;
        std::vector<block_matrix<Matrix, SymmGroup> > t2(left.aux_dim());
        
        for (std::size_t b = 0; b < left.aux_dim(); ++b)
        {
            gemm(conjugate(transpose(bra_tensor.data_)), left.data_[b], t1);
            gemm(t1, ket_tensor.data_, t2[b]);
        }
        
        Boundary<Matrix, SymmGroup> ret(ket_tensor.col_dim(), bra_tensor.col_dim(), mpo.col_dim());
        
        typedef typename Index<SymmGroup>::basis_iterator bit;
        for (std::size_t b1 = 0; b1 < left.aux_dim(); ++b1)
            for (std::size_t b2 = 0; b2 < mpo.col_dim(); ++b2)
                for (bit u = ket_tensor.col_dim().basis_begin(); !u.end(); ++u)
                    for (bit d = bra_tensor.col_dim().basis_begin(); !d.end(); ++d)
                        for (bit s1 = ket_tensor.site_dim().basis_begin(); !s1.end(); ++s1)
                            for (bit s2 = bra_tensor.site_dim().basis_begin(); !s2.end(); ++s2) {
//                                cout << (*u).first << " " << (*d).first << endl;
//                                cout << calculate_index(ket_tensor.phys_i ^ ket_tensor.right_i,
//                                                        -(*s1) ^ *u).first << " ";
//                                cout << calculate_index(bra_tensor.phys_i ^ bra_tensor.right_i,
//                                                        -(*s2) ^ *d).first << endl;
//                                ret.data_[b2](*u, *d);
//                                t2[b1](calculate_index(ket_tensor.phys_i ^ ket_tensor.right_i,
//                                                       -(*s1) ^ *u),
//                                       calculate_index(bra_tensor.phys_i ^ bra_tensor.right_i,
//                                                       -(*s2) ^ *d));
                                ret.data_[b2](*u, *d) +=
                                t2[b1](calculate_index(ket_tensor.phys_i ^ ket_tensor.right_i,
                                                       -(*s1) ^ *u),
                                       calculate_index(bra_tensor.phys_i ^ bra_tensor.right_i,
                                                       -(*s2) ^ *d))
                                * mpo(b1, b2, *s1, *s2);
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
        MPSTensor<Matrix, SymmGroup>
        r_bra_tensor = bra_tensor.get_reflected(),
        r_ket_tensor = ket_tensor.get_reflected();
        MPOTensor<Matrix, SymmGroup>
        r_mpo = mpo.get_reflected();
        
        return overlap_mpo_left_step(r_bra_tensor, r_ket_tensor, right, r_mpo);
    }
    
    template<class Matrix, class SymmGroup>
    static MPSTensor<Matrix, SymmGroup>
    site_hamil(MPSTensor<Matrix, SymmGroup> const & ket_tensor,
               Boundary<Matrix, SymmGroup> const & left,
               Boundary<Matrix, SymmGroup> const & right,
               MPOTensor<Matrix, SymmGroup> const & mpo)
    {
        ket_tensor.make_right_paired();
        
        std::vector<MPSTensor<Matrix, SymmGroup> >
        t(left.aux_dim(),
          MPSTensor<Matrix, SymmGroup>(ket_tensor.site_dim(),
                                       left.lower_dim(),
                                       ket_tensor.col_dim(), false));
        
        for (std::size_t b = 0; b < left.aux_dim(); ++b) {
            gemm(transpose(left.data_[b]), ket_tensor.data_, t[b].data_);
            t[b].cur_storage = RightPaired;
            t[b].make_left_paired();
        }
        
        MPSTensor<Matrix, SymmGroup> ret(ket_tensor.site_dim(),
                                         left.lower_dim(), right.lower_dim(), false);
        ret.make_left_paired();
        ret.data_ = transpose(ret.data_);
        
        typedef typename Index<SymmGroup>::basis_iterator bit;
        
        for (std::size_t b1 = 0; b1 < left.aux_dim(); ++b1)
            for (std::size_t b2 = 0; b2 < right.aux_dim(); ++b2) {
                block_matrix<Matrix, SymmGroup> tmp;
                gemm(t[b1].data_, right.data_[b2], tmp);
                tmp = transpose(tmp);
                
                for (bit l = left.lower_dim().basis_begin(); !l.end(); ++l)
                    for (bit s1 = ket_tensor.site_dim().basis_begin(); !s1.end(); ++s1)
                        for (bit s2 = ket_tensor.site_dim().basis_begin(); !s2.end(); ++s2) {
                            std::pair<typename SymmGroup::charge, std::size_t>
                            i1 = calculate_index(ket_tensor.site_dim() ^ ket_tensor.left_i, *s2 ^ *l),
                            i2 = calculate_index(ket_tensor.site_dim() ^ ket_tensor.left_i, *s1 ^ *l);
                            
                            typename Matrix::value_type mpo_v = mpo(b1, b2, *s1, *s2);

                            for (typename Index<SymmGroup>::const_iterator sectors = right.lower_dim().begin();
                                 sectors != right.lower_dim().end(); ++sectors)
                            {
                                std::size_t ret_block = ret.data_.left_basis().position(sectors->first);
                                std::size_t tmp_block = tmp.left_basis().position(sectors->first);
                                
                                Matrix &m1 = ret.data_[ret_block], m2 = tmp[tmp_block];
                                
                                for (std::size_t elem = 0; elem < sectors->second; ++elem)
                                    m1(elem, i1.second) += m2(elem, i2.second) * mpo_v;
                            }
                        }
            }
        
        ret.data_ = conjugate(transpose(ret.data_));
        
        return ret;
    }
};

#endif
