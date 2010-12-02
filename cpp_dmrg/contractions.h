#ifndef CONTRACTIONS_H
#define CONTRACTIONS_H

#include "mpstensor.h"
#include "mpotensor.h"

struct contraction {
    template<class Matrix, class SymmGroup>
    static block_matrix<Matrix, SymmGroup>
    overlap_left_step(MPSTensor<Matrix, SymmGroup> & bra_tensor,
                      MPSTensor<Matrix, SymmGroup> & ket_tensor,
                      block_matrix<Matrix, SymmGroup> & left,
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
    overlap_right_step(MPSTensor<Matrix, SymmGroup> & bra_tensor,
                       MPSTensor<Matrix, SymmGroup> & ket_tensor,
                       block_matrix<Matrix, SymmGroup> & right,
                       block_matrix<Matrix, SymmGroup> * localop = NULL)
    {
        if (localop != NULL)
            throw std::runtime_error("Not implemented!");
        
        assert(right.left_basis() == bra_tensor.right_i);
        assert(right.right_basis() == ket_tensor.right_i);
        
        bra_tensor.make_right_paired();
        ket_tensor.make_left_paired();
        
        block_matrix<Matrix, SymmGroup> t1, t2 = conjugate_transpose(bra_tensor.data_), t3 = transpose(right);
        gemm(ket_tensor.data_, t3, t1);
        reshape_left_to_right(ket_tensor.phys_i, ket_tensor.left_i, right.left_basis(),
                              t1, t3);
        
        gemm(t3, t2, t1);
        return transpose(t1);
    }
    
//    template<class Matrix, class SymmGroup>
//    static block_matrix<Matrix, SymmGroup>
//    overlap_mpo_left_step(MPSTensor<Matrix, SymmGroup> & bra_tensor,
//                          MPSTensor<Matrix, SymmGroup> & ket_tensor,
//                          block_matrix<Matrix, SymmGroup> & left,
//                          MPOTensor<Matrix, SymmGroup> & mpo)
//    {
//        block_matrix<Matrix, SymmGroup> t1(bra_tensor.left_i * mpo.phys_i * mpo.right_i, mpo.phys_i * ket_tensor.left_i);
//        
//    }
};

#endif
