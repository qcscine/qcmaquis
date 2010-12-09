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
        
        block_matrix<Matrix, SymmGroup> t1, t2 = conjugate_transpose(bra_tensor.data_), t3 = transpose(right);
        gemm(ket_tensor.data_, t3, t1);
        reshape_left_to_right(ket_tensor.phys_i, ket_tensor.left_i, right.left_basis(),
                              t1, t3);
        
        gemm(t3, t2, t1);
        return transpose(t1);
    }
    
    template<class Matrix, class SymmGroup>
    static block_matrix<Matrix, SymmGroup>
    mpo_left(MPOTensor<Matrix, SymmGroup> const & mpo,
             block_matrix<Matrix, SymmGroup> const & left,
             Index<SymmGroup> const & i2,
             Index<SymmGroup> const & o3)
    {   
        /*
         i2  i1
         |   |
         O-s-O - o1
         |   |
         o3   o2
         */
        
        Index<SymmGroup>
        i1 = mpo.phys_i, o1 = mpo.right_i, o2 = mpo.phys_i, s = mpo.left_i;
        
        block_matrix<Matrix, SymmGroup> ret(i1*i2, o1*o2*o3);
        
        typedef typename Index<SymmGroup>::basis_iterator bit;
        
        for (bit ii1 = i1.basis_begin(); !ii1.end(); ++ii1)
            for (bit ii2 = i2.basis_begin(); !ii2.end(); ++ii2)
                for (bit io1 = o1.basis_begin(); !io1.end(); ++io1)
                    for (bit io2 = o2.basis_begin(); !io2.end(); ++io2)
                        for (bit io3 = o3.basis_begin(); !io3.end(); ++io3)
                            for (bit is = s.basis_begin(); !is.end(); ++is) {
//                                ret(calculate_index(i1 ^ i2, *ii1 ^ *ii2),
//                                    calculate_index(o1 ^ o2 ^ o3, *io1 ^ *io2 ^ *io3));
//                                left(calculate_index(i2 ^ s, *ii2 ^ *is),
//                                     calculate_index(_(o3), _(*io3)));
                                ret(calculate_index(i1 ^ i2, *ii1 ^ *ii2),
                                    calculate_index(o1 ^ o2 ^ o3, *io1 ^ *io2 ^ *io3))
                                +=
                                left(calculate_index(i2 ^ s, *ii2 ^ *is),
                                     calculate_index(_(o3), _(*io3)))
                                * mpo(*is, *io1, *ii1, *io2);
                            }
        
        return ret;
    }
    
    template<class Matrix, class SymmGroup>
    static block_matrix<Matrix, SymmGroup>
    overlap_mpo_left_step(MPSTensor<Matrix, SymmGroup> const & bra_tensor,
                          MPSTensor<Matrix, SymmGroup> const & ket_tensor,
                          block_matrix<Matrix, SymmGroup> const & left,
                          MPOTensor<Matrix, SymmGroup> const & mpo)
    {
        bra_tensor.make_left_paired();
        ket_tensor.make_left_paired();
        
        block_matrix<Matrix, SymmGroup> left_mpo = mpo_left(mpo, left,
                                                            ket_tensor.left_i,
                                                            bra_tensor.left_i),
        braconj = conjugate(bra_tensor.data_);
        
        left_mpo = transpose(left_mpo);
        
        block_matrix<Matrix, SymmGroup> t, t2(ket_tensor.right_i * mpo.right_i,
                                              mpo.phys_i * bra_tensor.left_i);
        gemm(left_mpo, ket_tensor.data_, t);
        
        typedef typename Index<SymmGroup>::basis_iterator bit;
        
        /*
         / - O - o1         / - O - i1
         |   |              |   |
         O - O - i1 ===>>>  O - O - i2
         |   |              |   |
         i3  i2             o2  o1
         */
        
        for (bit o1_i1 = ket_tensor.right_i.basis_begin(); !o1_i1.end(); ++o1_i1)
            for (bit i1_i2 = mpo.right_i.basis_begin(); !i1_i2.end(); ++i1_i2)
                for (bit i2_o1 = mpo.phys_i.basis_begin(); !i2_o1.end(); ++i2_o1)
                    for (bit i3_o2 = bra_tensor.left_i.basis_begin(); !i3_o2.end(); ++i3_o2)
                        t2(calculate_index(ket_tensor.right_i ^ mpo.right_i, *o1_i1 ^ *i1_i2),
                           calculate_index(mpo.phys_i ^ bra_tensor.left_i, *i2_o1 ^ *i3_o2))
                        =
                        t(calculate_index(mpo.right_i ^ mpo.phys_i ^ bra_tensor.left_i,
                                          *i1_i2 ^ *i2_o1 ^ *i3_o2),
                          calculate_index(_(ket_tensor.right_i), _(*o1_i1)));
        
        gemm(t2, bra_tensor.data_, t);
        return t;
    }
    
    template<class Matrix, class SymmGroup>
    static block_matrix<Matrix, SymmGroup>
    overlap_mpo_right_step(MPSTensor<Matrix, SymmGroup> const & bra_tensor,
                           MPSTensor<Matrix, SymmGroup> const & ket_tensor,
                           block_matrix<Matrix, SymmGroup> right,
                           MPOTensor<Matrix, SymmGroup> const & mpo)
    {
        MPSTensor<Matrix, SymmGroup>
        r_bra_tensor = bra_tensor.get_reflected(),
        r_ket_tensor = ket_tensor.get_reflected();
        MPOTensor<Matrix, SymmGroup>
        r_mpo = mpo.get_reflected();
        
        right = overlap_mpo_left_step(r_bra_tensor, r_ket_tensor, right, r_mpo);
        
        return right;
    }
    
    template<class Matrix, class SymmGroup>
    static MPSTensor<Matrix, SymmGroup>
    site_hamil(MPSTensor<Matrix, SymmGroup> const & ket_tensor,
               block_matrix<Matrix, SymmGroup> left_mpo,
               block_matrix<Matrix, SymmGroup> const & right,
               MPOTensor<Matrix, SymmGroup> const & mpo)
    {
        ket_tensor.make_left_paired();
        
        left_mpo = transpose(left_mpo);
        
        block_matrix<Matrix, SymmGroup> t, t2(ket_tensor.right_i * mpo.right_i,
                                              mpo.phys_i * ket_tensor.left_i);
        gemm(left_mpo, ket_tensor.data_, t);
        
        typedef typename Index<SymmGroup>::basis_iterator bit;
        
        for (bit o1_i1 = ket_tensor.right_i.basis_begin(); !o1_i1.end(); ++o1_i1)
            for (bit i1_i2 = mpo.right_i.basis_begin(); !i1_i2.end(); ++i1_i2)
                for (bit i2_o1 = mpo.phys_i.basis_begin(); !i2_o1.end(); ++i2_o1)
                    for (bit i3_o2 = ket_tensor.left_i.basis_begin(); !i3_o2.end(); ++i3_o2)
                        t2(calculate_index(ket_tensor.right_i ^ mpo.right_i, *o1_i1 ^ *i1_i2),
                           calculate_index(mpo.phys_i ^ ket_tensor.left_i, *i2_o1 ^ *i3_o2))
                        =
                        t(calculate_index(mpo.right_i ^ mpo.phys_i ^ ket_tensor.left_i,
                                          *i1_i2 ^ *i2_o1 ^ *i3_o2),
                          calculate_index(_(ket_tensor.right_i), _(*o1_i1)));
        t2 = transpose(t2);
        gemm(t2, right, t);
        
        MPSTensor<Matrix, SymmGroup> ret(ket_tensor.site_dim(),
                                         ket_tensor.row_dim(), ket_tensor.col_dim());
        ret.data_ = conjugate(t);
        ret.cur_storage = LeftPaired;
        
//        cout << "Left dimensions: " << endl;
//        cout << "Ket: " << ket_tensor.site_dim()*ket_tensor.row_dim() << endl;
//        cout << "t: " << t.left_basis() << endl;
//        
//        cout << "right dimensions: " << endl;
//        cout << "Ket: " << ket_tensor.col_dim() << endl;
//        cout << "t: " << t.right_basis() << endl;
        
        return ret;
    }
};

#endif
