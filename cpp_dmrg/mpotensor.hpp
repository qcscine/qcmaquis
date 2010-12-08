#include "mpotensor.h"

#include "reshapes.h"

template<class Matrix, class SymmGroup>
MPOTensor<Matrix, SymmGroup>::MPOTensor(Index<SymmGroup> const & sd,
                                        Index<SymmGroup> const & ld,
                                        Index<SymmGroup> const & rd)
: phys_i(sd)
, left_i(ld)
, right_i(rd)
, data_(sd*ld, sd*rd)
, cur_storage(LeftUp)
, cur_normalization(U)
{
    data_.fill(drand48);
}

template<class Matrix, class SymmGroup>
typename MPOTensor<Matrix, SymmGroup>::scalar_type & 
MPOTensor<Matrix, SymmGroup>::operator()(MPOTensor<Matrix, SymmGroup>::access_type const & left_index,
                                         MPOTensor<Matrix, SymmGroup>::access_type const & right_index,
                                         MPOTensor<Matrix, SymmGroup>::access_type const & ket_index,
                                         MPOTensor<Matrix, SymmGroup>::access_type const & bra_index)
{
    return data_(calculate_index(phys_i ^ left_i,
                                 ket_index ^ left_index),
                 calculate_index(phys_i ^ right_i,
                                 bra_index ^ right_index));
}

template<class Matrix, class SymmGroup>
void MPOTensor<Matrix, SymmGroup>::reflect()
{
    typedef typename Index<SymmGroup>::basis_iterator bit;
    
    block_matrix<Matrix, SymmGroup> t(phys_i*right_i, phys_i*left_i);
    
    for (bit l1 = phys_i.basis_begin(); !l1.end(); ++l1)
        for (bit l2 = left_i.basis_begin(); !l2.end(); ++l2)
            for (bit r1 = phys_i.basis_begin(); !r1.end(); ++r1)
                for (bit r2 = right_i.basis_begin(); !r2.end(); ++r2)
                    t(calculate_index(phys_i ^ right_i,
                                      *l1 ^ *r2),
                      calculate_index(phys_i ^ left_i,
                                      *r1 ^ *l2))
                    =
                    data_(calculate_index(phys_i ^ left_i,
                                          *l1 ^ *l2),
                          calculate_index(phys_i ^ right_i,
                                          *r1 ^ *r2));
    
    swap(data_, t);
    swap(left_i, right_i);
}

template<class Matrix, class SymmGroup>
void MPOTensor<Matrix, SymmGroup>::multiply_by_scalar(scalar_type v)
{
    data_ *= v;
}

template<class Matrix, class SymmGroup>
Index<SymmGroup> MPOTensor<Matrix, SymmGroup>::row_dim() const
{
    return left_i;
}

template<class Matrix, class SymmGroup>
Index<SymmGroup> MPOTensor<Matrix, SymmGroup>::col_dim() const
{
    return right_i;
}
