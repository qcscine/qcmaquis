#ifndef MPOTENSOR_H
#define MPOTENSOR_H

#include "block_matrix.h"
#include "indexing.h"

#include <iostream>

enum MPOStorageLayout { LeftUp, LeftDown };

/*
 LeftUp:
 
      i1
      |
 i2 - O - o2
      |
      o1
 
 */

template<class Matrix, class SymmGroup>
class MPOTensor
{
public:
    typedef typename Matrix::value_type scalar_type;
    typedef double real_type;
    typedef std::pair<typename SymmGroup::charge, std::size_t> access_type;
    
    // the constructor asumes that the upper and lower physical dimension is the same
    MPOTensor(Index<SymmGroup> const & sd = Index<SymmGroup>(),
              Index<SymmGroup> const & ld = Index<SymmGroup>(),
              Index<SymmGroup> const & rd = Index<SymmGroup>());
    
    Index<SymmGroup> site_bra_dim() const;
    Index<SymmGroup> site_ket_dim() const;
    Index<SymmGroup> row_dim() const;
    Index<SymmGroup> col_dim() const;
    
    scalar_type & operator()(access_type const & left_index,
                             access_type const & right_index,
                             access_type const & ket_index,
                             access_type const & bra_index);
    scalar_type const & operator()(access_type const & left_index,
                                   access_type const & right_index,
                                   access_type const & ket_index,
                                   access_type const & bra_index) const;
    
    friend struct contraction;
    
private:
    void make_leftup_paired();
    
    block_matrix<Matrix, SymmGroup> data_;
    Index<SymmGroup> phys_i, left_i, right_i;
    MPOStorageLayout cur_storage;
    Indicator cur_normalization;
};
  
#include "mpotensor.hpp"

#endif
