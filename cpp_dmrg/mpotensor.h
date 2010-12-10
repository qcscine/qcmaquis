#ifndef MPOTENSOR_H
#define MPOTENSOR_H

#include "block_matrix.h"
#include "indexing.h"

#include <iostream>

#include "general_matrix.hpp"

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
class Boundary
{
public:
    typedef typename Matrix::value_type scalar_type;
    typedef std::pair<typename SymmGroup::charge, std::size_t> access_type;
    
    Boundary(Index<SymmGroup> const & ud = Index<SymmGroup>(),
             Index<SymmGroup> const & ld = Index<SymmGroup>(),
             std::size_t ad = 1)
    : data_(ad, block_matrix<Matrix, SymmGroup>(ud, ld))
    , upper_i(ud), lower_i(ld)
    { }
    
    Index<SymmGroup> const & upper_dim() const { return upper_i; }
    Index<SymmGroup> const & lower_dim() const { return lower_i; }
    std::size_t aux_dim() const { return data_.size(); }
    
    scalar_type & operator()(std::size_t i, access_type j, access_type k) { return data_[i](j, k); }
    scalar_type const & operator()(std::size_t i, access_type j, access_type k) const { return data_[i](j, k); }
    
    friend struct contraction;
    
private:
    std::vector<block_matrix<Matrix, SymmGroup> > data_;
    Index<SymmGroup> upper_i, lower_i;
};

template<class Matrix, class SymmGroup>
class MPOTensor
{
public:
    typedef typename Matrix::value_type scalar_type;
    typedef double real_type;
    typedef std::pair<typename SymmGroup::charge, std::size_t> access_type;
    
    // the constructor asumes that the upper and lower physical dimension is the same
    MPOTensor(Index<SymmGroup> const & = Index<SymmGroup>(),
              std::size_t = 1,
              std::size_t = 1);
    
    Index<SymmGroup> const & site_bra_dim() const;
    Index<SymmGroup> const & site_ket_dim() const;
    std::size_t row_dim() const;
    std::size_t col_dim() const;
    
    scalar_type & operator()(std::size_t left_index,
                             std::size_t right_index,
                             access_type const & ket_index,
                             access_type const & bra_index);
    scalar_type const & operator()(std::size_t left_index,
                                   std::size_t right_index,
                                   access_type const & ket_index,
                                   access_type const & bra_index) const;
    
    void multiply_by_scalar(scalar_type);
    
    MPOTensor get_reflected() const;
    
private:
    blas::general_matrix<block_matrix<Matrix, SymmGroup> > data_;
    Index<SymmGroup> phys_i;
};
  
#include "mpotensor.hpp"

#endif
