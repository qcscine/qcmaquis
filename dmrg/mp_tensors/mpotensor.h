#ifndef MPOTENSOR_H
#define MPOTENSOR_H

#include "block_matrix/block_matrix.h"
#include "block_matrix/indexing.h"
#include "utils/function_objects.h"

#include <iostream>

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
//    , upper_i(ud), lower_i(ld)
    { }

    std::size_t aux_dim() const { return data_.size(); }
    
    scalar_type & operator()(std::size_t i, access_type j, access_type k) { return data_[i](j, k); }
    scalar_type const & operator()(std::size_t i, access_type j, access_type k) const { return data_[i](j, k); }
    
    friend struct contraction;
    
    std::vector<scalar_type> traces() const
    {
        std::vector<scalar_type> ret;
        std::transform(data_.begin(), data_.end(), back_inserter(ret),
                       utils::functor_trace());
        return ret;
    }
    
#ifdef HAVE_ALPS_HDF5
    void serialize(alps::hdf5::iarchive & ar)
    {
        ar >> alps::make_pvp("data", data_);
    }
    
    void serialize(alps::hdf5::oarchive & ar) const
    {
        ar << alps::make_pvp("data", data_);
    }
#endif
    
public:
    std::vector<block_matrix<Matrix, SymmGroup> > data_;
//    Index<SymmGroup> upper_i, lower_i;
};

template<class Matrix, class SymmGroup>
class MPOTensor
{
public:
    typedef typename Matrix::value_type scalar_type;
    typedef double real_type;
    typedef std::pair<typename SymmGroup::charge, std::size_t> access_type;
    
    // the constructor asumes that the upper and lower physical dimension is the same
    MPOTensor(std::size_t = 1,
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
    
    block_matrix<Matrix, SymmGroup> const & operator()(std::size_t left_index,
                                                       std::size_t right_index) const;
    block_matrix<Matrix, SymmGroup> & operator()(std::size_t left_index,
                                                 std::size_t right_index);
    
    void multiply_by_scalar(scalar_type);
    
private:
    std::vector<block_matrix<Matrix, SymmGroup> > data_;
    std::size_t left_i, right_i;
//    Index<SymmGroup> phys_i;
};

#include "mp_tensors/mpotensor.hpp"

#endif
