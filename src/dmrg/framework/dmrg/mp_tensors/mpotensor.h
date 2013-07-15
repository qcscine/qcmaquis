/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MPOTENSOR_H
#define MPOTENSOR_H

#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/indexing.h"
#include "utils/function_objects.h"

#include <iostream>
#include <set>


namespace MPOTensor_detail
{
    struct pair_cmp
    {
        bool operator()(std::pair<std::size_t, std::size_t> const & i,
                        std::pair<std::size_t, std::size_t> const & j) const
        {
            if (i.first < j.first)
                return true;
            else if (i.first > j.first)
                return false;
            else
                return i.second < j.second;
        }
    };
}

template <class Matrix, class SymmGroup> class column_iterator;
template <class Matrix, class SymmGroup> class compressor;
template <class Matrix, class SymmGroup> class MPOIndexer;

template<class Matrix, class SymmGroup>
class MPOTensor
{
    typedef std::pair<std::size_t, std::size_t> key_t;
    typedef block_matrix<Matrix, SymmGroup> value_t;
    typedef std::map<key_t, value_t, MPOTensor_detail::pair_cmp> data_t;
    
    typedef std::set<std::size_t> used_set_t;
    
public:
    typedef typename Matrix::value_type value_type;
    typedef typename maquis::traits::scalar_type<Matrix>::type scalar_type;
    typedef std::pair<typename SymmGroup::charge, std::size_t> access_type;
    
    MPOTensor(std::size_t = 1, std::size_t = 1);
    
    std::size_t row_dim() const;
    std::size_t col_dim() const;
    
    value_type & operator()(std::size_t left_index,
                            std::size_t right_index,
                            access_type const & ket_index,
                            access_type const & bra_index);
    value_type const & operator()(std::size_t left_index,
                                  std::size_t right_index,
                                  access_type const & ket_index,
                                  access_type const & bra_index) const;
    
    block_matrix<Matrix, SymmGroup> const & operator()(std::size_t left_index,
                                                       std::size_t right_index) const;
    block_matrix<Matrix, SymmGroup> & operator()(std::size_t left_index,
                                                 std::size_t right_index);
    
    void multiply_by_scalar(const scalar_type&);
    void divide_by_scalar(const scalar_type&);
#ifdef AMBIENT
    void persist() const;
#endif
    
    bool has(std::size_t left_index, std::size_t right_index) const;

    friend class column_iterator<Matrix, SymmGroup>;
    friend class compressor<Matrix, SymmGroup>;
    friend class MPOIndexer<Matrix, SymmGroup>;
    
private:
    data_t data_;
    
    std::size_t left_i, right_i;
};

#include "dmrg/mp_tensors/mpotensor.hpp"


#endif
