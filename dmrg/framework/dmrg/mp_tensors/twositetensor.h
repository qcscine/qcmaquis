/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef TWOSITETENSOR_H
#define TWOSITETENSOR_H

#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/mpotensor.h"
#include "dmrg/block_matrix/indexing.h"
#include "dmrg/block_matrix/multi_index.h"
#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"

#include <iostream>
#include <algorithm>

enum TwoSiteStorageLayout {TSRightPaired, TSLeftPaired, TSBothPaired};

template<class Matrix, class SymmGroup>
class TwoSiteTensor
{
public:
    typedef std::size_t size_type;
    typedef typename MultiIndex<SymmGroup>::index_id index_id;
    typedef typename MultiIndex<SymmGroup>::set_id set_id;
    
    TwoSiteTensor(MPSTensor<Matrix, SymmGroup> const & mps1,
                  MPSTensor<Matrix, SymmGroup> const & mps2);

    TwoSiteTensor(MPSTensor<Matrix, SymmGroup> const & twin_mps);

    Index<SymmGroup> const & site_dim() const;
    Index<SymmGroup> const & row_dim() const;
    Index<SymmGroup> const & col_dim() const;
    
    block_matrix<Matrix, SymmGroup> & data();
    block_matrix<Matrix, SymmGroup> const & data() const;
    
    template<class Matrix_, class SymmGroup_>
    friend std::ostream& operator<<(std::ostream&, TwoSiteTensor<Matrix_, SymmGroup_> const &);

    TwoSiteTensor<Matrix, SymmGroup> & operator << ( MPSTensor<Matrix, SymmGroup> const & rhs);
    
    friend struct contraction;
    friend struct compression;
    friend struct multigrid;

    void make_left_paired() const;
    void make_both_paired() const;
    void make_right_paired() const;
    
    MPSTensor<Matrix, SymmGroup> make_mps() const;
    std::pair<MPSTensor<Matrix, SymmGroup>,
              MPSTensor<Matrix, SymmGroup> > split_mps_l2r(std::size_t Mmax, double cutoff) const;
    std::pair<MPSTensor<Matrix, SymmGroup>,
              MPSTensor<Matrix, SymmGroup> > split_mps_r2l(std::size_t Mmax, double cutoff) const;
    
    void swap_with(TwoSiteTensor & b);

    friend void swap(TwoSiteTensor & a, TwoSiteTensor & b)
    {
        a.swap_with(b);
    }
   
    template<class Archive> void load(Archive & ar);
    template<class Archive> void save(Archive & ar) const;
    
private:
    MultiIndex<SymmGroup> midx;
    set_id left_paired;
    set_id right_paired;
    set_id both_paired;

    Index<SymmGroup> phys_i, phys_i_orig, left_i, right_i;
    mutable block_matrix<Matrix, SymmGroup> data_;
    mutable TwoSiteStorageLayout cur_storage;
    Indicator cur_normalization;
};

#include "twositetensor.hpp"

#include "ts_ops.h"

#endif
