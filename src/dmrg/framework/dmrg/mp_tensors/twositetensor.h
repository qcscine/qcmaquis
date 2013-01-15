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
              MPSTensor<Matrix, SymmGroup> > split_mps_l2r(std::size_t Mmax, double cutoff, Logger * iter_log = NULL) const;
    std::pair<MPSTensor<Matrix, SymmGroup>,
              MPSTensor<Matrix, SymmGroup> > split_mps_r2l(std::size_t Mmax, double cutoff, Logger * iter_log = NULL) const;
    
    void swap_with(TwoSiteTensor & b);

    friend void swap(TwoSiteTensor & a, TwoSiteTensor & b)
    {
        a.swap_with(b);
    }
    
#ifdef HAVE_ALPS_HDF5
    void load(alps::hdf5::archive & ar);
    void save(alps::hdf5::archive & ar) const;
#endif
    
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

#include "dmrg/mp_tensors/mpotensor.h"
template<class Matrix, class SymmGroup>
MPOTensor<Matrix, SymmGroup> make_twosite_mpo(MPOTensor<Matrix, SymmGroup> const & mpo1,
                                              MPOTensor<Matrix, SymmGroup> const & mpo2, Index<SymmGroup> const & phys_i)
{
    assert(mpo1.col_dim() == mpo2.row_dim());
    MPOTensor<Matrix, SymmGroup> mpo_big(mpo1.row_dim(), mpo2.col_dim());

    std::size_t b1, b2, b3;
    for( b1=0; b1 < mpo1.row_dim(); ++b1)
	for( b2=0; b2 < mpo1.col_dim(); ++b2)
	    for( b3=0; b3 < mpo2.col_dim(); ++b3) {
		if (! (mpo1.has(b1, b2) && mpo2.has(b2, b3)) )
		    continue;
		    block_matrix<Matrix, SymmGroup> tmp;
		    op_kron(phys_i, mpo1(b1,b2), mpo2(b2,b3), tmp);
		if (mpo_big.has(b1,b3))
		    mpo_big(b1,b3) += tmp;
		else
		    mpo_big(b1,b3) = tmp;
	    }

    return mpo_big;
}

template<class Matrix, class SymmGroup>
MPO<Matrix, SymmGroup> make_ts_cache_mpo(MPO<Matrix, SymmGroup> const & mpo_orig, Index<SymmGroup> const & site_dim)
{
    std::size_t L_ts = mpo_orig.length() - 1;
    MPO<Matrix, SymmGroup> r(L_ts);
    for(std::size_t p=0; p < L_ts; ++p)
    {
        r[p] = make_twosite_mpo(mpo_orig[p], mpo_orig[p+1], site_dim);
    }

    return r;

}

#endif
