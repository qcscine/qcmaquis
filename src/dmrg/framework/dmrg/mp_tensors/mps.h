/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MPS_H
#define MPS_H

#include "dmrg/utils/logger.h"

#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/mpotensor.h"

#include <limits>

template<class Matrix, class SymmGroup>
struct mps_initializer;

template<class Matrix, class SymmGroup>
class MPS
{
    typedef std::vector<MPSTensor<Matrix, SymmGroup> > data_t;
public:
    typedef std::size_t size_t;

    // reproducing interface of std::vector
    typedef typename data_t::size_type size_type;
    typedef typename data_t::value_type value_type;
    typedef typename data_t::iterator iterator;
    typedef typename data_t::const_iterator const_iterator;
    typedef typename MPSTensor<Matrix, SymmGroup>::scalar_type scalar_type;
    
    MPS();
    MPS(size_t L);  
    MPS(size_t L, size_t Mmax, Index<SymmGroup> phys,
        typename SymmGroup::charge right_end,
        mps_initializer<Matrix, SymmGroup> & init);
    
    size_t size() const { return data_.size(); }
    size_t length() const { return size(); }
    Index<SymmGroup> const & site_dim(size_t i) const { return data_[i].site_dim(); }
    Index<SymmGroup> const & row_dim(size_t i) const { return data_[i].row_dim(); }
    Index<SymmGroup> const & col_dim(size_t i) const { return data_[i].col_dim(); }
    
    value_type const & operator[](size_t i) const;
    value_type& operator[](size_t i);
    
    void resize(size_t L);
    
    const_iterator begin() const {return data_.begin();}
    const_iterator end() const {return data_.end();}
    const_iterator const_begin() const {return data_.begin();}
    const_iterator const_end() const {return data_.end();}
    iterator begin() {return data_.begin();}
    iterator end() {return data_.end();}
    
    size_t canonization(bool=false) const;
    void canonize(size_t center);
    
    void normalize_left();
    void normalize_right();
    
    void move_normalization_l2r(size_t p1, size_t p2);
    void move_normalization_r2l(size_t p1, size_t p2);
    
    std::string description() const;
    
    void grow_l2r_sweep(MPOTensor<Matrix, SymmGroup> const & mpo,
                        Boundary<Matrix, SymmGroup> const & left,
                        Boundary<Matrix, SymmGroup> const & right,
                        std::size_t l, double alpha,
                        double cutoff, std::size_t Mmax,
                        Logger &);
    void grow_r2l_sweep(MPOTensor<Matrix, SymmGroup> const & mpo,
                        Boundary<Matrix, SymmGroup> const & left,
                        Boundary<Matrix, SymmGroup> const & right,
                        std::size_t l, double alpha,
                        double cutoff, std::size_t Mmax,
                        Logger &);
    
    Boundary<Matrix, SymmGroup> left_boundary() const;
    Boundary<Matrix, SymmGroup> right_boundary() const;
    
#ifdef HAVE_ALPS_HDF5
    void load(alps::hdf5::archive & ar);
    void save(alps::hdf5::archive & ar) const;
#endif
    
private:
    
    data_t data_;
    mutable size_t canonized_i;
};

template<class Matrix, class SymmGroup>
struct mps_initializer
{
    virtual void operator()(MPS<Matrix, SymmGroup> & mps,
                            std::size_t Mmax,
                            Index<SymmGroup> const & phys,
                            typename SymmGroup::charge right_end) = 0;
};

template<class Matrix, class SymmGroup>
MPS<Matrix, SymmGroup> join(MPS<Matrix, SymmGroup> const & a,
                            MPS<Matrix, SymmGroup> const & b)
{
    assert( a.length() == b.length() );
    
    MPS<Matrix, SymmGroup> ret = a;
    for (std::size_t p = 0; p < a.length(); ++p)
        ret[p] = join(a[p], b[p]);
    return ret;
}

#include "dmrg/mp_tensors/mps.hpp"

#endif
