/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MPS_H
#define MPS_H

#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/mpotensor.h"
#include "dmrg/mp_tensors/boundary.h"

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
    void canonize(size_t center, DecompMethod method = DefaultSolver());
    
    void normalize_left();
    void normalize_right();
    
    void move_normalization_l2r(size_t p1, size_t p2, DecompMethod method=DefaultSolver());
    void move_normalization_r2l(size_t p1, size_t p2, DecompMethod method=DefaultSolver());
    
    std::string description() const;
   
    template<class OtherMatrix>
    truncation_results grow_l2r_sweep(MPOTensor<Matrix, SymmGroup> const & mpo,
                                      Boundary<OtherMatrix, SymmGroup> const & left,
                                      Boundary<OtherMatrix, SymmGroup> const & right,
                                      std::size_t l, double alpha,
                                      double cutoff, std::size_t Mmax);
    template<class OtherMatrix>
    truncation_results grow_r2l_sweep(MPOTensor<Matrix, SymmGroup> const & mpo,
                                      Boundary<OtherMatrix, SymmGroup> const & left,
                                      Boundary<OtherMatrix, SymmGroup> const & right,
                                      std::size_t l, double alpha,
                                      double cutoff, std::size_t Mmax);
    
    Boundary<Matrix, SymmGroup> left_boundary() const;
    Boundary<Matrix, SymmGroup> right_boundary() const;
    
    void apply(block_matrix<Matrix, SymmGroup> const&, size_type);
    void apply(block_matrix<Matrix, SymmGroup> const&, block_matrix<Matrix, SymmGroup> const&, size_type);
    
    friend void swap(MPS& a, MPS& b)
    {
        using std::swap;
        swap(a.data_, b.data_);
        swap(a.canonized_i, b.canonized_i);
    }
    
private:
    
    data_t data_;
    mutable size_t canonized_i;
};

template<class Matrix, class SymmGroup>
void load(std::string const& dirname, MPS<Matrix, SymmGroup> & mps);
template<class Matrix, class SymmGroup>
void save(std::string const& dirname, MPS<Matrix, SymmGroup> const& mps);

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
                            MPS<Matrix, SymmGroup> const & b,
                            double alpha=1., double beta=1.)
{
    assert( a.length() == b.length() );
    
    MPSTensor<Matrix, SymmGroup> aright=a[a.length()-1], bright=b[a.length()-1];
    aright.multiply_by_scalar(alpha);
    bright.multiply_by_scalar(beta);

    MPS<Matrix, SymmGroup> ret(a.length());
    ret[0] = join(a[0],b[0],l_boundary_f);
    ret[a.length()-1] = join(aright,bright,r_boundary_f);
    for (std::size_t p = 1; p < a.length()-1; ++p)
        ret[p] = join(a[p], b[p]);
    return ret;
}

#include "dmrg/mp_tensors/mps.hpp"

#endif
