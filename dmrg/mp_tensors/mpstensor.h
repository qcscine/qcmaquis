/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MPSTENSOR_H
#define MPSTENSOR_H

#include "block_matrix/block_matrix.h"
#include "block_matrix/indexing.h"

//#include "dense_matrix/matrix_algorithms.hpp"

#include <iostream>
#include <algorithm>

enum MPSStorageLayout { LeftPaired, RightPaired };
// these are actually used in several places
enum Indicator { Unorm, Lnorm, Rnorm };
enum DecompMethod { QR, SVD };

template<class Matrix, class SymmGroup>
class MPSTensor
{
public:
    typedef typename Matrix::value_type scalar_type;
    typedef typename Matrix::value_type value_type;
    typedef double real_type;
    typedef double magnitude_type;
    typedef std::size_t size_type;
    
    MPSTensor(Index<SymmGroup> const & sd = Index<SymmGroup>(),
              Index<SymmGroup> const & ld = Index<SymmGroup>(),
              Index<SymmGroup> const & rd = Index<SymmGroup>(),
              bool fillrand = true);
        
    Index<SymmGroup> const & site_dim() const;
    Index<SymmGroup> const & row_dim() const;
    Index<SymmGroup> const & col_dim() const;
    bool isobccompatible(Indicator) const;
    
    // these are not const because after a numerical test
    // they may update the status
    bool isleftnormalized(bool test = false);
    bool isrightnormalized(bool test = false);
    bool isnormalized(bool test = false);
    
    block_matrix<Matrix, SymmGroup> normalize_left(DecompMethod method = QR,
                                                   bool multiplied = true,
                                                   double truncation = 0,
                                                   Index<SymmGroup> bond_dim = Index<SymmGroup>());
    block_matrix<Matrix, SymmGroup> normalize_right(DecompMethod method = QR,
                                                    bool multiplied = true,
                                                    double truncation = 0,
                                                    Index<SymmGroup> bond_dim = Index<SymmGroup>());
    
    void multiply_from_left(block_matrix<Matrix, SymmGroup> const &);
    void multiply_from_right(block_matrix<Matrix, SymmGroup> const &);
    void multiply_by_scalar(scalar_type);
    
    scalar_type scalar_overlap(MPSTensor const &) const;
    real_type scalar_norm() const;
    
    // this is completely useless in C++, only exists for consistency with Python
    MPSTensor copy() const;
    
    block_matrix<Matrix, SymmGroup> & data();
    block_matrix<Matrix, SymmGroup> const & data() const;
    
    std::vector<block_matrix<Matrix, SymmGroup> > to_list() const;
    
    template<class Matrix_, class SymmGroup_>
    friend std::ostream& operator<<(std::ostream&, MPSTensor<Matrix_, SymmGroup_> const &);
    
    friend struct contraction;
    friend struct compression;
    friend struct multigrid;
    
    // math functions: these are not part of the Python code, but required by IETL
    MPSTensor const & operator*=(scalar_type);
    MPSTensor const & operator/=(scalar_type);
    
    MPSTensor const & operator+=(MPSTensor const &);
    MPSTensor const & operator-=(MPSTensor const &);
    
    void make_left_paired() const;
    void make_right_paired() const;
    
    void swap_with(MPSTensor & b);
    friend void swap(MPSTensor & a, MPSTensor & b)
    {
        a.swap_with(b);
    }
    
    void inplace_conjugate();
    
    /* friend
    MPSTensor join(MPSTensor const & m1,
                   MPSTensor const & m2)
    {
        m1.make_left_paired();
        m2.make_left_paired();
        
        MPSTensor<Matrix, SymmGroup> ret;
        
        for (std::size_t b = 0; b < m1.data_.n_blocks(); ++b) {
            if (m2.data_.has_block(m1.data_.left_basis()[b], m1.data_.right_basis()[b])) {
                Matrix nb = blas::join(m1.data_[b],
                                       m2.data_(m1.data_.left_basis()[b].first, m1.data_.right_basis()[b].first));
                
                ret.data_.insert_block(nb, m1.data_.left_basis()[b].first, m1.data_.right_basis()[b].first);
            } else
                ret.data_.insert_block(m1.data_[b], m1.data_.left_basis()[b].first, m1.data_.right_basis()[b].first);
        }
        
        for (std::size_t b = 0; b < m2.data_.n_blocks(); ++b) {
            if (m1.data_.has_block(m2.data_.left_basis()[b], m2.data_.right_basis()[b])) // those should've been taken care of above.
                continue;
                
            ret.data_.insert_block(m2.data_[b], m2.data_.left_basis()[b].first, m2.data_.right_basis()[b].first);
        }
        
        ret.right_i = ret.data_.right_basis();
        
        ret.left_i = m1.left_i;
        for (typename Index<SymmGroup>::const_iterator it = m2.left_i.begin();
             it != m2.left_i.end(); ++it) {
            if (ret.left_i.has(it->first))
                ret.left_i[ret.left_i.position(it->first)].second += it->second;
            else
                ret.left_i.insert(*it);
        }
        
        ret.phys_i = m1.phys_i;
        
        return ret;
    } */
    
#ifdef HAVE_ALPS_HDF5
    void load(alps::hdf5::archive & ar);
    void save(alps::hdf5::archive & ar) const;
#endif
    
    void check_equal(MPSTensor<Matrix, SymmGroup> const &) const;
    bool reasonable() const;
    
private:
    Index<SymmGroup> phys_i, left_i, right_i;
    mutable block_matrix<Matrix, SymmGroup> data_;
    mutable MPSStorageLayout cur_storage;
    Indicator cur_normalization;
};

// this is also required by IETL
template<class Matrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup> operator*(typename MPSTensor<Matrix, SymmGroup>::scalar_type t,
                                       MPSTensor<Matrix, SymmGroup> m)
{
    m *= t;
    return m;
}
template<class Matrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup> operator*(MPSTensor<Matrix, SymmGroup> m,
                                       typename MPSTensor<Matrix, SymmGroup>::scalar_type t)
{
    m *= t;
    return m;
}
template<class Matrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup> operator/(MPSTensor<Matrix, SymmGroup> m,
                                       typename MPSTensor<Matrix, SymmGroup>::scalar_type t)
{
    m *= 1/t;
    return m;
}

template<class Matrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup> operator-(MPSTensor<Matrix, SymmGroup> m,
                                       MPSTensor<Matrix, SymmGroup> const & m2)
{
    m -= m2;
    return m;
}
template<class Matrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup> operator+(MPSTensor<Matrix, SymmGroup> m,
                                       MPSTensor<Matrix, SymmGroup> const & m2)
{
    m += m2;
    return m;
}

template<class Matrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup> operator-(MPSTensor<Matrix, SymmGroup> m)
{
    m *= -1;
    return m;
}

#include "mp_tensors/mpstensor.hpp"

#endif
