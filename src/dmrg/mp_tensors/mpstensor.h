#ifndef MPSTENSOR_H
#define MPSTENSOR_H

#include "p_block_matrix/p_block_matrix.h"
#include "p_block_matrix/indexing.h"

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
    
    p_block_matrix<Matrix, SymmGroup> normalize_left(DecompMethod method = QR,
                                                   bool multiplied = true,
                                                   double truncation = 0,
                                                   Index<SymmGroup> bond_dim = Index<SymmGroup>());
    p_block_matrix<Matrix, SymmGroup> normalize_right(DecompMethod method = QR,
                                                    bool multiplied = true,
                                                    double truncation = 0,
                                                    Index<SymmGroup> bond_dim = Index<SymmGroup>());
    
    void multiply_from_left(p_block_matrix<Matrix, SymmGroup> const &);
    void multiply_from_right(p_block_matrix<Matrix, SymmGroup> const &);
    void multiply_by_scalar(scalar_type);
    
    scalar_type scalar_overlap(MPSTensor const &) const;
    real_type scalar_norm() const;
    
    // this is completely useless in C++, only exists for consistency with Python
    MPSTensor copy() const;
    
    p_block_matrix<Matrix, SymmGroup> & data();
    p_block_matrix<Matrix, SymmGroup> const & data() const;
    
    template<class Matrix_, class SymmGroup_>
    friend std::ostream& operator<<(std::ostream&, MPSTensor<Matrix_, SymmGroup_> const &);
    
    friend struct contraction;
    
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
    
#ifdef HAVE_ALPS_HDF5
    void serialize(alps::hdf5::iarchive & ar);
    void serialize(alps::hdf5::oarchive & ar) const;
#endif
    
private:
    Index<SymmGroup> phys_i, left_i, right_i;
    mutable p_block_matrix<Matrix, SymmGroup> data_;
    mutable MPSStorageLayout cur_storage;
    Indicator cur_normalization;
};

// this is also required by IETL
template<class Matrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup> operator*(typename MPSTensor<Matrix, SymmGroup>::scalar_type t, MPSTensor<Matrix, SymmGroup> m)
{
    m *= t;
    return m;
}

#include "mp_tensors/mpstensor.hpp"

#endif
