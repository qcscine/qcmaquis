/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 * 
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 * 
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#ifndef MPSTENSOR_H
#define MPSTENSOR_H

#include <iostream>
#include <algorithm>

#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/indexing.h"
//#include "solver.h"

#ifdef ENABLE_SU2
#include "dmrg/mp_tensors/contractions/non-abelian/lq_scale.hpp"
#endif

enum boundary_flag_t {no_boundary_f,l_boundary_f,r_boundary_f};
enum MPSStorageLayout { LeftPaired, RightPaired };
// these are actually used in several places
enum Indicator { Unorm, Lnorm, Rnorm };
enum DecompMethod {QR, SVD}; 

static DecompMethod DefaultSolver() {return QR;} // QR or SVD

template<class Matrix, class SymmGroup>
class TwoSiteTensor;

template<class Matrix, class SymmGroup>
class MPSTensor
{
public:
    typedef typename maquis::traits::scalar_type<Matrix>::type scalar_type;
    typedef typename maquis::traits::real_type<Matrix>::type real_type;
    typedef typename Matrix::value_type value_type;
    typedef double magnitude_type; // should become future (todo: Matthias, 30.04.12 / scalar-value types)
    typedef std::size_t size_type;
    
    MPSTensor(Index<SymmGroup> const & sd = Index<SymmGroup>(),
              Index<SymmGroup> const & ld = Index<SymmGroup>(),
              Index<SymmGroup> const & rd = Index<SymmGroup>(),
              bool fillrand = true,
              value_type val = 0);

    MPSTensor(Index<SymmGroup> const& sd,
              Index<SymmGroup> const& ld,
              Index<SymmGroup> const& rd,
              block_matrix<Matrix, SymmGroup> const& block,
              MPSStorageLayout layout);

    Index<SymmGroup> const & site_dim() const;
    Index<SymmGroup> const & row_dim() const;
    Index<SymmGroup> const & col_dim() const;
    bool isobccompatible(Indicator) const;
    std::size_t num_elements() const;
    #ifdef USE_AMBIENT
    void print_distribution() const;
    #endif
    
    void replace_right_paired(block_matrix<Matrix, SymmGroup> const &, Indicator =Unorm);
    void replace_left_paired(block_matrix<Matrix, SymmGroup> const &, Indicator =Unorm);
    
    // these are not const because after a numerical test
    // they may update the status
    bool isleftnormalized(bool test = false) const;
    bool isrightnormalized(bool test = false) const;
    bool isnormalized(bool test = false) const;
    
    block_matrix<Matrix, SymmGroup> normalize_left(DecompMethod method = DefaultSolver(),
                                                   bool multiplied = true,
                                                   double truncation = 0,
                                                   Index<SymmGroup> bond_dim = Index<SymmGroup>());
    block_matrix<Matrix, SymmGroup> normalize_right(DecompMethod method = DefaultSolver(),
                                                    bool multiplied = true,
                                                    double truncation = 0,
                                                    Index<SymmGroup> bond_dim = Index<SymmGroup>());
    
    void shift_aux_charges(typename SymmGroup::charge);
    
    template <class OtherMatrix>
    void multiply_from_left(block_matrix<OtherMatrix, SymmGroup> const &);
    template <class OtherMatrix>
    void multiply_from_right(block_matrix<OtherMatrix, SymmGroup> const &);
    void multiply_by_scalar(const scalar_type&);
    void divide_by_scalar(const scalar_type&);
    
    scalar_type scalar_overlap(MPSTensor const &) const;
    real_type scalar_norm() const;
    
    // this is completely useless in C++, only exists for consistency with Python
    MPSTensor copy() const;
    
    block_matrix<Matrix, SymmGroup> & data();
    block_matrix<Matrix, SymmGroup> const & data() const;
    block_matrix<Matrix, SymmGroup> const & const_data() const;
    
    std::vector<block_matrix<Matrix, SymmGroup> > to_list() const;
    
    template<class Matrix_, class SymmGroup_>
    friend std::ostream& operator<<(std::ostream&, MPSTensor<Matrix_, SymmGroup_> const &);
    
    // math functions: these are not part of the Python code, but required by IETL
    MPSTensor const & operator*=(const scalar_type&);
    MPSTensor const & operator/=(const scalar_type&);
    
    MPSTensor const & operator+=(MPSTensor const &);
    MPSTensor const & operator-=(MPSTensor const &);
    
    void make_left_paired() const;
    void make_right_paired() const;
    
    void conjugate_inplace();
    void swap_with(MPSTensor & b);
    friend void swap(MPSTensor& a, MPSTensor& b){
        a.swap_with(b); 
    }
    template<class Matrix_, class SymmGroup_>
    friend MPSTensor<Matrix_, SymmGroup_> join(MPSTensor<Matrix_, SymmGroup_> const &, MPSTensor<Matrix_, SymmGroup_> const &, boundary_flag_t);
    
    template<class Archive> void load(Archive & ar);
    template<class Archive> void save(Archive & ar) const;
    template <class Archive> void serialize(Archive & ar, const unsigned int version);

    void check_equal(MPSTensor<Matrix, SymmGroup> const &) const;
    bool reasonable() const;
    bool num_check() const; // checks for nan or inf
    
    Index<SymmGroup> phys_i, left_i, right_i;
private:
    mutable block_matrix<Matrix, SymmGroup> data_;
    mutable MPSStorageLayout cur_storage;
    Indicator cur_normalization;
};

// this is also required by IETL
template<class Matrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup> operator*(const typename MPSTensor<Matrix, SymmGroup>::scalar_type& t,
                                       MPSTensor<Matrix, SymmGroup> m)
{
    m *= t;
    return m;
}
template<class Matrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup> operator*(MPSTensor<Matrix, SymmGroup> m,
                                       const typename MPSTensor<Matrix, SymmGroup>::scalar_type& t)
{
    m *= t;
    return m;
}
template<class Matrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup> operator/(MPSTensor<Matrix, SymmGroup> m,
                                       const typename MPSTensor<Matrix, SymmGroup>::scalar_type& t)
{
    m /= t;
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
    m *= typename MPSTensor<Matrix, SymmGroup>::scalar_type(-1.0);
    return m;
}


template<class Matrix, class SymmGroup>
std::size_t size_of(MPSTensor<Matrix, SymmGroup> const & m)
{
    return size_of(m.data());
}

#include "dmrg/mp_tensors/mpstensor.hpp"

#endif
