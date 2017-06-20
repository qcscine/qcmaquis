/*****************************************************************************
 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations
 *
 * ALPS Libraries
 *
 * Copyright (C) 2017-2017 by Alberto Baiardi <alberto.baiardi@sns.it>
 *
 * This software is part of the ALPS libraries, published under the ALPS
 * Library License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 *
 * You should have received a copy of the ALPS Library License along with
 * the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
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

#ifndef PARTIAL_OVERLAP_H
#define PARTIAL_OVERLAP_H

#include <sstream>
#include <algorithm>
#include <numeric>

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/block_matrix/symmetry.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpstensor.h"

#include <ietl/traits.h>
#include <ietl/fmatrix.h>
#include <ietl/ietl2lapack.h>

#include "utils/traits.hpp"

#include <boost/ptr_container/ptr_vector.hpp>

template<class Matrix, class SymmGroup>
class partial_overlap
{
public:
    // Types definition
    typedef typename boost::ptr_vector<Matrix> vector_overlap ;
    typedef typename MPS<Matrix,SymmGroup>::MPS MPSWave ;
    typedef typename MPSTensor<Matrix,SymmGroup>::MPSTensor MPSTensor ;
    typedef typename block_matrix<Matrix,SymmGroup>::block_matrix bmatrix ;
    typedef int basis_type ;
    typedef typename std::vector<basis_type> basis_vector ;
    typedef typename std::size_t dim_type ;
    typedef typename SymmGroup::charge charge ;
    typedef typename Matrix::value_type value_type ;
    // Constructors
    partial_overlap(const dim_type& L , const dim_type& m );
    partial_overlap(const dim_type& L , const dim_type& m , const basis_vector& basis);
    partial_overlap(const MPSWave& MPS , const basis_vector& basis);
    // Methods
    basis_type get_basis(const dim_type& site) ;
    void update(const MPSWave& MPS, const dim_type& l, const int& direction) ;
    value_type overlap(const dim_type& i);
private:
    // Private attributes
    vector_overlap data_left_ , data_right_ ;
    dim_type lattice_L_ ;
    basis_vector basis_ ;
    // Local variables
    enum Modality { Left , Right };
    static const charge identity = SymmGroup::IdentityCharge ;
    // Private functions
    void multiply (const MPSWave& MPS, const Matrix& B, const basis_type& sigma,
                   const dim_type& l, const dim_type& m1, const dim_type& m2,
                   const Modality& mod, Matrix& output) ;
    void multiply_first (const MPSWave& MPS, const basis_type& sigma, const dim_type& l,
                         const dim_type& m1, const dim_type& m2, const Modality& mod,
                         Matrix& output) ;
    void extract(const bmatrix& bm, const basis_type& sigma, const std::size_t m1, const std::size_t m2, Matrix& output) ;
};

// +------------+
//  CONSTRUCTORS
// +------------+

template<class Matrix, class SymmGroup>
partial_overlap<Matrix,SymmGroup>::partial_overlap(const partial_overlap<Matrix,SymmGroup>::dim_type& L ,
                                                   const partial_overlap<Matrix,SymmGroup>::dim_type& m)
                                                   : lattice_L_(L)
{
    for (dim_type k = 0; k < L; ++k) {
        data_left_.push_back(new Matrix(m,m)) ;
        data_right_.push_back(new Matrix(m,m)) ;
        basis_.push_back(0);
    }

};

template<class Matrix, class SymmGroup>
partial_overlap<Matrix,SymmGroup>::partial_overlap(const partial_overlap<Matrix,SymmGroup>::dim_type& L ,
                                                   const partial_overlap<Matrix,SymmGroup>::dim_type& m ,
                                                   const partial_overlap<Matrix,SymmGroup>::basis_vector& basis): lattice_L_(L)
{
    assert( L == basis.size()) ;
    for (dim_type k = 0; k < L; ++k) {
        data_left_.push_back(new Matrix(m, m)) ;
        data_right_.push_back(new Matrix(m,m)) ;
        basis_.push_back(basis[k]);
    }

};

template<class Matrix, class SymmGroup>
partial_overlap<Matrix,SymmGroup>::partial_overlap(const partial_overlap<Matrix,SymmGroup>::MPSWave& MPS ,
                                                   const partial_overlap<Matrix,SymmGroup>::basis_vector& basis)
{
    dim_type L = MPS.size() ;
    if (basis.size() == 0)
        for (int i = 0 ; i < L ; ++i)
            basis_.push_back(0) ;
    else
        for (int i = 0 ; i < L ; ++i)
            basis_.push_back(basis[i]) ;
    lattice_L_ = L ;
    dim_type m1, m2 ;
    MPSTensor MPSTns ;
    for (dim_type i = 0; i < L; ++i) {
        Matrix *tmp ;
        if (i == 0) {
            // Left Update
            MPSTns = MPS[i];
            MPSTns.make_left_paired();
            m1 = MPSTns.row_dim().size_of_block(identity);
            m2 = MPSTns.col_dim().size_of_block(identity);
            tmp = new Matrix(m2, m2);
            multiply_first(MPS, basis_[i], i, m1, m2, Left, *tmp);
            data_left_.push_back(tmp);
            // Right update
            MPSTns = MPS[L-1-i];
            MPSTns.make_left_paired();
            m1 = MPSTns.row_dim().size_of_block(identity);
            m2 = MPSTns.col_dim().size_of_block(identity);
            tmp = new Matrix(m1, m1);
            multiply_first(MPS, basis_[L-1], L-1, m1, m2, Right, *tmp);
            data_right_.push_back(tmp);
        } else {
            // Left update
            MPSTns = MPS[i] ;
            MPSTns.make_left_paired() ;
            m1 = MPSTns.row_dim().size_of_block(identity) ;
            m2 = MPSTns.col_dim().size_of_block(identity) ;
            tmp = new Matrix(m1,m2) ;
            multiply(MPS, data_left_[i-1], basis_[i], i, m1, m2, Left, *tmp) ;
            if (i != lattice_L_-1) {
                data_left_.push_back(tmp);
            } else {
                Matrix* tmp2;
                tmp2 = new Matrix(1,1);
                (*tmp2)(0,0) = (*tmp)(0,0);
                data_left_.push_back(tmp2);
            }
            // Right update
            MPSTns = MPS[L-1-i] ;
            MPSTns.make_left_paired() ;
            m1 = MPSTns.row_dim().size_of_block(identity) ;
            m2 = MPSTns.col_dim().size_of_block(identity) ;
            tmp = new Matrix(m1,m2) ;
            multiply(MPS, data_right_[i-1], basis_[L-1-i], L-1-i, m1, m2, Right, *tmp) ;
            if (i != lattice_L_-1){
                data_right_.push_back(tmp);
            } else {
                Matrix* tmp2 ;
                tmp2 = new Matrix(1,1);
                (*tmp2)(0,0) = (*tmp)(0,0) ;
                data_right_.push_back(tmp2);
            }
        }
    }
};

// +-------+
//  METHODS
// +-------+

template<class Matrix, class SymmGroup>
typename partial_overlap<Matrix, SymmGroup>::basis_type partial_overlap<Matrix, SymmGroup>::get_basis(const dim_type &site)
{
    return basis_[site] ;
}

template<class Matrix, class SymmGroup>
void partial_overlap<Matrix, SymmGroup>::update(const partial_overlap<Matrix, SymmGroup>::MPSWave& MPS,
                                                const dim_type& l,
                                                const int& direction)
{
    // Get the index to update in the left and right overlaps
    dim_type i1, i2 ;
    Modality mod ;
    if (direction == 1)
        mod = Left;
    else if (direction == -1)
        mod = Right;
    // Set the sites where the
    if (mod == Left) {
        i1 = l;
        i2 = l+1 ;
    } else if (mod == Right) {
        i1 = lattice_L_-1 - l ;
        i2 = lattice_L_-1 - (l+1) ;
    }
    MPSTensor MPSTns ;
    dim_type m1, m2 ;
    //
    MPSTns = MPS[l] ;
    MPSTns.make_left_paired() ;
    m1 = MPSTns.row_dim().size_of_block(identity);
    m2 = MPSTns.col_dim().size_of_block(identity);
    std::cout << m1 << " " << m2 << std::endl ;
    Matrix *tmp ;
    if (i1 == 0) {
        if (mod == Left) {
            tmp = new Matrix(m2, m2);
            multiply_first(MPS, basis_[0], 0, m1, m2, mod, *tmp);
            data_left_[0] = *tmp;
        } else if (mod == Right) {
            tmp = new Matrix(m1, m1);
            multiply_first(MPS, basis_[l], l, m1, m2, mod, *tmp);
            data_right_[0] = *tmp;
        }
    } else {
        if (mod == Left) {
            tmp = new Matrix(m1, m2);
            multiply(MPS, data_left_[i1-1], basis_[l], l, m1, m2, mod, *tmp);
            if (i1 != lattice_L_-1 ) {
                data_left_[i1] = *tmp;
            } else {
                Matrix* tmp2 ;
                tmp2 = new Matrix(1,1);
                (*tmp2)(0,0) = (*tmp)(0,0) ;
                data_left_[i1] = *tmp2 ;
            }
        } else if (mod == Right) {
            tmp = new Matrix(m1, m2);
            multiply(MPS, data_right_[i1-1], basis_[l], l, m1, m2, mod, *tmp);
            if (i1 != lattice_L_-1 ) {
                data_right_[i1] = *tmp;
            } else {
                Matrix* tmp2 ;
                tmp2 = new Matrix(1,1);
                (*tmp2)(0,0) = (*tmp)(0,0) ;
                data_right_[i1] = *tmp2 ;
            }
        }
    }
};

template<class Matrix, class SymmGroup>
typename partial_overlap<Matrix, SymmGroup>::value_type partial_overlap<Matrix, SymmGroup>::overlap(const dim_type &i)
{
    // Check data consistency
    assert (i < lattice_L_) ;
    value_type result = 0 ;
    if (i == 0) {
        Matrix a = data_left_[0];
        assert (a.num_cols() == a.num_rows()) ;
        for (int k = 0; k < a.num_rows(); ++k)
            result += a(k, k);
    } else if (i == lattice_L_ ) {
        Matrix a = data_right_[lattice_L_-1] ;
        assert (a.num_cols() == a.num_rows()) ;
        for (int k = 0; k < a.num_rows() ; ++k)
            result += a(k,k) ;
    } else {
        Matrix a = data_left_[i] ;
        Matrix b = data_right_[lattice_L_-i] ;
        assert (a.num_rows() == b.num_rows() && a.num_cols() == b.num_cols()) ;
        for (int k1 = 0; k1 < a.num_rows() ; ++k1)
            for (int k2 = 0; k2 < a.num_cols() ; ++k2)
                result += a(k1,k2)*b(k1,k2) ;
    }
    return result ;
};

// +---------------+
//  PRIVATE METHODS
// +---------------+

template<class Matrix, class SymmGroup>
void partial_overlap<Matrix, SymmGroup>::multiply(const partial_overlap<Matrix,SymmGroup>::MPSWave& MPS,
                                                  const Matrix& B,
                                                  const basis_type& sigma,
                                                  const dim_type& l ,
                                                  const dim_type& m1 ,
                                                  const dim_type& m2,
                                                  const partial_overlap<Matrix, SymmGroup>::Modality& mod,
                                                  Matrix& output)
{
    // Initialization
    MPSTensor MPSTns = MPS[l] ;
    MPSTns.make_left_paired() ;
    block_matrix<Matrix, SymmGroup> bm = MPSTns.data() ;
    Matrix *tmp ;
    tmp = new Matrix(m1, m2) ;
    extract(bm, sigma, m1, m2, *tmp) ;
    if (mod == Left)
        gemm(transpose(B), *tmp, output);
    else if (mod == Right)
        gemm(*tmp, transpose(B), output);
}

template<class Matrix, class SymmGroup>
void partial_overlap<Matrix, SymmGroup>::multiply_first(const partial_overlap<Matrix,SymmGroup>::MPSWave& MPS,
                                                        const basis_type& sigma,
                                                        const dim_type& l ,
                                                        const dim_type& m1,
                                                        const dim_type& m2,
                                                        const partial_overlap<Matrix, SymmGroup>::Modality& mod,
                                                        Matrix& output)
{
    MPSTensor MPSTns = MPS[l] ;
    MPSTns.make_left_paired() ;
    block_matrix<Matrix, SymmGroup> bm = MPSTns.data() ;
    Matrix *tmp ;
    if (mod == Left) {
        tmp = new Matrix(1,m2) ;
        extract(bm, sigma, m1, m2, *tmp) ;
        assert (m1 == 1) ;
        for (int i = 0; i < m2; ++i)
            output(i,i) = (*tmp)(0,i) ;
    } else if (mod == Right) {
        tmp = new Matrix(m1,1) ;
        extract(bm, sigma, m1, m2, *tmp) ;
        assert (m2 == 1) ;
        for (int i = 0; i < m1; ++i)
            output(i,i) = (*tmp)(i,0) ;
    }
    return ;
}

template<class Matrix, class SymmGroup>
void partial_overlap<Matrix, SymmGroup>::extract(const bmatrix& bm,
                                                 const basis_type& sigma,
                                                 const std::size_t m1,
                                                 const std::size_t m2,
                                                 Matrix& output)
{
    std::size_t offset = m1*sigma ;
    for (int i = 0; i < m1; ++i)
        for (int j = 0 ; j < m2 ; ++j)
            output(i,j) = bm[0](i+offset,j);
};

#endif