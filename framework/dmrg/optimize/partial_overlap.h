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
    // Constructors
    partial_overlap(const dim_type& L , const dim_type& m );
    partial_overlap(const dim_type& L , const dim_type& m , const basis_vector& basis);
    partial_overlap(const MPSWave& MPS );
    partial_overlap(const MPSWave& MPS , const basis_vector& basis);
    // Methods
    basis_type get_basis(const dim_type& site) ;
    void update_left(const MPSTensor& MPST, const dim_type& l) ;
    void update_right(const MPSTensor& MPST, const dim_type& r) ;
private:
    // Private attributes
    vector_overlap data_left_ , data_right_ ;
    dim_type lattice_L_ ;
    basis_vector basis_ ;
    // Local variables
    enum Modality { Left , Right };
    charge identity = SymmGroup::IdentityCharge ;
    // Private functions
    Matrix multiply (const MPSWave& MPS, const Matrix& B, const basis_type& sigma,
                     const dim_type& l, const Modality& mod) ;
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
                                                   const partial_overlap<Matrix,SymmGroup>::basis_vector& basis): lattice_l_(L)
{
    assert( L == basis.size()) ;
    for (dim_type k = 0; k < L; ++k) {
        data_left_.push_back(new Matrix(m, m)) ;
        data_right_.push_back(new Matrix(m,m)) ;
        basis_.push_back(basis[k]);
    }

};

template<class Matrix, class SymmGroup>
partial_overlap<Matrix,SymmGroup>::partial_overlap(const partial_overlap<Matrix,SymmGroup>::MPSWave& MPS)
{
    dim_type L = MPS.size() ;
    lattice_L_ = L ;
    dim_type m1, m2 ;
    for (dim_type i = 0; i < L; ++i) {
        basis_.push_back(0);
        // Left overlap
        Modality mod = Left;
        m1 = MPS[i].row_dim().size_of_block(identity) ;
        m2 = MPS[i].col_dim().size_of_block(identity) ;
        if (i == 0) {
            //Matrix *tmp ;
            //tmp = new Matrix(m2,m2);
            Matrix tmp(1,1) ;
            multiply_first(MPS, basis_[i], i, m1, m2, mod, tmp);
            //data_left_.push_back(tmp) ;
            std::cout << "Pippo" << std::endl ;
        } else {
            data_left_.push_back(new Matrix(m1, m2));
            data_left_[i] = multiply(MPS, data_left_[i-1], basis_[i], i, mod);
        }
    }
};

template<class Matrix, class SymmGroup>
partial_overlap<Matrix,SymmGroup>::partial_overlap(const partial_overlap<Matrix,SymmGroup>::MPSWave& MPS ,
                                                   const partial_overlap<Matrix,SymmGroup>::basis_vector& basis)
{
    dim_type L = MPS.size() ;
    assert(L == basis.size) ;
    dim_type m1, m2 ;
    for (dim_type i = 0; i < L; ++i){
        m1 = MPS.row_dim() ;
        m2 = MPS.col_dim() ;
        // ALB TODO here check if the matrix is square. Don't know if this check is general enough
        assert(m1 == m2);
        data_.push_back(new Matrix(m1, m1));
        basis_push_back(basis[i]);
        if (i == 0)
            *(data_[i]) = MPS[0].data()[basis_[i]] ;
        else
            *(data_[i]) = multiply(MPS, *(data_[i-1]), basis_[i], i);
    }
    lattice_L_ = L ;
};

// +-------+
//  METHODS
// +-------+

template<class Matrix, class SymmGroup>
partial_overlap<Matrix, SymmGroup>::basis_type partial_overlap<Matrix, SymmGroup>::get_basis(const dim_type &site)
{
    return basis_[site] ;
}

template<class Matrix, class SymmGroup>
void partial_overlap<Matrix, SymmGroup>::update_left(const partial_overlap<Matrix, SymmGroup>::MPSTensor& MPST,
                                                     const dim_type& l)
{
    basis_type sigma = basis_[l] ;
    Matrix result, tmp = MPST.data()[sigma] ;
    ietl::mult(tmp, &(data_left_[l]), result) ;
    data_left_[l+1] = *result ;
};

template<class Matrix, class SymmGroup>
void partial_overlap<Matrix, SymmGroup>::update_right(const partial_overlap<Matrix, SymmGroup>::MPSTensor &MPST,
                                                      const dim_type &r)
{
    // ALB TODO Write this part of code
};

// +---------------+
//  PRIVATE METHODS
// +---------------+

template<class Matrix, class SymmGroup>
Matrix partial_overlap<Matrix, SymmGroup>::multiply(const partial_overlap<Matrix,SymmGroup>::MPSWave& MPS,
                                                    const Matrix& B,
                                                    const basis_type& sigma,
                                                    const dim_type& l ,
                                                    const partial_overlap<Matrix, SymmGroup>::Modality& mod)
{
    // Initialization
    block_matrix<Matrix, SymmGroup> bm = MPS[l].data() ;
    std::size_t m1 = MPS[l].row_dim().size_of_block(identity) ;
    std::size_t m2 = MPS[l].col_dim().size_of_block(identity) ;
    Matrix result, tmp ;
    extract(bm, sigma, m1, m2, tmp) ;
    if (mod == Left) {
        gemm(tmp, B, result);
    } else if (mod == Right) {
        gemm(B, tmp, result);
    }
    return result ;
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
    // Initialization
    block_matrix<Matrix, SymmGroup> bm = MPS[l].data() ;
    Matrix *tmp ;
    tmp = new Matrix(1,m2) ;
    extract(bm, sigma, m1, m2, *tmp) ;
    std::cout << output << std::endl ;
    if (mod == Left) {
        assert (m1 == 1) ;
        for (int i = 0; i < m2; ++i)
            for (int j = 0; j < m2; ++j)
                output(i,j) = (*tmp)(0,i)*(*tmp)(0,j) ;
    } else if (mod == Right) {
        assert (m2 == 1) ;
        for (int i = 0; i < m1; ++i)
            for (int j = 0; j < m1; ++j)
                output(i,j) = (*tmp)(i,0)*(*tmp)(j,0) ;
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
            std::copy(bm[0].row(offset+i).first, bm[0].row(offset+i).second, output.row(i).first);
};


#endif