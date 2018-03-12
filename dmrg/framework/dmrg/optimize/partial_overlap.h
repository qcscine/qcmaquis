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
    typedef          int                                            basis_type ;
    typedef typename std::vector<basis_type>                        basis_vector ;
    typedef typename block_matrix<Matrix,SymmGroup>::block_matrix   bmatrix ;
    typedef typename SymmGroup::charge                              charge ;
    typedef typename std::size_t                                    dim_type ;
    typedef typename MPSTensor<Matrix,SymmGroup>::MPSTensor         MPSTensor ;
    typedef typename MPS<Matrix,SymmGroup>::MPS                     MPSWave ;
    typedef typename Matrix::value_type                             value_type ;
    typedef typename boost::ptr_vector<Matrix>                      vector_overlap ;
    // Constructors
    partial_overlap(void) ;
    partial_overlap(const MPSWave& MPS , const basis_vector& basis);
    partial_overlap(const MPSWave& MPS , const MPSWave& MPSOther);
    // Methods
    bool is_defined(void) ;
    bool is_general(void) ;
    basis_type get_basis(const dim_type& site) ;
    void print(void) ;
    void update(const MPSWave& MPS, const dim_type& l, const int& direction) ;
    value_type overlap(const dim_type& i);
    value_type overlap(const MPSTensor& MPSTns , const dim_type& i) ;
    value_type overlap(const MPSTensor& MPSTns , const MPSTensor& MPSOther, const dim_type& i) ;
    value_type overlap(const MPSTensor& MPSTns , const dim_type& i, const dim_type& j) ;
    value_type overlap(const MPSTensor& MPSTns , const MPSTensor& MPSOther, const dim_type& i, const dim_type& j) ;
private:
    // Private attributes
    vector_overlap data_left_ , data_right_ ;
    dim_type lattice_L_ ;
    basis_vector basis_ , phys_sizes_ ;
    bool is_defined_ , is_general_mps_ ;
    // Local variables
    enum Modality { Left , Right };
    static const charge identity = SymmGroup::IdentityCharge ;
    // Private function
    Matrix *multiply (const MPSTensor& MPSTns,  const Matrix& input, const basis_type& sigma, const Modality& mod,
                      bool modality) ;
    Matrix *multiply (const MPSTensor& MPSTns1, const MPSTensor& MPSTns2, const Matrix& input, const Modality& mod, bool modality) ;
    Matrix *multiply (const MPSTensor& MPSTns, const Matrix& input, const basis_type& sigma1, const basis_type& sigma2,
                      const Modality& mod, bool modality) ;
    Matrix* multiply_first (const MPSTensor& MPSTns, const basis_type& sigma, const Modality& mod,
                            bool modality) ;
    Matrix* multiply_first (const MPSTensor& MPSTns, const MPSTensor& MPSOther, const Modality& mod, bool modality) ;
    Matrix* multiply_first (const MPSTensor& MPSTns, const basis_type& sigma1, const basis_type& sigma2,
                            const Modality& mod, bool modality) ;
    void extract(const bmatrix& bm, const basis_type& sigma, const std::size_t m1, const std::size_t m2, Matrix& output) ;
    void extract(const bmatrix& bm, const basis_type& sigma1, const basis_type& sigma2, const std::size_t m1,
                 const std::size_t m2, Matrix& output) ;
    void multiply (const MPSWave& MPS, const dim_type& l, const Modality& mod, bool modality, vector_overlap& lst) ;
    void multiply (const MPSWave& MPS, const MPSWave& MPSOther, const dim_type& l, const Modality& mod,
                   bool modality, vector_overlap& lst) ;
    void multiply_first (const MPSWave& MPS, const Modality& mod, bool modality, vector_overlap& lst) ;
    void multiply_first (const MPSWave& MPS, const MPSWave& MPSOther, const Modality& mod, bool modality, vector_overlap& lst) ;
};

// +------------+
//  CONSTRUCTORS
// +------------+

template<class Matrix, class SymmGroup>
partial_overlap<Matrix,SymmGroup>::partial_overlap(void) :
    lattice_L_(0) , is_defined_(false) , is_general_mps_(false) { } ;

template<class Matrix, class SymmGroup>
partial_overlap<Matrix,SymmGroup>::partial_overlap(const MPSWave& MPS, const basis_vector& basis)
{
    dim_type L = MPS.size() ;
    if (basis.size() == 0) {
        is_defined_ = false ;
        is_general_mps_ = false ;
        return ;
    } else {
        is_defined_ = true ;
        for (int i = 0; i < L; ++i)
            basis_.push_back(basis[i]);
        lattice_L_ = L;
        MPSTensor MPSTns;
        //
        for (dim_type i = 0; i < L; ++i) {
            phys_sizes_.push_back(MPS.site_dim(i).size_of_block(identity)) ;
            if (i == 0) {
                multiply_first(MPS, Left, true, data_left_);
                multiply_first(MPS, Right, true, data_right_);
            } else {
                multiply(MPS, i, Left, true, data_left_);
                multiply(MPS, i, Right, true, data_right_);
            }
        }
    }
};

template<class Matrix, class SymmGroup>
partial_overlap<Matrix,SymmGroup>::partial_overlap(const MPSWave& MPS, const MPSWave& MPSOther)
{
    is_defined_ = true ;
    dim_type L1  = MPS.size() , L2 = MPSOther.size() ;
    assert (L1 == L2) ;
    lattice_L_ = L1 ;
    //
    multiply_first(MPS, MPSOther, Left, true, data_left_);
    multiply_first(MPS, MPSOther, Right, true, data_right_);
    for (dim_type i = 1; i < lattice_L_ ; ++i) {
        multiply(MPS, MPSOther, i, Left, true, data_left_);
        multiply(MPS, MPSOther, i, Right, true, data_right_);
        phys_sizes_.push_back(MPS.site_dim(i).size_of_block(identity)) ;
    }
};

// +-------+
//  METHODS
// +-------+

template<class Matrix, class SymmGroup>
bool partial_overlap<Matrix, SymmGroup>::is_defined(void)
{
    return is_defined_ ;
}

template<class Matrix, class SymmGroup>
bool partial_overlap<Matrix, SymmGroup>::is_general(void)
{
    return is_general_mps_ ;
}

template<class Matrix, class SymmGroup>
typename partial_overlap<Matrix, SymmGroup>::basis_type partial_overlap<Matrix, SymmGroup>::get_basis(const dim_type &site)
{
    return basis_[site] ;
}

template<class Matrix, class SymmGroup>
void partial_overlap<Matrix, SymmGroup>::update(const MPSWave& MPS,
                                                const dim_type& l,
                                                const int& direction)
{
    // Get the index to update in the left and right overlaps
    dim_type i1l, i2l;
    dim_type i1r, i2r;
    Modality mod ;
    if (direction == 1)
        mod = Left;
    else if (direction == -1)
        mod = Right;
    // Set the sites where the
    if (mod == Left) {
        i1l = l;
        i2l = l + 1;
        i1r = lattice_L_ - 1 - i2l ;
        i2r = lattice_L_ - 1 - i1l ;
    } else {
        i1l = l - 1 ;
        i2l = l ;
        i1r = lattice_L_ - 1 - i2l ;
        i2r = lattice_L_ - 1 - i1l ;
    }
    if (i1l == 0) {
        multiply_first(MPS, Left, false, data_left_);
    } else if (i1l > 0 && i1l < lattice_L_)
        multiply(MPS, i1l, Left, false, data_left_);
    //
    if (i2l == 0) {
        multiply_first(MPS, Left, false, data_left_);
    } else if (i2l > 0 && i2l < lattice_L_) {
        multiply(MPS, i2l, Left, false, data_left_);
    }
    //
    if (i1r == 0) {
        multiply_first(MPS, Right, false, data_right_);
    } else if (i1r > 0 && i1r < lattice_L_) {
        multiply(MPS, i1r, Right, false, data_right_);
    }
    //
    if (i2r == 0) {
        multiply_first(MPS, Right, false, data_right_);
    } else if (i2r > 0 && i2r < lattice_L_) {
        multiply(MPS, i2r, Right, false, data_right_);
    }
    //
};

// +-------------------------------+
//  METHODS TO COMPUTE THE OVERLAPS
// +-------------------------------+

template<class Matrix, class SymmGroup>
typename partial_overlap<Matrix, SymmGroup>::value_type partial_overlap<Matrix, SymmGroup>::overlap(const dim_type &i)
{
    // Check data consistency
    assert (i < lattice_L_) ;
    value_type result = 0 ;
    Matrix a = data_left_[i] ;
    if (i == lattice_L_-1) {
        assert (a.num_cols() == a.num_rows() == 1) ;
        result = a(0,0);
    } else {
        Matrix b = data_right_[lattice_L_-2-i] ;
        assert (a.num_rows() == b.num_rows() && a.num_cols() == b.num_cols()) ;
        for (int k1 = 0; k1 < a.num_rows() ; ++k1)
            for (int k2 = 0; k2 < a.num_cols() ; ++k2)
                result += a(k1,k2)*b(k1,k2) ;
    }
    return result ;
};

template<class Matrix, class SymmGroup>
typename partial_overlap<Matrix, SymmGroup>::value_type partial_overlap<Matrix, SymmGroup>::overlap(const MPSTensor& MPSTns,
                                                                                                    const dim_type &i)
{
    // Check data consistency and declaration
    MPSTns.make_left_paired() ;
    assert (i < lattice_L_) ;
    value_type result = 0 ;
    dim_type indx = lattice_L_-2-i ;
    Matrix *tmp , *tmp2 ;
    if ( i == lattice_L_-1 ) {
        tmp2 = this->multiply(MPSTns, data_left_[i-1], basis_[i], Left, true) ;
        assert ((*tmp2).num_rows() == (*tmp2).num_cols() == 1) ;
        result = (*tmp2)(0,0) ;
    } else {
        if (i == 0)
            tmp2 = this->multiply_first(MPSTns, basis_[0], Left, true);
        else
            tmp2 = this->multiply(MPSTns, data_left_[i-1], basis_[i], Left, true);
        assert((*tmp2).num_rows() == data_right_[indx].num_cols() &&
               (*tmp2).num_cols() == data_right_[indx].num_rows());
        for (int k1 = 0; k1 < (*tmp2).num_rows(); ++k1)
            for (int k2 = 0; k2 < (*tmp2).num_cols(); ++k2) {
                result += (*tmp2)(k1, k2)*data_right_[indx](k2, k1);
        }
    }
    return result ;
};

template<class Matrix, class SymmGroup>
typename partial_overlap<Matrix, SymmGroup>::value_type partial_overlap<Matrix, SymmGroup>::overlap(const MPSTensor& MPSTns,
                                                                                                    const MPSTensor& MPSOther,
                                                                                                    const dim_type &i)
{
    // Check data consistency and declaration
    assert (i < lattice_L_) ;
    dim_type indx = lattice_L_-2-i ;
    dim_type NM1 = MPSTns.site_dim.size_of_block(identity) ;
    dim_type NM2 = MPSOther.site_dim.size_of_block(identity) ;
    assert (NM1 == NM2) ;
    value_type result = 0 ;
    Matrix *tmp ;
    if ( i == lattice_L_-1 ) {
        //assert (m2 == 1);
        tmp = multiply(MPSTns, MPSOther, data_left_[i-1], NM1, Left, true) ;
        assert ((*tmp).num_rows() == (*tmp).num_cols() == 1) ;
        result = (*tmp)(0,0) ;
    } else {
        if (i == 0) {
            tmp = multiply_first(MPSTns, MPSOther, NM1, Left, true);
        } else {
            tmp = multiply(MPSTns, data_left_[i-1], basis_[i], Left, true);
        }
        assert((*tmp).num_rows() == data_right_[indx].num_cols() &&
               (*tmp).num_cols() == data_right_[indx].num_rows());
        for (int k1 = 0; k1 < (*tmp).num_rows(); ++k1)
            for (int k2 = 0; k2 < (*tmp).num_cols(); ++k2) {
                result += (*tmp)(k1, k2)*data_right_[indx](k2, k1);
            }
    }
    return result ;
};

template<class Matrix, class SymmGroup>
typename partial_overlap<Matrix, SymmGroup>::value_type partial_overlap<Matrix, SymmGroup>::overlap(const MPSTensor& MPSTns,
                                                                                                    const dim_type& i1,
                                                                                                    const dim_type& i2)
{
    // Check data consistency and declaration
    assert ( i1 < lattice_L_ && i2 < lattice_L_ ) ;
    value_type result = 0 ;
    dim_type indx = lattice_L_-2-i2 ;
    Matrix *tmp , *tmp2 ;
    if ( i2 == lattice_L_-1 ) {
        //assert (m2 == 1);
        tmp2 = this->multiply(MPSTns, data_left_[i1-1], basis_[i1], basis_[i2], Left, true) ;
        assert ((*tmp2).num_rows() == (*tmp2).num_cols() == 1) ;
        result = (*tmp2)(0,0) ;
    } else {
        if (i1 == 0) {
            //assert (m1 == 1);
            tmp2 = this->multiply_first(MPSTns, basis_[i1], basis_[i2], Left, true);
        } else {
            tmp2 = this->multiply(MPSTns, data_left_[i1-1], basis_[i1], basis_[i2], Left, true);
        }
        assert((*tmp2).num_rows() == data_right_[indx].num_cols() &&
               (*tmp2).num_cols() == data_right_[indx].num_rows());
        for (int k1 = 0; k1 < (*tmp2).num_rows(); ++k1)
            for (int k2 = 0; k2 < (*tmp2).num_cols(); ++k2) {
                result += (*tmp2)(k1, k2)*data_right_[indx](k2, k1);
            }
    }
    return result ;
};

// +-----------------------+
//  MULTIPLICATION ROUTINES
// +-----------------------+

template<class Matrix, class SymmGroup>
void partial_overlap<Matrix, SymmGroup>::multiply(const MPSWave& MPS,
                                                  const dim_type& l ,
                                                  const Modality& mod,
                                                  bool modality,
                                                  vector_overlap& data_)
{
    // Initialization
    MPSTensor MPSTns ;
    Matrix *output ;
    basis_type sigma ;
    //
    if (mod == Left) {
        MPSTns = MPS[l] ;
        sigma = basis_[l] ;
    } else {
        MPSTns = MPS[lattice_L_-1 - l] ;
        sigma = basis_[lattice_L_-1 -l ] ;
    }
    MPSTns.make_left_paired();
    output = this->multiply(MPSTns, data_[l-1], sigma, mod, modality) ;
    // Finalization
    if (modality)
        data_.push_back(output);
    else
        data_[l] = *output ;
}


template<class Matrix, class SymmGroup>
void partial_overlap<Matrix, SymmGroup>::multiply(const MPSWave& MPS,
                                                  const MPSWave& MPSOther,
                                                  const dim_type& l ,
                                                  const Modality& mod,
                                                  bool modality,
                                                  vector_overlap& data_)
{
    // Initialization
    MPSTensor MPSTns1, MPSTns2 ;
    Matrix *output ;
    //
    if (mod == Left) {
        MPSTns1 = MPS[l] ;
        MPSTns2 = MPSOther[l] ;
    } else {
        MPSTns1 = MPS[lattice_L_-1 - l] ;
        MPSTns2 = MPSOther[lattice_L_-1 -l ] ;
    }
    output = multiply(MPSTns1, MPSTns2, data_[l-1], mod, modality) ;
    // Finalization
    if (modality)
        data_.push_back(output);
    else
        data_[l] = *output ;
}

template<class Matrix, class SymmGroup>
Matrix* partial_overlap<Matrix, SymmGroup>::multiply(const MPSTensor& MPSTns,
                                                     const Matrix& input ,
                                                     const basis_type& sigma,
                                                     const Modality& mod,
                                                     bool modality)
{
    // Initialization
    Matrix *tmp , *result , *result2;
    bmatrix bm, output ;
    bm = MPSTns.data() ;
    std::size_t m1 = MPSTns.row_dim().size_of_block(identity) ;
    std::size_t m2 = MPSTns.col_dim().size_of_block(identity) ;
    std::size_t m1in = input.num_rows();
    std::size_t m2in = input.num_cols();
    tmp     = new Matrix(m1, m2) ;
    extract(bm, sigma, m1, m2, *tmp) ;
    // Matrix multiplication calculation
    if (mod == Left) {
        result2  = new Matrix(m2, m2in) ;
        gemm(transpose(*tmp), input, *result2);
    }  else {
        result2  = new Matrix(m1in, m1) ;
        gemm(input, transpose(*tmp), *result2);
    }
    return result2 ;
}


template<class Matrix, class SymmGroup>
Matrix* partial_overlap<Matrix, SymmGroup>::multiply(const MPSTensor& MPSTns1,
                                                     const MPSTensor& MPSTns2,
                                                     const Matrix& input ,
                                                     const Modality& mod,
                                                     bool modality)
{
    // Initialization
    Matrix *result_buffer, *result_intermediate , *result_final;
    Matrix *tmp1, *tmp2 ;
    bmatrix bm1, bm2, output ;
    MPSTns1.make_left_paired() ;
    MPSTns2.make_left_paired() ;
    bm1 = MPSTns1.data() ;
    bm2 = MPSTns2.data() ;
    size_t m1_1 = MPSTns1.row_dim().size_of_block(identity) ;
    size_t m2_1 = MPSTns1.col_dim().size_of_block(identity) ;
    size_t m1_2 = MPSTns2.row_dim().size_of_block(identity) ;
    size_t m2_2 = MPSTns2.col_dim().size_of_block(identity) ;
    size_t m1in = input.num_rows();
    size_t m2in = input.num_cols();
    size_t NM1  = MPSTns1.site_dim().size_of_block(identity) ;
    size_t NM2  = MPSTns2.site_dim().size_of_block(identity) ;
    assert (NM1 == NM2) ;
    size_t NM = NM1 ;
    // Matrix multiplication calculation
    if (mod == Left) {
        tmp1 = new Matrix(m1_1, m2_1) ;
        tmp2 = new Matrix(m1_2, m2_2) ;
        result_intermediate  = new Matrix(m2_1, m2in) ;
        result_final         = new Matrix(m2_1, m2_2) ;
        result_buffer        = new Matrix(m2_1, m2_2) ;
        extract(bm1, 0, m1_1, m2_1, *tmp1) ;
        extract(bm2, 0, m1_2, m2_2, *tmp2) ;
        gemm(transpose(*tmp1), input, *result_intermediate);
        gemm(*result_intermediate, *tmp2, *result_final) ;
        for (int k = 1 ; k < NM ; k++){
            extract(bm1, k, m1_1, m2_1, *tmp1) ;
            extract(bm2, k, m1_2, m2_2, *tmp2) ;
            gemm(transpose(*tmp1), input, *result_intermediate);
            gemm(*result_intermediate, *tmp2, *result_buffer) ;
            *result_final += *result_buffer ;
        }
    }  else {
        tmp1 = new Matrix(m1_1, m2_1) ;
        tmp2 = new Matrix(m1_2, m2_2) ;
        result_intermediate  = new Matrix(m1_1, m1in) ;
        result_final         = new Matrix(m1_1, m1_2) ;
        result_buffer        = new Matrix(m1_2, m1_2) ;
        extract(bm1, 0, m1_1, m2_1, *tmp1) ;
        extract(bm2, 0, m1_2, m2_2, *tmp2) ;
        gemm(*tmp1, input, *result_intermediate);
        gemm(*result_intermediate, transpose(*tmp2), *result_final) ;
        for (int k = 1 ; k < NM ; k++){
            extract(bm1, k, m1_1, m2_1, *tmp1) ;
            extract(bm2, k, m1_2, m2_2, *tmp2) ;
            gemm(*tmp1, input, *result_intermediate);
            gemm(*result_intermediate, transpose(*tmp2), *result_buffer) ;
            *result_final += *result_buffer ;
        }
    }
    return result_final ;
}

template<class Matrix, class SymmGroup>
Matrix* partial_overlap<Matrix, SymmGroup>::multiply(const MPSTensor& MPSTns,
                                                     const Matrix& input ,
                                                     const basis_type& sigma1,
                                                     const basis_type& sigma2,
                                                     const Modality& mod,
                                                     bool modality)
{
    // Initialization
    MPSTns.make_left_paired() ;
    Matrix *tmp , *result , *result2;
    bmatrix bm, output ;
    bm = MPSTns.data() ;
    std::size_t m1 = MPSTns.row_dim().size_of_block(identity) ;
    std::size_t m2 = MPSTns.col_dim().size_of_block(identity) ;
    std::size_t m1in = input.num_rows();
    std::size_t m2in = input.num_cols();
    tmp     = new Matrix(m1, m2) ;
    extract(bm, sigma1, sigma2, m1, m2, *tmp) ;
    // Matrix multiplication calculation
    if (mod == Left) {
        result2 = new Matrix(m2, m2in) ;
        gemm(transpose(*tmp), input, *result2);
    }  else {
        result2 = new Matrix(m1in, m1) ;
        gemm(input, transpose(*tmp), *result2);
    }
    return result2 ;
}

template<class Matrix, class SymmGroup>
void partial_overlap<Matrix, SymmGroup>::multiply_first(const MPSWave& MPS,
                                                        const Modality& mod,
                                                        bool modality,
                                                        vector_overlap& data_)
{
    MPSTensor MPSTns ;
    basis_type sigma ;
    Matrix *output ;
    if (mod == Left) {
        MPSTns = MPS[0];
        sigma = basis_[0];
    } else {
        MPSTns = MPS[lattice_L_-1];
        sigma  = basis_[lattice_L_-1];
    }
    MPSTns.make_left_paired() ;
    output = multiply_first(MPSTns, sigma, mod, modality) ;
    if (modality)
        data_.push_back(output);
    else
        data_[0] = *output ;
    return ;
}

template<class Matrix, class SymmGroup>
void partial_overlap<Matrix, SymmGroup>::multiply_first(const MPSWave& MPS,
                                                        const MPSWave& MPSOther,
                                                        const Modality& mod,
                                                        bool modality,
                                                        vector_overlap& data_)
{
    MPSTensor MPSTns1, MPSTns2 ;
    Matrix *output ;
    if (mod == Left) {
        MPSTns1 = MPS[0] ;
        MPSTns2 = MPSOther[0] ;
    } else {
        MPSTns1 = MPS[lattice_L_-1];
        MPSTns2 = MPSOther[lattice_L_-1];
    }
    MPSTns1.make_left_paired() ;
    MPSTns2.make_left_paired() ;
    output = multiply_first(MPSTns1, MPSTns2, mod, modality) ;
    if (modality)
        data_.push_back(output);
    else
        data_[0] = *output ;
    return ;
}

template<class Matrix, class SymmGroup>
Matrix* partial_overlap<Matrix, SymmGroup>::multiply_first(const MPSTensor& MPSTns,
                                                           const basis_type& sigma,
                                                           const Modality& mod,
                                                           bool modality)
{
    bmatrix bm = MPSTns.data() ;
    size_t m1  = MPSTns.row_dim().size_of_block(identity) ;
    size_t m2  = MPSTns.col_dim().size_of_block(identity) ;
    Matrix *tmp, *tmp2;
    if (mod == Left) {
        tmp  = new Matrix(1,m2) ;
        tmp2 = new Matrix(m2,1) ;
        extract(bm, sigma, 1, m2, *tmp) ;
        for (size_t i = 0; i < m2; i++)
            (*tmp2)(i,0) = (*tmp)(0,i) ;
    } else if (mod == Right) {
        tmp  = new Matrix(m1,1) ;
        tmp2 = new Matrix(1,m1) ;
        extract(bm, sigma, m1, 1, *tmp) ;
        for (size_t i = 0; i < m1; i++)
            (*tmp2)(0,i) = (*tmp)(i,0) ;
    }
    return tmp2 ;
}

template<class Matrix, class SymmGroup>
Matrix* partial_overlap<Matrix, SymmGroup>::multiply_first(const MPSTensor& MPSTns,
                                                           const MPSTensor& MPSOther,
                                                           const Modality& mod,
                                                           bool modality)
{
    bmatrix bm1 = MPSTns.data() ;
    bmatrix bm2 = MPSOther.data() ;
    size_t m1_1 = MPSTns.row_dim().size_of_block(identity) ;
    size_t m2_1 = MPSTns.col_dim().size_of_block(identity) ;
    size_t m1_2 = MPSOther.row_dim().size_of_block(identity) ;
    size_t m2_2 = MPSOther.col_dim().size_of_block(identity) ;
    size_t NM1  = MPSTns.site_dim().size_of_block(identity) ;
    size_t NM2  = MPSOther.site_dim().size_of_block(identity) ;
    assert (m1_1 == m2_1 == 1 && NM1 == NM2) ;
    size_t NM = NM1 ;
    Matrix *tmp1, *tmp2, *output , *output_buffer ;
    if (mod == Left) {
        tmp1          = new Matrix(m1_1 , m2_1) ;
        tmp2          = new Matrix(m1_2 , m2_2) ;
        output        = new Matrix(m2_1 , m2_2) ;
        output_buffer = new Matrix(m2_1 , m2_2) ;
        extract(bm1, 0, m1_1, m2_1, *tmp1) ;
        extract(bm2, 0, m1_2, m2_2, *tmp2) ;
        gemm(transpose(*tmp1),*tmp2,*output);
        for (int i = 0; i < NM; ++i) {
            tmp1          = new Matrix(m1_1 , m2_1) ;
            extract(bm1, i, m1_1, m2_1, *tmp1) ;
            extract(bm2, i, m1_2, m2_2, *tmp2) ;
            gemm(transpose(*tmp1) , *tmp2 , *output_buffer) ;
            *output += *output_buffer ;
        }
    } else if (mod == Right) {
        tmp1          = new Matrix(m1_1 , m2_1) ;
        tmp2          = new Matrix(m1_2 , m2_2) ;
        output        = new Matrix(m1_1 , m1_2) ;
        output_buffer = new Matrix(m1_1 , m1_2) ;
        extract(bm1, 0, m1_1, m2_1, *tmp1) ;
        extract(bm2, 0, m1_2, m2_2, *tmp2) ;
        gemm(*tmp1,transpose(*tmp2),*output);
        for (int i = 0; i < NM; ++i) {
            extract(bm1, i, m1_1, m2_1, *tmp1) ;
            extract(bm2, i, m1_2, m2_2, *tmp2) ;
            gemm(*tmp1 , transpose(*tmp2) , *output_buffer) ;
            *output += *output_buffer ;
        }
    }
    return output ;
}

template<class Matrix, class SymmGroup>
Matrix* partial_overlap<Matrix, SymmGroup>::multiply_first(const MPSTensor& MPSTns,
                                                           const basis_type& sigma1,
                                                           const basis_type& sigma2,
                                                           const Modality& mod,
                                                           bool modality)
{
    MPSTns.make_left_paired() ;
    bmatrix bm = MPSTns.data() ;
    std::size_t m1 = MPSTns.row_dim().size_of_block(identity) ;
    std::size_t m2 = MPSTns.col_dim().size_of_block(identity) ;
    Matrix *tmp , *tmp2 ;
    if (mod == Left) {
        tmp  = new Matrix(1,m2) ;
        tmp2 = new Matrix(m2,1) ;
        extract(bm, sigma1, sigma2, 1, m2, *tmp) ;
        for (size_t i = 0; i < m2; i++)
            (*tmp2)(i,0) = (*tmp)(0,i) ;
    } else if (mod == Right) {
        tmp  = new Matrix(m1,1) ;
        tmp2 = new Matrix(1,m1) ;
        extract(bm, sigma1, sigma2, m1, 1, *tmp) ;
        for (size_t i = 0; i < m1; i++)
            (*tmp2)(0,i) = (*tmp)(i,0) ;
    }
    return tmp2 ;
}

// +---------------+
//  PRIVATE METHODS
// +---------------+

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
            output(i, j) = bm[0](i+offset, j);
};

template<class Matrix, class SymmGroup>
void partial_overlap<Matrix, SymmGroup>::extract(const bmatrix& bm,
                                                 const basis_type& sigma1,
                                                 const basis_type& sigma2,
                                                 const std::size_t m1,
                                                 const std::size_t m2,
                                                 Matrix& output)
{
    dim_type NMax = phys_sizes_[sigma2] ;
    std::size_t offset = m1*NMax*sigma1 + m1*sigma2 ;
    for (int i = 0; i < m1; ++i)
        for (int j = 0 ; j < m2 ; ++j) {
            output(i, j) = bm[0](i+offset,j);
        }
};

template<class Matrix, class SymmGroup>
void partial_overlap<Matrix, SymmGroup>::print(void)
{
    std::size_t i, j, k ;
    std::cout << " ------------------------------------ " << std::endl ;
    std::cout << " STATUS OF THE PARTIAL OVERLAP OBJECT " << std::endl ;
    std::cout << " ------------------------------------ " << std::endl ;
    Matrix m1 ;
    for (i = 0; i < lattice_L_ ; ++i) {
        m1 = data_left_[i];
        j = m1.num_rows();
        k = m1.num_cols();
        std::cout << " Block number left - " << i << " is " << j << "x" << k << std::endl;
    }
    for (i = 0; i < lattice_L_ ; ++i){
        m1 = data_right_[i] ;
        j  = m1.num_rows() ;
        k  = m1.num_cols() ;
        std::cout << " Block number right - " << i << " is " << j << "x" << k << std::endl ;
    }
    std::cout << " ------------------------------------ " << std::endl ;
};

#endif
