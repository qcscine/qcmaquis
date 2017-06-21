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
    partial_overlap(const MPSWave& MPS , const basis_vector& basis);
    // Methods
    basis_type get_basis(const dim_type& site) ;
    void print(void) ;
    void update(const MPSWave& MPS, const dim_type& l, const int& direction) ;
    value_type overlap(const dim_type& i);
    //value_type overlap(const MPSTensor& MPSTns , const dim_type& i) ;
private:
    // Private attributes
    vector_overlap data_left_ , data_right_ ;
    dim_type lattice_L_ ;
    basis_vector basis_ ;
    // Local variables
    enum Modality { Left , Right };
    static const charge identity = SymmGroup::IdentityCharge ;
    // Private function
    void multiply (const MPSWave& MPS, const basis_type& sigma, const dim_type& l,
                   const Modality& mod, vector_overlap& lst, bool modality) ;
    void multiply_first (const MPSWave& MPS, const basis_type& sigma, const Modality& mod, vector_overlap& lst, bool modality) ;
    void extract(const bmatrix& bm, const basis_type& sigma, const std::size_t m1, const std::size_t m2, Matrix& output) ;
};

// +------------+
//  CONSTRUCTORS
// +------------+

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
    MPSTensor MPSTns ;
    // Please note the order of the
    for (dim_type i = 0; i < L; ++i) {
        if (i == 0) {
            multiply_first(MPS, basis_[0],   Left,  data_left_, true);
            multiply_first(MPS, basis_[L-1], Right, data_right_, true);
        } else {
            multiply(MPS, basis_[i], i, Left, data_left_, true) ;
            multiply(MPS, basis_[L-1-i], i, Right, data_right_, true) ;
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
    if (i1l == 0)
        multiply_first(MPS, basis_[0], Left, data_left_, false);
    else if (i1l > 0 && i1l < lattice_L_)
        multiply(MPS, basis_[i1l], i1l, Left, data_left_, false);
    //
    if (i2l == 0)
        multiply_first(MPS, basis_[i2l], Left, data_left_, false);
    else if (i2l > 0 && i2l < lattice_L_)
        multiply(MPS, basis_[i2l], i2l, Left, data_left_, false);
    //
    if (i1r == 0) {
        multiply_first(MPS, basis_[i1r], Right, data_right_, false);
    } else if (i1r > 0 && i1r < lattice_L_) {
        multiply(MPS, basis_[i1r], i1r, Right, data_right_, false);
    }
    //
    if (i2r == 0) {
        multiply_first(MPS, basis_[i2r], Right, data_right_, false);
    } else if (i2r > 0 && i2r < lattice_L_) {
        multiply(MPS, basis_[i2r], i2r, Right, data_right_, false);
    }
    std::cout << "Finalized"  << std::endl ;
    //
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
    } else if (i == lattice_L_-1 ) {
        Matrix a = data_left_[i] ;
        assert (a.num_cols() == 1 && a.num_rows() == 1) ;
        result = a(0,0);
    } else {
        Matrix a = data_left_[i] ;
        Matrix b = data_right_[lattice_L_-1-i] ;
        assert (a.num_rows() == b.num_rows() && a.num_cols() == b.num_cols()) ;
        for (int k1 = 0; k1 < a.num_rows() ; ++k1)
            for (int k2 = 0; k2 < a.num_cols() ; ++k2)
                result += a(k1,k2)*b(k1,k2) ;
    }
    return result ;
};

//template<class Matrix, class SymmGroup>
//typename partial_overlap<Matrix, SymmGroup>::value_type partial_overlap<Matrix, SymmGroup>::overlap(const partial_overlap<Matrix, SymmGroup>::MPSTensor& MPSTns,
//                                                                                                    const dim_type &i)
//{
//    // Check data consistency
//    assert (i < lattice_L_) ;
//    value_type result = 0 ;
//    // Extract the block of the TrivialGroup symmetry from the block matrix associated to the
//    // MPS tensor
//    MPSTns.make_left_paired() ;
//    block_matrix<Matrix, SymmGroup> bm = MPSTns.data() ;
//    std::size_t m1 = MPSTns.row_dim().size_of_block(identity) ;
//    std::size_t m2 = MPSTns.col_dim().size_of_block(identity) ;
//    Matrix *tmp ;
//    tmp = new Matrix(m1, m2) ;
//    extract(bm, sigma, m1, m2, *tmp) ;
//    // Starting site of the
//    if (i == 0) {
//        Matrix a = data_right_[lattice_L_-1] ;
//        assert (tmp.num_rows() == a.num_rows() && tmp.num_cols() == a.num_cols() ) ;
//        for (int k1 = 0; k1 < a.num_rows(); ++k1)
//            for (int k2 = 0; k2 < a.num_cols(); ++k2)
//                result += a(k1, k2)*tmp(k1, k2);
//    } else if (i == lattice_L_-1 ) {
//        Matrix a = data_left_[lattice_L_-2] ;
//        assert (m1 == a.num_cols() && m2 = 1) ;
//        for (int k1 = 0; k1 < a.num_cols() ; ++k1)
//            result += a(1,k1)*tmp(k1,1) ;
//    } else {
//        Matrix a = data_left_[i-1] ;
//        Matrix b = data_right_[lattice_L_-1-i] ;
//        tmp2 = new Matrix(a.num_rows(),m2) ;
//        gemm(transpose(*tmp), a, *tmp2) ;
//        for (int k1 = 0; k1 < (*tmp2).num_rows() ; ++k1)
//            for (int k2 = 0; k2 < (*tmp2).num_cols() ; ++k2)
//                result += (*tmp2)(k1,k2)*b(k1,k2) ;
//    }
//    return result ;
//};

// +---------------+
//  PRIVATE METHODS
// +---------------+
//
// Methods that takes in input:
// 1) an MPS
// 2) a site where the optimization is centered
// 3) the direction of the optimization
// and returns in output the updated list_ with the new overlap
//

template<class Matrix, class SymmGroup>
void partial_overlap<Matrix, SymmGroup>::multiply(const partial_overlap<Matrix,SymmGroup>::MPSWave& MPS,
                                                  const basis_type& sigma,
                                                  const dim_type& l ,
                                                  const partial_overlap<Matrix, SymmGroup>::Modality& mod,
                                                  vector_overlap& lst_,
                                                  bool modality)
{
    // Initialization
    MPSTensor MPSTns ;
    Matrix *tmp , *result ;
    // If we are in the right modality, takes the MPS site starting from the end
    if (mod == Left)
        MPSTns = MPS[l] ;
    else
        MPSTns = MPS[lattice_L_-1 - l] ;
    MPSTns.make_left_paired() ;
    block_matrix<Matrix, SymmGroup> bm = MPSTns.data() ;
    std::size_t m1 = MPSTns.row_dim().size_of_block(identity) ;
    std::size_t m2 = MPSTns.col_dim().size_of_block(identity) ;
    std::size_t m1in = lst_[l-1].num_rows();
    std::size_t m2in = lst_[l-1].num_cols();
    std::cout << m1 << m2 << m1in << m2in << std::endl ;
    tmp = new Matrix(m1, m2) ;
    extract(bm, sigma, m1, m2, *tmp) ;
    if (mod == Left) {
        result = new Matrix(m2, m2in);
        gemm(transpose(*tmp), lst_[l-1], *result);
    }  else {
        result = new Matrix(m1in, m1);
        gemm(lst_[l-1], transpose(*tmp), *result);
    }
    // Finalization
    std::cout << "Finalize" << std::endl ;
    if (modality)
        lst_.push_back(result);
    else
        lst_[l] = *result ;
}

template<class Matrix, class SymmGroup>
void partial_overlap<Matrix, SymmGroup>::multiply_first(const partial_overlap<Matrix,SymmGroup>::MPSWave& MPS,
                                                        const basis_type& sigma,
                                                        const partial_overlap<Matrix, SymmGroup>::Modality& mod,
                                                        vector_overlap& lst_,
                                                        bool modality)
{
    MPSTensor MPSTns ;
    if (mod == Left)
        MPSTns = MPS[0] ;
    else
        MPSTns = MPS[lattice_L_-1] ;
    MPSTns.make_left_paired() ;
    block_matrix<Matrix, SymmGroup> bm = MPSTns.data() ;
    std::size_t m1 = MPSTns.row_dim().size_of_block(identity) ;
    std::size_t m2 = MPSTns.col_dim().size_of_block(identity) ;
    Matrix *tmp, *output ;
    if (mod == Left) {
        tmp = new Matrix(1,m2) ;
        output = new Matrix(m2,m2) ;
        extract(bm, sigma, m1, m2, *tmp) ;
        assert (m1 == 1) ;
        for (int i = 0; i < m2; ++i)
            (*output)(i,i) = (*tmp)(0,i) ;
    } else if (mod == Right) {
        tmp = new Matrix(m1,1) ;
        output = new Matrix(m1,m1) ;
        extract(bm, sigma, m1, m2, *tmp) ;
        assert (m2 == 1) ;
        for (int i = 0; i < m1; ++i)
            (*output)(i,i) = (*tmp)(i,0) ;
    }
    if (modality)
        lst_.push_back(output);
    else
        lst_[0] = *output ;
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

template<class Matrix, class SymmGroup>
void partial_overlap<Matrix, SymmGroup>::print(void)
{
    std::size_t i, j, k ;
    std::cout << " ------------------------------------ " << std::endl ;
    std::cout << " STATUS OF THE PARTIAL OVERLAP OBJECT " << std::endl ;
    std::cout << " ------------------------------------ " << std::endl ;
    for (i = 0; i < lattice_L_ ; ++i){
        Matrix m1, m2 ;
        m1 = data_left_[i] ;
        j  = m1.num_rows() ;
        k  = m1.num_cols() ;
        std::cout << " Block number left - " << i << " is " << j << "x" << k << std::endl ;
        m1 = data_left_[lattice_L_-1-i] ;
        j  = m1.num_rows() ;
        k  = m1.num_cols() ;
        std::cout << " Block number right - " << i << " is " << j << "x" << k << std::endl ;
    }
    std::cout << " ------------------------------------ " << std::endl ;
};

#endif