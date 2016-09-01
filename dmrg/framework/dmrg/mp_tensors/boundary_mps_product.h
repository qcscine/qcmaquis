/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2016 Laboratory of Physical Chemistry, ETH Zurich
 *               2016-2016 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

#ifndef BOUNDARY_MPS_PRODUCT_H
#define BOUNDARY_MPS_PRODUCT_H

#include "dmrg/mp_tensors/boundary.h"
#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/mpotensor.h"

template<class Matrix, class OtherMatrix, class SymmGroup, class Gemm>
class BoundaryMPSProduct 
{
public:
    typedef typename maquis::traits::scalar_type<Matrix>::type scalar_type;
    typedef typename Matrix::value_type value_type;
    typedef typename MPOTensor<Matrix, SymmGroup>::index_type index_type;

    BoundaryMPSProduct(Boundary<OtherMatrix, SymmGroup> const & left_,
                       MPOTensor<Matrix, SymmGroup> const & mpo_) : left(left_), mpo(mpo_), data_(left_.aux_dim()) { }

    std::size_t aux_dim() const { 
        return data_.size(); 
    }

    void resize(size_t n){
        if(n < data_.size()) 
            return data_.resize(n);
        data_.reserve(n);
        for(int i = data_.size(); i < n; ++i)
            data_.push_back(block_matrix<Matrix, SymmGroup>());
    }
    
    void multiply (index_type b1, MPSTensor<Matrix, SymmGroup> const & mps);
    
    block_matrix<Matrix, SymmGroup> & operator[](std::size_t k) { return data_[k]; }
    block_matrix<Matrix, SymmGroup> const & operator[](std::size_t k) const { return data_[k]; }

private:
    std::vector<block_matrix<Matrix, SymmGroup> > data_;
    Boundary<OtherMatrix, SymmGroup> const & left;
    MPOTensor<Matrix, SymmGroup> const & mpo;
};

template <class Matrix, class OtherMatrix, class SymmGroup, class Gemm>
void BoundaryMPSProduct<Matrix, OtherMatrix, SymmGroup, Gemm>::multiply(index_type b1, MPSTensor<Matrix, SymmGroup> const & mps)
{
    if (mpo.herm_info.left_skip(b1))
    {
        parallel::guard group(scheduler(b1), parallel::groups_granularity);
        typename Gemm::gemm_trim_left()(left[mpo.herm_info.left_conj(b1)], mps.data(), data_[b1]);
    }
    else {
        parallel::guard group(scheduler(b1), parallel::groups_granularity);
        typename Gemm::gemm_trim_left()(transpose(left[b1]), mps.data(), data_[b1]);
    }

}


#endif
