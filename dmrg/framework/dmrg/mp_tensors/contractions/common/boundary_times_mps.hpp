/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *                    Laboratory for Physical Chemistry, ETH Zurich
 *               2014-2014 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

#ifndef ENGINE_COMMON_MPS_TIMES_BOUNDDARY_H
#define ENGINE_COMMON_MPS_TIMES_BOUNDDARY_H

#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/mpotensor.h"
#include "dmrg/mp_tensors/reshapes.h"
#include "dmrg/block_matrix/indexing.h"


namespace contraction {
    namespace common {

    template <class Matrix, class OtherMatrix, class SymmGroup>
    typename boost::enable_if<symm_traits::HasSU2<SymmGroup>, std::vector<typename Matrix::value_type> >::type
    conjugate_phases(block_matrix<Matrix, SymmGroup> const & bm,
                     MPOTensor<OtherMatrix, SymmGroup> const & mpo,
                     size_t k, bool left, bool forward)
    {
        typedef typename Matrix::value_type value_type;
        typename SymmGroup::subcharge S = (left) ? mpo.left_spin(k).get() : mpo.right_spin(k).get();

        std::vector<value_type> ret(bm.n_blocks());

        for (size_t b = 0; b < bm.n_blocks(); ++b)
        {
            value_type scale = ::SU2::conjugate_correction<typename Matrix::value_type, SymmGroup>
                                 (bm.basis().left_charge(b), bm.basis().right_charge(b), S);
            if (forward)
                scale *= (left) ? mpo.herm_info.left_phase(mpo.herm_info.left_conj(k))
                                    : mpo.herm_info.right_phase(mpo.herm_info.right_conj(k));
            else
                scale *= (left) ? mpo.herm_info.left_phase(k)
                                    : mpo.herm_info.right_phase(k);
            ret[b] = scale;
        }
        return ret;
    }

    template <class Matrix, class OtherMatrix, class SymmGroup>
    typename boost::disable_if<symm_traits::HasSU2<SymmGroup>, std::vector<typename Matrix::value_type> >::type
    conjugate_phases(block_matrix<Matrix, SymmGroup> const & bm,
                     MPOTensor<OtherMatrix, SymmGroup> const & mpo,
                     size_t k, bool left, bool forward)
    {
        return std::vector<typename Matrix::value_type>(bm.n_blocks(), 1.);
    }

    template <class Matrix, class SymmGroup>
    typename boost::enable_if<symm_traits::HasSU2<SymmGroup> >::type recover_conjugate(block_matrix<Matrix, SymmGroup> & bm,
                                                                                       MPOTensor<Matrix, SymmGroup> const & mpo,
                                                                                       size_t k, bool left, bool forward)
    {
        typedef typename Matrix::value_type value_type;
        std::vector<value_type> scales = conjugate_phases(bm, mpo, k, left, forward);

        for (size_t b = 0; b < bm.n_blocks(); ++b)
            bm[b] *= scales[b];
    }

    template <class Matrix, class SymmGroup>
    typename boost::disable_if<symm_traits::HasSU2<SymmGroup> >::type recover_conjugate(block_matrix<Matrix, SymmGroup> & bm,
                                                                                        MPOTensor<Matrix, SymmGroup> const & mpo,
                                                                                        size_t b, bool left, bool forward)
    { }

    template<class Matrix, class OtherMatrix, class SymmGroup, class Gemm>
    class BoundaryMPSProduct
    {
    public:
        typedef typename maquis::traits::scalar_type<Matrix>::type scalar_type;
        typedef typename Matrix::value_type value_type;
        typedef typename MPOTensor<Matrix, SymmGroup>::index_type index_type;

        BoundaryMPSProduct(MPSTensor<Matrix, SymmGroup> const & mps_,
                           Boundary<OtherMatrix, SymmGroup> const & left_,
                           MPOTensor<Matrix, SymmGroup> const & mpo_) : mps(mps_), left(left_), mpo(mpo_), data_(left_.aux_dim())
        {
            parallel::scheduler_permute scheduler(mpo.placement_l, parallel::groups_granularity);

            int loop_max = left.aux_dim();
            mps.make_right_paired();
            omp_for(int b1, parallel::range(0,loop_max), {

                // exploit single use sparsity (delay multiplication until the object is used)
                if (mpo.num_row_non_zeros(b1) == 1) continue;

                // exploit hermiticity if available
                if (mpo.herm_info.left_skip(b1))
                {
                    parallel::guard group(scheduler(b1), parallel::groups_granularity);

                    std::vector<value_type> scales = conjugate_phases(left[mpo.herm_info.left_conj(b1)], mpo, b1, true, false);
                    typename Gemm::gemm_trim_left()(left[mpo.herm_info.left_conj(b1)], mps.data(), data_[b1], scales);
                }
                else {
                    parallel::guard group(scheduler(b1), parallel::groups_granularity);
                    typename Gemm::gemm_trim_left()(transpose(left[b1]), mps.data(), data_[b1]);
                }
            });
        }

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

        void multiply (index_type b1);

        block_matrix<Matrix, SymmGroup> & operator[](std::size_t k) { return data_[k]; }
        block_matrix<Matrix, SymmGroup> const & operator[](std::size_t k) const { return data_[k]; }

        block_matrix<Matrix, SymmGroup> const & at(std::size_t k, block_matrix<Matrix, SymmGroup> & storage) const
        {
            if (mpo.num_row_non_zeros(k) == 1)
            {
                if (mpo.herm_info.left_skip(k))
                {
                    //parallel::guard group(scheduler(b1), parallel::groups_granularity);
                    //typename Gemm::gemm_trim_left()(left[mpo.herm_info.left_conj(k)], mps.data(), storage);
                    std::vector<value_type> scales = conjugate_phases(left[mpo.herm_info.left_conj(k)], mpo, k, true, false);
                    typename Gemm::gemm_trim_left()(left[mpo.herm_info.left_conj(k)], mps.data(), storage, scales);
                }
                else {
                    //parallel::guard group(scheduler(b1), parallel::groups_granularity);
                    typename Gemm::gemm_trim_left()(transpose(left[k]), mps.data(), storage);
                }

                return storage;
            }

            else
                return data_[k];
        }

    private:
        std::vector<block_matrix<Matrix, SymmGroup> > data_;

        MPSTensor<Matrix, SymmGroup> const & mps;
        Boundary<OtherMatrix, SymmGroup> const & left;
        MPOTensor<Matrix, SymmGroup> const & mpo;
    };


    template <class Matrix, class OtherMatrix, class SymmGroup, class Gemm>
    void BoundaryMPSProduct<Matrix, OtherMatrix, SymmGroup, Gemm>::multiply(index_type b1)
    {
    }

    template<class Matrix, class OtherMatrix, class SymmGroup, class Gemm>
    class MPSBoundaryProduct
    {
    public:
        typedef typename maquis::traits::scalar_type<Matrix>::type scalar_type;
        typedef typename Matrix::value_type value_type;
        typedef typename MPOTensor<Matrix, SymmGroup>::index_type index_type;

        MPSBoundaryProduct(MPSTensor<Matrix, SymmGroup> const & mps_,
                           Boundary<OtherMatrix, SymmGroup> const & right_,
                           MPOTensor<Matrix, SymmGroup> const & mpo_) : mps(mps_), right(right_), mpo(mpo_), data_(right_.aux_dim())
                                                                      , pop_(right_.aux_dim(), 0)
        {
            parallel::scheduler_permute scheduler(mpo.placement_r, parallel::groups_granularity);

            int loop_max = right.aux_dim();
            mps.make_left_paired();
            omp_for(int b2, parallel::range(0,loop_max), {

                // exploit single use sparsity (delay multiplication until the object is used)
                if (mpo.num_col_non_zeros(b2) == 1) continue;

                // exploit hermiticity if available
                if (mpo.herm_info.right_skip(b2))
                {
                    parallel::guard group(scheduler(b2), parallel::groups_granularity);
                    block_matrix<typename maquis::traits::transpose_view<Matrix>::type, SymmGroup> trv = transpose(right[mpo.herm_info.right_conj(b2)]);
                    std::vector<value_type> scales = conjugate_phases(trv, mpo, b2, false, true);

                    typename Gemm::gemm_trim_right()(mps.data(), trv, data_[b2], scales);
                }
                else {
                    parallel::guard group(scheduler(b2), parallel::groups_granularity);
                    typename Gemm::gemm_trim_right()(mps.data(), right[b2], data_[b2]);
                }
            });
        }

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

        void multiply (index_type b2);

        block_matrix<Matrix, SymmGroup> & operator[](index_type k) { return data_[k]; }
        block_matrix<Matrix, SymmGroup> const & operator[](index_type k) const { return data_[k]; }

        DualIndex<SymmGroup> basis_at(index_type k) const
        {
            if (mpo.num_col_non_zeros(k) == 1)
            {
                if (mpo.herm_info.right_skip(k))
                {
                    //parallel::guard group(scheduler(b2), parallel::groups_granularity);
                    block_matrix<typename maquis::traits::transpose_view<Matrix>::type, SymmGroup> trv = transpose(right[mpo.herm_info.right_conj(k)]);
                    //return typename Gemm::gemm_trim_right()(mps.data(), trv);
                    return SU2::gemm_trim_right_pretend(mps.data(), trv);
                }
                else {
                    //parallel::guard group(scheduler(b2), parallel::groups_granularity);
                    //return typename Gemm::gemm_trim_right_pretend()(mps.data(), right[k]);
                    return SU2::gemm_trim_right_pretend(mps.data(), right[k]);
                }
            }
            else
                return data_[k].basis();
        }

        //block_matrix<Matrix, SymmGroup> const & at(std::size_t k, block_matrix<Matrix, SymmGroup> & storage) const
        block_matrix<Matrix, SymmGroup> const & at(index_type k) const
        {
            if (mpo.num_col_non_zeros(k) == 1)
            {
                if (!pop_[k])
                {
                    if (mpo.herm_info.right_skip(k))
                    {
                        //parallel::guard group(scheduler(b2), parallel::groups_granularity);
                        //typename Gemm::gemm_trim_right()(mps.data(), transpose(right[mpo.herm_info.right_conj(k)]), storage);
                        block_matrix<typename maquis::traits::transpose_view<Matrix>::type, SymmGroup> trv = transpose(right[mpo.herm_info.right_conj(k)]);
                        std::vector<value_type> scales = conjugate_phases(trv, mpo, k, false, true);
                        //typename Gemm::gemm_trim_right()(mps.data(), trv, storage, scales);
                        typename Gemm::gemm_trim_right()(mps.data(), trv, data_[k], scales);
                    }
                    else {
                        //parallel::guard group(scheduler(b2), parallel::groups_granularity);
                        //typename Gemm::gemm_trim_right()(mps.data(), right[k], storage);
                        typename Gemm::gemm_trim_right()(mps.data(), right[k], data_[k]);
                    }
                    pop_[k] = 1;
                }

                return data_[k];
            }
            else
                return data_[k];
        }

        void free(index_type b1) const
        {
            for (index_type b2 = 0; b2 < mpo.col_dim(); ++b2)
                if (mpo.num_col_non_zeros(b2) == 1)
                    if (mpo.has(b1,b2))
                    {
                        data_[b2].clear();
                        break;
                    }
        }

    private:
        mutable std::vector<block_matrix<Matrix, SymmGroup> > data_;
        mutable std::vector<char> pop_;

        MPSTensor<Matrix, SymmGroup> const & mps;
        Boundary<OtherMatrix, SymmGroup> const & right;
        MPOTensor<Matrix, SymmGroup> const & mpo;
    };

    template <class Matrix, class OtherMatrix, class SymmGroup, class Gemm>
    void MPSBoundaryProduct<Matrix, OtherMatrix, SymmGroup, Gemm>::multiply(index_type b2)
    {
    }

    } // namespace common
} // namespace contraction

#endif
