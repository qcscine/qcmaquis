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

#ifndef CONTRACTIONS_SU2_SITE_HAMIL_HPP
#define CONTRACTIONS_SU2_SITE_HAMIL_HPP

namespace contraction {

    // forward declarations
    template<class Matrix, class OtherMatrix, class SymmGroup>
    MPSTensor<Matrix, SymmGroup>
    site_hamil_lbtm(MPSTensor<Matrix, SymmGroup> ket_tensor,
                    Boundary<OtherMatrix, SymmGroup> const & left,
                    Boundary<OtherMatrix, SymmGroup> const & right,
                    MPOTensor<Matrix, SymmGroup> const & mpo);

    template<class Matrix, class OtherMatrix, class SymmGroup>
    MPSTensor<Matrix, SymmGroup>
    site_hamil_rbtm(MPSTensor<Matrix, SymmGroup> ket_tensor,
                    Boundary<OtherMatrix, SymmGroup> const & left,
                    Boundary<OtherMatrix, SymmGroup> const & right,
                    MPOTensor<Matrix, SymmGroup> const & mpo);
    // *************************************************************


    template<class Matrix, class OtherMatrix, class SymmGroup>
    MPSTensor<Matrix, SymmGroup>
    Engine<Matrix, OtherMatrix, SymmGroup, typename boost::enable_if<symm_traits::HasSU2<SymmGroup> >::type>::
    site_hamil2(MPSTensor<Matrix, SymmGroup> ket_tensor,
                Boundary<OtherMatrix, SymmGroup> const & left,
                Boundary<OtherMatrix, SymmGroup> const & right,
                MPOTensor<Matrix, SymmGroup> const & mpo)
    {
        return site_hamil_lbtm(ket_tensor, left, right, mpo);
    }

    // *************************************************************
    // specialized variants

    template<class Matrix, class OtherMatrix, class SymmGroup>
    MPSTensor<Matrix, SymmGroup>
    site_hamil_lbtm(MPSTensor<Matrix, SymmGroup> ket_tensor,
                    Boundary<OtherMatrix, SymmGroup> const & left,
                    Boundary<OtherMatrix, SymmGroup> const & right,
                    MPOTensor<Matrix, SymmGroup> const & mpo)
    {
        typedef typename SymmGroup::charge charge;
        typedef typename MPOTensor<Matrix, SymmGroup>::index_type index_type;

        //std::vector<block_matrix<Matrix, SymmGroup> > t
        BoundaryMPSProduct<Matrix, OtherMatrix, SymmGroup, ::SU2::SU2Gemms> t
            = common::boundary_times_mps<Matrix, OtherMatrix, SymmGroup, ::SU2::SU2Gemms>(ket_tensor, left, mpo);

        Index<SymmGroup> const & physical_i = ket_tensor.site_dim(),
                               & left_i = ket_tensor.row_dim();
        Index<SymmGroup> right_i = ket_tensor.col_dim(),
                         out_left_i = physical_i * left_i;

        common_subset(out_left_i, right_i);
        ProductBasis<SymmGroup> out_left_pb(physical_i, left_i);
        ProductBasis<SymmGroup> in_right_pb(physical_i, right_i,
                                boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                        -boost::lambda::_1, boost::lambda::_2));

        MPSTensor<Matrix, SymmGroup> ret;
        ret.phys_i = ket_tensor.site_dim(); ret.left_i = ket_tensor.row_dim(); ret.right_i = ket_tensor.col_dim();
        index_type loop_max = mpo.col_dim();

#ifdef USE_AMBIENT
        {
            block_matrix<Matrix, SymmGroup> empty;
            swap(ket_tensor.data(), empty); // deallocating mpstensor before exiting the stack
        }
        parallel::sync();
        ContractionGrid<Matrix, SymmGroup> contr_grid(mpo, left.aux_dim(), mpo.col_dim());

        parallel_for(index_type b2, parallel::range<index_type>(0,loop_max), {
            SU2::lbtm_kernel(b2, contr_grid, left, t, mpo, ket_tensor.data().basis(), right_i, out_left_i, in_right_pb, out_left_pb);
        });
        omp_for(index_type b2, parallel::range<index_type>(0,loop_max), {
            contr_grid.multiply_column(b2, right[b2]);
        });
        t.clear();
        parallel::sync();

        swap(ret.data(), contr_grid.reduce());
        parallel::sync();

#else
        omp_for(index_type b2, parallel::range<index_type>(0,loop_max), {
            ContractionGrid<Matrix, SymmGroup> contr_grid(mpo, 0, 0);
            block_matrix<Matrix, SymmGroup> tmp;

            typename MPOTensor<OtherMatrix, SymmGroup>::col_proxy cp = mpo.column(b2);
            index_type num_ops = std::distance(cp.begin(), cp.end());
            if (num_ops > 3) {
                SU2::lbtm_kernel_rp(b2, contr_grid, left, t, mpo, ket_tensor.data().basis(), right_i, out_left_i, in_right_pb, out_left_pb, right[b2].basis());
                block_matrix<Matrix, SymmGroup> tmp2;
                reshape_right_to_left_new(physical_i, left_i, right_i, contr_grid(0,0), tmp2);

                contr_grid(0,0).clear();
                ::SU2::gemm_trim(tmp2, right[b2], tmp);

                for (std::size_t k = 0; k < tmp.n_blocks(); ++k)
                    if (!out_left_i.has(tmp.basis().left_charge(k)))
                        tmp.remove_block(k--);
            }

            else {
                SU2::lbtm_kernel(b2, contr_grid, left, t, mpo, ket_tensor.data().basis(), right_i, out_left_i, in_right_pb, out_left_pb);
                ::SU2::gemm_trim(contr_grid(0,0), right[b2], tmp);

                contr_grid(0,0).clear();
            }

            parallel_critical
            for (std::size_t k = 0; k < tmp.n_blocks(); ++k)
                ret.data().match_and_add_block(tmp[k], tmp.basis().left_charge(k), tmp.basis().right_charge(k));
        });
#endif
        return ret;
    }

    template<class Matrix, class OtherMatrix, class SymmGroup>
    MPSTensor<Matrix, SymmGroup>
    site_hamil_rbtm(MPSTensor<Matrix, SymmGroup> ket_tensor,
                    Boundary<OtherMatrix, SymmGroup> const & left,
                    Boundary<OtherMatrix, SymmGroup> const & right,
                    MPOTensor<Matrix, SymmGroup> const & mpo)
    {
        typedef typename SymmGroup::charge charge;
        typedef typename MPOTensor<Matrix, SymmGroup>::index_type index_type;

        //std::vector<block_matrix<Matrix, SymmGroup> > t
        MPSBoundaryProduct<Matrix, OtherMatrix, SymmGroup, ::SU2::SU2Gemms> t
            = common::mps_times_boundary<Matrix, OtherMatrix, SymmGroup, ::SU2::SU2Gemms>(ket_tensor, right, mpo);

        Index<SymmGroup> const & physical_i = ket_tensor.site_dim(),
                                 right_i = ket_tensor.col_dim();
        Index<SymmGroup> left_i = ket_tensor.row_dim(),
                         out_right_i = adjoin(physical_i) * right_i;

        common_subset(out_right_i, left_i);
        ProductBasis<SymmGroup> in_left_pb(physical_i, left_i);
        ProductBasis<SymmGroup> out_right_pb(physical_i, right_i,
                                             boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                                 -boost::lambda::_1, boost::lambda::_2));
        block_matrix<Matrix, SymmGroup> collector;
        MPSTensor<Matrix, SymmGroup> ret;
        ret.phys_i = ket_tensor.site_dim(); ret.left_i = ket_tensor.row_dim(); ret.right_i = ket_tensor.col_dim();

        index_type loop_max = mpo.row_dim();
        omp_for(index_type b1, parallel::range<index_type>(0,loop_max), {

            block_matrix<Matrix, SymmGroup> tmp, tmp2;
            SU2::rbtm_kernel(b1, tmp, right, t, mpo, ket_tensor.data().basis(), left_i, out_right_i, in_left_pb, out_right_pb);

            ::SU2::gemm_trim(transpose(left[b1]), tmp, tmp2);
            tmp.clear();

            parallel_critical
            for (std::size_t k = 0; k < tmp2.n_blocks(); ++k)
                collector.match_and_add_block(tmp2[k], tmp2.basis().left_charge(k), tmp2.basis().right_charge(k));
        });

        reshape_right_to_left_new(physical_i, left_i, right_i, collector, ret.data());
        return ret;
    }

} // namespace contraction

#endif
