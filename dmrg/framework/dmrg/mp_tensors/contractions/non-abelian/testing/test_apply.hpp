/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
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

#ifndef CONTRACTIONS_SU2_TESTING_APPLY_HPP
#define CONTRACTIONS_SU2_TESTING_APPLY_HPP

#include <boost/tuple/tuple_io.hpp>

#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/mpotensor.h"

namespace SU2 { // Forward declarations

    double mod_coupling(int two_ja, int two_jb, int two_jc,
                        int two_jd, int two_je, int two_jf,
                        int two_jg, int two_jh, int two_ji);

    template<class Matrix1, class Matrix2, class Matrix3, class SymmGroup>
    void gemm(block_matrix<Matrix1, SymmGroup> const & A,
              block_matrix<Matrix2, SymmGroup> const & B,
              block_matrix<Matrix3, SymmGroup> & C);
}

namespace contraction {
namespace SU2 {

    template <class SymmGroup>
    bool column_check(typename SymmGroup::charge c, Index<SymmGroup> const & ind)
    {
        int count = 0;
        for (std::size_t p = 0; p < ind.size(); ++p)
            if (ind[p].first == c) ++count;

        return (count == 1);
    }

    template <typename T> int sgn(T val) {
        return (T(0) < val) - (val < T(0));
    }

    template<class Matrix, class OtherMatrix, class SymmGroup>
    block_matrix<OtherMatrix, SymmGroup>
    apply_operator(MPSTensor<Matrix, SymmGroup> const & bra_tensor,
                   MPSTensor<Matrix, SymmGroup> const & ket_tensor,
                   block_matrix<OtherMatrix, SymmGroup> const & left,
                   MPOTensor<Matrix, SymmGroup> const & mpo,
                   std::vector<int> const & config,
                   bool debug = false)
    {
        typedef typename SymmGroup::charge charge;
        typedef typename MPOTensor<Matrix, SymmGroup>::index_type index_type;
        typedef typename DualIndex<SymmGroup>::const_iterator const_iterator;
        typedef typename Index<SymmGroup>::const_iterator ci;

        assert(ket_tensor.phys_i == bra_tensor.phys_i);

        bra_tensor.make_left_paired();
        ket_tensor.make_right_paired();

        MPOTensor_detail::const_term_descriptor<Matrix, SymmGroup> access = mpo.at(0,0);
        block_matrix<Matrix, SymmGroup> const & W = access.op;
        block_matrix<OtherMatrix, SymmGroup> ret;

        block_matrix<OtherMatrix, SymmGroup> t1;
        ::SU2::gemm(left, ket_tensor.data(), t1);


        Index<SymmGroup> const & left_i = ket_tensor.row_dim();
        Index<SymmGroup> const & right_i = ket_tensor.col_dim();
        Index<SymmGroup> const & phys_i = ket_tensor.site_dim();

        ProductBasis<SymmGroup> out_left_pb(phys_i, left_i);
        ProductBasis<SymmGroup> in_right_pb(ket_tensor.site_dim(), right_i,
                                boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                        -boost::lambda::_1, boost::lambda::_2));

        //maquis::cout << "BraLeft: " << bra_tensor.data().left_basis() << std::endl;
        //maquis::cout << "outLeft: " << (phys_i*left_i) << std::endl;
        for (std::size_t k = 0; k < left.n_blocks(); ++k) {

            charge lc = left.basis().lc(k);
            charge mc = left.basis().rc(k);

            std::pair<const_iterator, const_iterator>
              er = std::equal_range(ket_tensor.data().basis().begin(), ket_tensor.data().basis().end(),
                boost::make_tuple(mc, SymmGroup::IdentityCharge, 0, 0), dual_index_detail::gt_row<SymmGroup>());

            assert(er.first + 1 == er.second);
            for (const_iterator it = er.first; it != er.second; ++it)
            {
                std::size_t matched_block = std::distance(ket_tensor.data().basis().begin(), it);
                //Matrix T_block(num_rows(left[k]), num_cols(ket_tensor[matched_block]));
                //gemm(A[k], B[matched_block], T_block);
                charge rc = ket_tensor.data().basis().rc(matched_block);
                charge mc1 = ket_tensor.data().basis().lc(matched_block);
                assert (mc == mc1);
                size_t t_pos = t1.basis().position(lc, rc);
                assert(t_pos == k);

                for (size_t w_block = 0; w_block < W.basis().size(); ++w_block)
                {
                    charge phys_in = W.basis().lc(w_block);
                    charge phys_out = W.basis().rc(w_block);

                    charge free_rc = SymmGroup::fuse(rc, phys_in);
                    if (!right_i.has(free_rc))
                        continue;

                    charge new_rc = SymmGroup::fuse(lc, phys_out);
                    if (!bra_tensor.col_dim().has(new_rc))
                        continue;

                    int i  = lc[1], ip = new_rc[1];
                    int j  = mc[1], jp  = free_rc[1];
                    int two_sp = std::abs(i - ip), two_s  = std::abs(j - jp);
                    int a = std::abs(i-j), k = std::abs(std::abs(phys_in[1])-std::abs(phys_out[1]));

                    for (int ap = std::abs(a-k); ap <= std::abs(a+k); ap+=2)
                    {
                    if (ap != std::abs(ip-jp)) continue;
                    if (ap >= 3) continue;

                    double coupling_coeff = ::SU2::mod_coupling(j, two_s, jp, a,k,ap, i, two_sp, ip);
                    coupling_coeff *= pow(ip+1., 0.5) * pow(j+1., 0.5);
                    coupling_coeff *= pow(i+1., -0.5) * pow(jp+1., -0.5);
                    coupling_coeff *= access.scale * W[w_block](0,0);

                    if (debug) {
                        std::cout << j << "," << two_s << "," << jp << " | " << a << "," << k << "," << ap << " | "
                                  << i << "," << two_sp << "," << ip << " | " << phys_in << phys_out
                                  << "  " << new_rc[0] << "," << lc[0] << "," << mc[0] << "," << free_rc[0]
                                  << std::right << std::setw(8) << "cc: " << std::setw(12) << coupling_coeff << std::endl;
                    }

                    // T Access
                    size_t right_offset = in_right_pb(phys_in, free_rc);
                    size_t left_offset = out_left_pb(phys_out, lc);
                    size_t ldim = t1.basis().ls(t_pos);
                    size_t rdim = right_i.size_of_block(free_rc);
                    //maquis::cout << "access " << mc << " + " << phys_in<< "|" << free_rc << " -- "
                    //    << lc << " + " << phys_out << "|" << new_rc
                    //    << " T(" << lc << "," << rc << "): +" << right_offset
                    //    << "(" << ldim << "x" << rdim << ")|" << t1.basis().rs(t_pos) << " -> +"
                    //    << left_offset << "(" << t1.basis().ls(t_pos) << "x" << rdim << ")|" << bra_tensor.data().left_basis().size_of_block(new_rc)
                    //    << " * " << j << two_s << jp << a << k << ap << i << two_sp << ip << std::endl;
                    Matrix T_cp(ldim, rdim, 0);
                    for (size_t c=0; c<rdim; ++c)
                        std::transform(t1[t_pos].col(right_offset+c).first,
                                       t1[t_pos].col(right_offset+c).second, T_cp.col(c).first, boost::lambda::_1*coupling_coeff);

                    // Bra Access
                    size_t bra_index = bra_tensor.col_dim().position(new_rc);
                    assert(column_check(new_rc, bra_tensor.col_dim()));
                    Matrix const & bra_block = bra_tensor.data()[bra_index];
                    // mc + phys_out = new_rc
                    //size_t left_offset = out_left_pb(phys_out, lc);
                    ldim = bra_tensor.row_dim().size_of_block(lc);
                    rdim = bra_tensor.col_dim()[bra_index].second;
                    Matrix bra_cp(ldim, rdim, 0);
                    for (size_t row=0; row<ldim; ++row)
                        std::copy(bra_block.row(left_offset+row).first,
                                  bra_block.row(left_offset+row).second, bra_cp.row(row).first);

                    // Multiply
                    Matrix prod(num_cols(bra_cp), num_cols(T_cp), 0);
                    gemm(transpose(bra_cp), T_cp, prod);

                    // Check-in
                    ret.match_and_add_block(prod, new_rc, free_rc);
                    }
                }
            }
        }
        if(debug) maquis::cout << std::endl;
        return ret;
    }
} // namespace SU2
} // namespace contractions

#endif
