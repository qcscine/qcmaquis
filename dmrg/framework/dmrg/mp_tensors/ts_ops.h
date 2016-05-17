/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2015 Institute for Theoretical Physics, ETH Zurich
 *                    Laboratory for Physical Chemistry, ETH Zurich
 *               2013-2015 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *
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

#ifndef MAQUIS_DMRG_TS_OPS_H
#define MAQUIS_DMRG_TS_OPS_H

#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/mpotensor.h"
#include "dmrg/block_matrix/indexing.h"
#include "dmrg/block_matrix/multi_index.h"
#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"

template <class MPOMatrix, class MPSMatrix, class SymmGroup>
MPOTensor<MPSMatrix, SymmGroup> make_twosite_mpo(MPOTensor<MPOMatrix, SymmGroup> const & mpo1,
                                                 MPOTensor<MPOMatrix, SymmGroup> const & mpo2,
                                                 Index<SymmGroup> const & phys_i1,
                                                 Index<SymmGroup> const & phys_i2)
{
    using MPOTensor_detail::const_term_descriptor;
    using boost::tuples::get;
    assert(mpo1.col_dim() == mpo2.row_dim());
    bool shared = (mpo1.get_operator_table() == mpo2.get_operator_table());

    typedef typename MPOTensor<MPOMatrix, SymmGroup>::index_type index_type;
    typedef typename MPOTensor<MPOMatrix, SymmGroup>::row_proxy row_proxy;
    typedef typename OPTable<MPOMatrix, SymmGroup>::tag_type tag_type;
    typedef typename OPTable<MPOMatrix, SymmGroup>::op_t op_t;
    typedef typename MPSMatrix::value_type value_type;
    typedef std::vector<boost::tuple<index_type, index_type, tag_type, value_type> > prempo_t;

    KronHandler<MPOMatrix, SymmGroup> kron_handler(mpo1.get_operator_table());
    prempo_t prempo;

    index_type b1, b2, b3;
    for (b1 = 0; b1 < mpo1.row_dim(); ++b1) {
        for (b3 = 0; b3 < mpo2.col_dim(); ++b3) {

            row_proxy row1 = mpo1.row(b1);

            op_t b3_op;
            std::pair<tag_type, value_type> scaled_tag;
            std::set<index_type> summands;

            for (typename row_proxy::const_iterator it = row1.begin(); it != row1.end(); ++it)
                if (mpo2.has(it.index(), b3))
                    summands.insert(it.index());

            if (summands.size() > 1 || (!shared && summands.size() > 0))
            {
                for (typename std::set<index_type>::const_iterator it = summands.begin(); it != summands.end(); ++it) {
                    index_type b2 = *it; 
                    const_term_descriptor<MPOMatrix, SymmGroup> p1 = mpo1.at(b1,b2), p2 = mpo2.at(b2,b3);

                    op_t product;
                    op_kron(phys_i1, phys_i2, p1.op, p2.op, product);
                    b3_op += product * p1.scale * p2.scale;
                }
                tag_detail::remove_empty_blocks(b3_op);
                scaled_tag = kron_handler.get_kronecker_table()->checked_register(b3_op);
                prempo.push_back(boost::make_tuple(b1, b3, scaled_tag.first, scaled_tag.second));
            }
            else if (summands.size() == 1)
            {
                index_type b2 = *summands.begin();
                const_term_descriptor<MPOMatrix, SymmGroup> p1 = mpo1.at(b1,b2), p2 = mpo2.at(b2,b3);
                scaled_tag.first = kron_handler.get_kron_tag(phys_i1, phys_i2, mpo1.tag_number(b1,b2), mpo2.tag_number(b2,b3));
                scaled_tag.second = p1.scale * p2.scale;
                prempo.push_back(boost::make_tuple(b1, b3, scaled_tag.first, scaled_tag.second));
            }

        } // b3
    } // b1

    MPOTensor<MPSMatrix, SymmGroup> mpo_big_tag(mpo1.row_dim(), mpo2.col_dim(), prempo, kron_handler.get_kronecker_table(),
                                                mpo1.herm_info * mpo2.herm_info);
    #ifdef MAQUIS_OPENMP
    #pragma omp critical
    #endif
    maquis::cout << "TSMPOTensor: " << mpo1.row_dim() << "x" << mpo2.col_dim() << ",  " << prempo.size() 
                 << " operators, " << kron_handler.get_kronecker_table()->size() << " tags\n";
    return mpo_big_tag;
}

template<class MPOMatrix, class MPSMatrix, class SymmGroup>
void make_ts_cache_mpo(MPO<MPOMatrix, SymmGroup> const & mpo_orig,
                       MPO<MPSMatrix, SymmGroup> & mpo_out, MPS<MPSMatrix, SymmGroup> const & mps)
{
    std::size_t L_ts = mpo_orig.length() - 1;
    mpo_out.resize(L_ts);

    omp_for(size_t p, parallel::range<size_t>(0,L_ts), {
        mpo_out[p] = make_twosite_mpo<MPOMatrix, MPSMatrix>(mpo_orig[p], mpo_orig[p+1], mps[p].site_dim(), mps[p+1].site_dim());
    });
        
    std::size_t ntags=0;
    for (int p=0; p<mpo_out.length(); ++p) {
        ntags += mpo_out[p].get_operator_table()->size();
    }
    maquis::cout << "Total number of tags: " << ntags << std::endl;
}

#endif
