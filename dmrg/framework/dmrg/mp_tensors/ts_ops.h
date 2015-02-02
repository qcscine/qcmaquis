/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

template<class MPOMatrix, class MPSMatrix, class SymmGroup>
MPOTensor<MPSMatrix, SymmGroup> make_twosite_mpo_shared(MPOTensor<MPOMatrix, SymmGroup> const & mpo1,
                                                        MPOTensor<MPOMatrix, SymmGroup> const & mpo2,
                                                        Index<SymmGroup> const & phys_i1,
                                                        Index<SymmGroup> const & phys_i2)
{
    using MPOTensor_detail::term_descriptor;
    using boost::tuples::get;
    assert(mpo1.col_dim() == mpo2.row_dim());

    typedef typename MPOTensor<MPOMatrix, SymmGroup>::index_type index_type;
    typedef typename MPOTensor<MPOMatrix, SymmGroup>::row_proxy row_proxy;
    typedef typename OPTable<MPSMatrix, SymmGroup>::tag_type tag_type;
    typedef typename OPTable<MPSMatrix, SymmGroup>::op_t op_t;
    typedef typename MPSMatrix::value_type value_type;
    typedef std::vector<boost::tuple<index_type, index_type, tag_type, value_type> > prempo_t;

    // Use a separate operator table for every twosite MPO Tensor - too much overhead to globally synchronize
    // a single table
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

            if (summands.size() > 1)
            {
                for (typename std::set<index_type>::const_iterator it = summands.begin(); it != summands.end(); ++it) {
                    index_type b2 = *it; 
                    term_descriptor<MPSMatrix, SymmGroup, true> p1 = mpo1.at(b1,b2), p2 = mpo2.at(b2,b3);

                    op_t product;
                    op_kron(phys_i1, phys_i2, p1.op(), p2.op(), product, mpo1.left_spin(b1), mpo1.right_spin(b2), mpo2.right_spin(b3));
                    b3_op += product * p1.scale() * p2.scale();
                }
                tag_detail::remove_empty_blocks(b3_op);
                scaled_tag = kron_handler.get_kronecker_table()->checked_register(b3_op);
                prempo.push_back(boost::make_tuple(b1, b3, scaled_tag.first, scaled_tag.second));
            }
            else if (summands.size() == 1)
            {
                index_type b2 = *summands.begin();
                term_descriptor<MPSMatrix, SymmGroup, true> p1 = mpo1.at(b1,b2), p2 = mpo2.at(b2,b3);
                scaled_tag.first = kron_handler.get_kron_tag(phys_i1, phys_i2, mpo1.tag_number(b1,b2), mpo2.tag_number(b2,b3),
                                                             mpo1.left_spin(b1), mpo1.right_spin(b2), mpo2.right_spin(b3));
                scaled_tag.second = p1.scale() * p2.scale();
                prempo.push_back(boost::make_tuple(b1, b3, scaled_tag.first, scaled_tag.second));
            }

        } // b3
    } // b1

    MPOTensor<MPSMatrix, SymmGroup> mpo_big_tag(mpo1.row_dim(), mpo2.col_dim(), prempo, kron_handler.get_kronecker_table(),
                                                mpo1.row_spin_dim(), mpo2.col_spin_dim());

    #ifdef MAQUIS_OPENMP
    #pragma omp critical
    #endif
    maquis::cout << "TSMPOTensor: " << mpo1.row_dim() << "x" << mpo2.col_dim() << ",  " << prempo.size() 
                 << " operators, " << kron_handler.get_kronecker_table()->size() << " tags\n";
    return mpo_big_tag;
}

template<class MPOMatrix, class MPSMatrix, class SymmGroup>
MPOTensor<MPSMatrix, SymmGroup> make_twosite_mpo(MPOTensor<MPOMatrix, SymmGroup> const & mpo1,
                                                 MPOTensor<MPOMatrix, SymmGroup> const & mpo2,
                                                 Index<SymmGroup> const & phys_i1,
                                                 Index<SymmGroup> const & phys_i2)
{
    using MPOTensor_detail::term_descriptor;
    using boost::tuples::get;
    assert(mpo1.col_dim() == mpo2.row_dim());

    typedef typename MPOTensor<MPOMatrix, SymmGroup>::index_type index_type;
    typedef typename MPOTensor<MPOMatrix, SymmGroup>::row_proxy row_proxy;
    typedef typename OPTable<MPSMatrix, SymmGroup>::tag_type tag_type;
    typedef typename OPTable<MPSMatrix, SymmGroup>::op_t op_t;
    typedef typename MPSMatrix::value_type value_type;
    typedef std::vector<boost::tuple<index_type, index_type, tag_type, value_type> > prempo_t;

    prempo_t prempo;
    typename MPOTensor<MPSMatrix, SymmGroup>::op_table_ptr op_table(new OPTable<MPSMatrix, SymmGroup>);

    index_type b1, b2, b3;
    for (b1 = 0; b1 < mpo1.row_dim(); ++b1) {
        for (b3 = 0; b3 < mpo2.col_dim(); ++b3) {

            row_proxy row1 = mpo1.row(b1);

            op_t b3_op;
            std::pair<tag_type, value_type> scaled_tag;
            bool commit = false;

            for (typename row_proxy::const_iterator it = row1.begin(); it != row1.end(); ++it) {
                index_type b2 = it.index();
                if (mpo2.has(b2, b3))
                {
                    commit = true;
                    term_descriptor<MPSMatrix, SymmGroup, true> p1 = mpo1.at(b1,b2), p2 = mpo2.at(b2,b3);

                    op_t product;
                    op_kron(phys_i1, phys_i2, p1.op(), p2.op(), product, mpo1.left_spin(b1), mpo1.right_spin(b2), mpo2.right_spin(b3));
                    b3_op += product * p1.scale() * p2.scale();
                }
            }
            if (commit) {
            tag_detail::remove_empty_blocks(b3_op);
            tag_type new_tag = op_table->register_op(b3_op);
            prempo.push_back(boost::make_tuple(b1, b3, new_tag, 1.0));
            }

        } // b3
    } // b1

    MPOTensor<MPSMatrix, SymmGroup> mpo_big(mpo1.row_dim(), mpo2.col_dim(), prempo, op_table, mpo1.row_spin_dim(), mpo2.col_spin_dim());

    #ifdef MAQUIS_OPENMP
    #pragma omp critical
    #endif
    maquis::cout << "TSMPOTensor: " << mpo1.row_dim() << "x" << mpo2.col_dim() << ",  " << prempo.size() 
                 << " operators, " << op_table->size() << " tags\n";
    return mpo_big;
}

template<class MPOMatrix, class MPSMatrix, class SymmGroup>
void make_ts_cache_mpo(MPO<MPOMatrix, SymmGroup> const & mpo_orig,
                       MPO<MPSMatrix, SymmGroup> & mpo_out, MPS<MPSMatrix, SymmGroup> const & mps)
{
    std::size_t L_ts = mpo_orig.length() - 1;
    mpo_out.resize(L_ts);

    omp_for(size_t p, parallel::range<size_t>(0,L_ts), {
        if (mpo_orig[p].get_operator_table() == mpo_orig[p+1].get_operator_table())
            mpo_out[p] = make_twosite_mpo_shared<MPOMatrix, MPSMatrix>(mpo_orig[p], mpo_orig[p+1], mps[p].site_dim(), mps[p+1].site_dim());
        else
            mpo_out[p] = make_twosite_mpo<MPOMatrix, MPSMatrix>(mpo_orig[p], mpo_orig[p+1], mps[p].site_dim(), mps[p+1].site_dim());
    });
        
    std::size_t ntags=0;
    for (int p=0; p<mpo_out.length(); ++p) {
        ntags += mpo_out[p].get_operator_table()->size();
    }
    maquis::cout << "Total number of tags: " << ntags << std::endl;
}

#endif
