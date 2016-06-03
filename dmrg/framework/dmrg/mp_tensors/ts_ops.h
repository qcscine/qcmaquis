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
#include "dmrg/block_matrix/site_operator.h"
#include "dmrg/block_matrix/site_operator_algorithms.h"


namespace ts_ops_detail
{
    template <class Integer>
    std::vector<Integer> allowed_spins(Integer left, Integer right, Integer k1, Integer k2)
    {
        std::vector<Integer> operator_spins;
        for (Integer s = std::abs(k1-k2); s <= std::abs(k1+k2); s+=2)
            operator_spins.push_back(s);

        // triangle condition for the operator action on input/output spins
        for (typename std::vector<Integer>::iterator it = operator_spins.begin(); it != operator_spins.end(); ++it)
            if ( !(right >= std::abs(*it-left)) || !(right <= std::abs(*it+left)) )
                operator_spins.erase(it--);

        return operator_spins;
    }

    template <class Integer, class Matrix, class SymmGroup>
    std::map<typename SymmGroup::subcharge, typename OPTable<Matrix, SymmGroup>::op_t>
    mpo_couple(std::set<Integer> const & summands, Integer b1, Integer b3, Index<SymmGroup> const & phys_i1, Index<SymmGroup> const & phys_i2,
               MPOTensor<Matrix, SymmGroup> const & mpo1, MPOTensor<Matrix, SymmGroup> const & mpo2)
    {
        using MPOTensor_detail::term_descriptor;
        typedef typename SymmGroup::subcharge spin_t;
        typedef typename OPTable<Matrix, SymmGroup>::op_t op_t;

        std::map<spin_t, typename OPTable<Matrix, SymmGroup>::op_t> ret;

        for (typename std::set<Integer>::const_iterator it1 = summands.begin(); it1 != summands.end(); ++it1) {
            Integer b2 = *it1;
            term_descriptor<Matrix, SymmGroup, true> p1 = mpo1.at(b1,b2), p2 = mpo2.at(b2,b3);

            std::vector<spin_t> op_spins = allowed_spins(mpo1.left_spin(b1).get(), mpo2.right_spin(b3).get(), p1.op().spin().get(), p2.op().spin().get());
            for (typename std::vector<spin_t>::const_iterator it2 = op_spins.begin(); it2 != op_spins.end(); ++it2)
            {
                SpinDescriptor<typename symm_traits::SymmType<SymmGroup>::type> prod_spin(*it2, mpo1.left_spin(b1).get(), mpo2.right_spin(b3).get());

                op_t product;
                op_kron(phys_i1, phys_i2, p1.op(), p2.op(), product, mpo1.left_spin(b1), mpo1.right_spin(b2), mpo2.right_spin(b3), prod_spin);
                ::tag_detail::remove_empty_blocks(product);
                ret[*it2] += product * p1.scale() * p2.scale();
            }
        }

        return ret;
    }

} // namespace ts_ops_detail

template <class MPOMatrix, class MPSMatrix, class SymmGroup>
MPOTensor<MPSMatrix, SymmGroup> make_twosite_mpo(MPOTensor<MPOMatrix, SymmGroup> const & mpo1,
                                                 MPOTensor<MPOMatrix, SymmGroup> const & mpo2,
                                                 Index<SymmGroup> const & phys_i1,
                                                 Index<SymmGroup> const & phys_i2)
{
    using MPOTensor_detail::term_descriptor;
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
            std::set<index_type> summands;

            for (typename row_proxy::const_iterator it = row1.begin(); it != row1.end(); ++it)
                if (mpo2.has(it.index(), b3))
                    summands.insert(it.index());

            if (summands.size() > 1 || (!shared && summands.size() > 0))
            {
                std::map<typename SymmGroup::subcharge, op_t> coupled_ops
                    = ts_ops_detail::mpo_couple(summands, b1, b3, phys_i1, phys_i2, mpo1, mpo2);

                for (typename std::map<typename SymmGroup::subcharge, op_t>::const_iterator it = coupled_ops.begin();
                        it != coupled_ops.end(); ++it)
                {
                    tag_type new_tag = kron_handler.get_kronecker_table()->register_op(it->second);
                    prempo.push_back(boost::make_tuple(b1, b3, new_tag, 1.0));
                }
            }
            else if (summands.size() == 1)
            {
                index_type b2 = *summands.begin();
                term_descriptor<MPOMatrix, SymmGroup, true> p1 = mpo1.at(b1,b2), p2 = mpo2.at(b2,b3);
                tag_type p_tag = kron_handler.get_kron_tag(phys_i1, phys_i2, mpo1.tag_number(b1,b2), mpo2.tag_number(b2,b3),
                                                             mpo1.left_spin(b1), mpo1.right_spin(b2), mpo2.right_spin(b3));
                value_type p_scale = p1.scale() * p2.scale();
                prempo.push_back(boost::make_tuple(b1, b3, p_tag, p_scale));
            }

        } // b3
    } // b1

    MPOTensor<MPSMatrix, SymmGroup> mpo_big_tag(mpo1.row_dim(), mpo2.col_dim(), prempo, kron_handler.get_kronecker_table(),
                                                mpo1.herm_info * mpo2.herm_info, mpo1.row_spin_dim(), mpo2.col_spin_dim());
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
