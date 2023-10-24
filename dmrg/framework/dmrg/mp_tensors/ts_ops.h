/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MAQUIS_DMRG_TS_OPS_H
#define MAQUIS_DMRG_TS_OPS_H

#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/mpotensor.h"
#include "dmrg/block_matrix/indexing.h"
#include "dmrg/block_matrix/multi_index.h"
#include "dmrg/block_matrix/site_operator.h"
#include "dmrg/block_matrix/site_operator_algorithms.h"
#include "dmrg/models/OperatorHandlers/KronHandler.h"

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
        using spin_t = typename SymmGroup::subcharge;
        using op_t = typename OPTable<Matrix, SymmGroup>::op_t;
        std::map<spin_t, typename OPTable<Matrix, SymmGroup>::op_t> ret;
        for (typename std::set<Integer>::const_iterator it1 = summands.begin(); it1 != summands.end(); ++it1)
        {
            auto b2 = *it1;
            term_descriptor<Matrix, SymmGroup, true> p1 = mpo1.at(b1,b2), p2 = mpo2.at(b2,b3);
            for (int ip1 = 0; ip1 < p1.size(); ++ip1) {
                for (int ip2 = 0; ip2 < p2.size(); ++ip2) {
                    std::vector<spin_t> op_spins = allowed_spins(mpo1.left_spin(b1).get(), mpo2.right_spin(b3).get(),
                                                                 p1.op(ip1).spin().get(), p2.op(ip2).spin().get());
                    for (auto it2 = op_spins.begin(); it2 != op_spins.end(); ++it2) {
                        SpinDescriptor<typename symm_traits::SymmType<SymmGroup>::type>
                            prod_spin(*it2, mpo1.left_spin(b1).get(), mpo2.right_spin(b3).get());
                        op_t product;
                        op_kron(phys_i1, phys_i2, p1.op(ip1), p2.op(ip2), product, mpo1.left_spin(b1),
                                mpo1.right_spin(b2), mpo2.right_spin(b3), prod_spin);
                        ::tag_detail::remove_empty_blocks(product);
                        ret[*it2] += product * p1.scale(ip1) * p2.scale(ip2);
                    }
                }
            }
        }
        return ret;
    }

} // namespace ts_ops_detail

/**
 * @brief Combines two (neighbouring) MPOs to generate a two-site MPO tensor.
 * @param mpo1 First MPO
 * @param mpo2 Second MPO
 * @param phys_i1 physical basis for the first MPO
 * @param phys_i2 physical basis for the second MPO
 * @return MPOTensor<MPSMatrix, SymmGroup> two-site MPO
 */
template <class MPOMatrix, class MPSMatrix, class SymmGroup>
MPOTensor<MPSMatrix, SymmGroup> make_twosite_mpo(MPOTensor<MPOMatrix, SymmGroup> const & mpo1,
                                                 MPOTensor<MPOMatrix, SymmGroup> const & mpo2,
                                                 Index<SymmGroup> const & phys_i1,
                                                 Index<SymmGroup> const & phys_i2)
{
    using index_type = typename MPOTensor<MPOMatrix, SymmGroup>::index_type;
    using row_proxy = typename MPOTensor<MPOMatrix, SymmGroup>::row_proxy;
    using tag_type = typename OPTable<MPOMatrix, SymmGroup>::tag_type;
    using op_t = typename OPTable<MPOMatrix, SymmGroup>::op_t;
    using value_type = typename MPSMatrix::value_type;
    using prempo_t = std::vector<boost::tuple<index_type, index_type, tag_type, value_type> >;
    //
    using MPOTensor_detail::term_descriptor;
    using boost::tuples::get;
    // Variable declaration
    assert(mpo1.col_dim() == mpo2.row_dim());
    KronHandler<MPOMatrix, SymmGroup> kron_handler(mpo1.get_operator_table());
    prempo_t prempo;
    index_type b1, b2, b3;
    // Main loop
    for (b1 = 0; b1 < mpo1.row_dim(); ++b1) {
        for (b3 = 0; b3 < mpo2.col_dim(); ++b3) {
            // This proxy tells how many b2 values are coupled to a given b1 value
            row_proxy row1 = mpo1.row(b1);
            op_t b3_op;
            std::set<index_type> summands;
            // Checks that the elements to which b1 is coupled are also coupled to b3
            for (typename row_proxy::const_iterator it = row1.begin(); it != row1.end(); ++it)
                if (mpo2.has(it.index(), b3))
                    summands.insert(it.index());
            // Evaluates the sum over b2 and returns the (b1, b3) element of the MPOTensor 
            auto coupled_ops = ts_ops_detail::mpo_couple(summands, b1, b3, phys_i1, phys_i2, mpo1, mpo2);
            for (auto it = coupled_ops.begin(); it != coupled_ops.end(); ++it) {
                tag_type new_tag = kron_handler.get_kronecker_table()->register_op(it->second);
                prempo.push_back(boost::make_tuple(b1, b3, new_tag, 1.0));
            }
        }
    }
    // Creates the final MPO
    MPOTensor<MPSMatrix, SymmGroup> mpo_big_tag(mpo1.row_dim(), mpo2.col_dim(), prempo, kron_handler.get_kronecker_table(),
                                                mpo1.herm_info * mpo2.herm_info, mpo1.row_spin_dim(), mpo2.col_spin_dim());
    // #ifdef MAQUIS_OPENMP
    // #pragma omp critical
    // #endif
    // maquis::cout << "TSMPOTensor: " << mpo1.row_dim() << "x" << mpo2.col_dim() << ",  " << prempo.size()
    //              << " operators, " << kron_handler.get_kronecker_table()->size() << " tags\n";
    return mpo_big_tag;
}

/**
 * @brief Generates the TwoSite MPO associated with a given MPO object.
 *
 * Note that the (iSite)-th element of the output MPO is the two-site MPO for sites (iSite, iSite+1)
 *
 * @param mpo_orig Input MPO
 * @param mpo_out Output MPO
 * @param mps Reference MPS (to get the physical dimension)
 */
template<class MPOMatrix, class MPSMatrix, class SymmGroup>
void make_ts_cache_mpo(MPO<MPOMatrix, SymmGroup> const & mpo_orig,
                       MPO<MPSMatrix, SymmGroup> & mpo_out, MPS<MPSMatrix, SymmGroup> const & mps)
{
    auto L_ts = mpo_orig.length() - 1;
    mpo_out.resize(L_ts);
    // Generates the two-site MPOs
    omp_for(size_t p, parallel::range<size_t>(0,L_ts), {
        mpo_out[p] = make_twosite_mpo<MPOMatrix, MPSMatrix>(mpo_orig[p], mpo_orig[p+1], mps[p].site_dim(), mps[p+1].site_dim());
    });
    // Calculates the overall number of tags
    std::size_t ntags=0;
    for (int p=0; p<mpo_out.length(); ++p)
        ntags += mpo_out[p].get_operator_table()->size();
    // maquis::cout << "Total number of tags: " << ntags << std::endl;
}

#endif
