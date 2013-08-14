/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2013-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *
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
MPOTensor<MPSMatrix, SymmGroup> make_twosite_mpo(MPOTensor<MPOMatrix, SymmGroup> const & mpo1,
                                                 MPOTensor<MPOMatrix, SymmGroup> const & mpo2,
                                                 Index<SymmGroup> const & phys_i,
                                                 bool global_table)
{
    using MPOTensor_detail::const_term_descriptor;

    assert(mpo1.col_dim() == mpo2.row_dim());

    typedef typename MPOTensor<MPOMatrix, SymmGroup>::index_type index_type;
    typedef typename MPOTensor<MPOMatrix, SymmGroup>::row_proxy row_proxy;
    typedef typename MPOTensor<MPOMatrix, SymmGroup>::col_proxy col_proxy;
    typedef typename OPTable<MPSMatrix, SymmGroup>::tag_type tag_type;
    typedef typename OPTable<MPSMatrix, SymmGroup>::op_t op_t;
    typedef typename MPSMatrix::value_type value_type;
    typedef std::map<index_type, op_t> op_map;
    typedef std::vector<boost::tuple<index_type, index_type, tag_type, value_type> > prempo_t;


    if (global_table) {
        KronHandler<MPOMatrix, SymmGroup> kron_handler(mpo1.get_operator_table());

        typedef std::map<index_type, std::pair<tag_type, value_type> > op_scale_map;

        prempo_t prempo;

        typename MPOTensor<MPOMatrix, SymmGroup>::op_table_ptr op_table = kron_handler.get_operator_table();

        index_type b1, b2, b3;
        for (b1=0; b1 < mpo1.row_dim(); ++b1) {

            op_map out_row;
            std::set<index_type> non_uniform;
            op_scale_map uniform_ops;

            row_proxy row1 = mpo1.row(b1);
            for (typename row_proxy::const_iterator it1 = row1.begin(); it1 != row1.end(); ++it1) {

                b2 = it1.index();
                row_proxy row2 = mpo2.row(b2);
                for (typename row_proxy::const_iterator it2 = row2.begin(); it2 != row2.end(); ++it2)
                {
                    b3 = it2.index();
                    assert((mpo1.has(b1, b2) && mpo2.has(b2, b3)));
                    tag_type kron_tag;

                    const_term_descriptor<MPSMatrix, SymmGroup> p1 = mpo1.at(b1,b2), p2 = mpo2.at(b2,b3);

                    // Compute the Kronecker product
                    kron_tag = kron_handler.get_kron_tag(phys_i, mpo1.tag_number(b1,b2), mpo2.tag_number(b2,b3));

                    if (!kron_handler.is_uniform(mpo1.tag_number(b1,b2)) ||
                        !kron_handler.is_uniform(mpo2.tag_number(b2,b3)) ||
                        uniform_ops.count(b3) > 0)
                    {
                        non_uniform.insert(b3);
                    }
                    else {
                        uniform_ops[b3].first = kron_tag;
                        uniform_ops[b3].second = (p1.scale * p2.scale);
                    }

                    block_matrix<MPSMatrix, SymmGroup> tmp_op;
                    tmp_op = kron_handler.get_op(kron_tag);
                    tmp_op *= (p1.scale * p2.scale);
                    out_row[b3] += tmp_op;
                }
            }

            for (b3 = 0; b3 < mpo2.col_dim(); ++b3) {
                if (non_uniform.count(b3) == 0 && uniform_ops.count(b3) == 0)
                    continue;

                if (non_uniform.count(b3) > 0) {
                    std::pair<tag_type, value_type> scaled_tag;
                    scaled_tag = kron_handler.get_kronecker_table()->checked_register(out_row[b3]);
                    prempo.push_back(boost::make_tuple(b1, b3, scaled_tag.first, scaled_tag.second));
                }
                else {
                    prempo.push_back(boost::make_tuple(b1, b3, uniform_ops[b3].first, uniform_ops[b3].second));
                }
            }

            /*
            #ifdef MAQUIS_OPENMP
            #pragma omp critical
            #endif
            for (typename op_map::iterator it = out_row.begin(); it != out_row.end(); ++it) {
                b3 = it->first;
                op_t & tmp = it->second;
                std::pair<tag_type, value_type> scaled_tag = kron_handler.get_kronecker_table()->checked_register(tmp);
                prempo.push_back(boost::make_tuple(b1, b3, scaled_tag.first, scaled_tag.second));
            }
            */
        } 



        #ifdef MAQUIS_OPENMP
        #pragma omp critical
        #endif
        maquis::cout << "TSMPOTensor: " << mpo1.row_dim() << "x" << mpo2.col_dim() << ",  " << prempo.size() 
                     << " operators, " << kron_handler.get_kronecker_table()->size() << " tags\n";

        using boost::tuples::get;
        MPOTensor<MPSMatrix, SymmGroup> mpo_big_tag(mpo1.row_dim(), mpo2.col_dim(), prempo, kron_handler.get_kronecker_table());

        return mpo_big_tag;

    }

    else {
        prempo_t prempo;

        typename MPOTensor<MPSMatrix, SymmGroup>::op_table_ptr op_table(new OPTable<MPSMatrix, SymmGroup>);

        index_type b1, b2, b3;
        for (b1=0; b1 < mpo1.row_dim(); ++b1) {

            op_map out_row;

            row_proxy row1 = mpo1.row(b1);
            for (typename row_proxy::const_iterator it1 = row1.begin(); it1 != row1.end(); ++it1) {

                b2 = it1.index();
                row_proxy row2 = mpo2.row(b2);
                for (typename row_proxy::const_iterator it2 = row2.begin(); it2 != row2.end(); ++it2)
                {
                    b3 = it2.index();
                    assert((mpo1.has(b1, b2) && mpo2.has(b2, b3)));
                    block_matrix<MPSMatrix, SymmGroup> product;

                    const_term_descriptor<MPSMatrix, SymmGroup> p1 = mpo1.at(b1,b2), p2 = mpo2.at(b2,b3);

                    op_kron(phys_i, p1.op, p2.op, product);
                    product *= (p1.scale * p2.scale);
                    out_row[b3] += product;
                }
            }

            for (typename op_map::iterator it = out_row.begin(); it != out_row.end(); ++it) {
                b3 = it->first;
                tag_type new_tag = op_table->register_op(it->second);
                prempo.push_back(boost::make_tuple(b1, b3, new_tag, 1.0));
            }
        } 

        using boost::tuples::get;
        MPOTensor<MPSMatrix, SymmGroup> mpo_big(mpo1.row_dim(), mpo2.col_dim(), prempo, op_table);

        return mpo_big;
    }
}

template<class MPOMatrix, class MPSMatrix, class SymmGroup>
void make_ts_cache_mpo(MPO<MPOMatrix, SymmGroup> const & mpo_orig,
                       MPO<MPSMatrix, SymmGroup> & mpo_out, Index<SymmGroup> const & site_dim)
{
    std::size_t L_ts = mpo_orig.length() - 1;
    mpo_out.resize(L_ts);

    bool global_table = true;    
    for (int p=0; p<L_ts && global_table; ++p)
        global_table = (mpo_orig[p].get_operator_table() == mpo_orig[0].get_operator_table());

    // For now until above function is parallel
    parallel_for(locale::compact(L_ts), locale p = 0; p < L_ts; ++p)
        mpo_out[p] = make_twosite_mpo<MPOMatrix, MPSMatrix>(mpo_orig[p], mpo_orig[p+1], site_dim, global_table);
        
    std::size_t ntags=0;
    for (int p=0; p<mpo_out.length(); ++p) {
        ntags += mpo_out[p].get_operator_table()->size();
    }
    maquis::cout << "Total number of tags: " << ntags << std::endl;

}

#endif
