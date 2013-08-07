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
                                                 MPOTensor<MPOMatrix, SymmGroup> const & mpo2, Index<SymmGroup> const & phys_i)
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

    if (mpo1.get_operator_table().get() == mpo2.get_operator_table().get()) {
        typedef std::map<index_type, std::pair<tag_type, value_type> > op_scale_map;

        prempo_t prempo;

        // TODO: This is ugly, improve interface... (maybe make this a MPO member or friend function)
        typename MPOTensor<MPSMatrix, SymmGroup>::op_table_ptr op_table = mpo1.get_operator_table();

        OPTable<MPSMatrix, SymmGroup> & kron_table = *op_table;

        /* new CSR code */
        index_type b1, b2, b3;
        for (b1=0; b1 < mpo1.row_dim(); ++b1) {

            op_map kron_sums;
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
                    block_matrix<MPSMatrix, SymmGroup> tmp;
                    tag_type tmptag;

                    const_term_descriptor<MPSMatrix, SymmGroup> p1 = mpo1.at(b1,b2), p2 = mpo2.at(b2,b3);

                    if (op_table->is_site_op(mpo1.tag_number(b1,b2)) || op_table->is_site_op(mpo2.tag_number(b2,b3))) {
                        non_uniform.insert(b3);
                    } else {
                        if (non_uniform.count(b3) == 0) {

                            // only needed if the operator_table is shared among mpotensors
                            #ifdef MAQUIS_OPENMP
                            #pragma omp critical
                            #endif
                            tmptag = kron_table.get_kron_tag(phys_i, mpo1.tag_number(b1,b2), mpo2.tag_number(b2,b3));

                            if (uniform_ops.count(b3) == 0) {
                                uniform_ops[b3].first = tmptag;
                                uniform_ops[b3].second = (p1.scale * p2.scale);

                            } else { /* Have to add kronecker products, very rare */
                                if (kron_table[uniform_ops[b3].first].left_basis() != kron_table[tmptag].left_basis() ||
                                    kron_table[uniform_ops[b3].first].right_basis() != kron_table[tmptag].right_basis())

                                    non_uniform.insert(b3);
                                else
                                    throw std::runtime_error("Need to get tag for a kronecker sum, which is not implemented\n");
                            }
                        }

                        else
                            uniform_ops.erase(b3);
                    }

                    // only needed if the operator_table is shared among mpotensors
                    #ifdef MAQUIS_OPENMP
                    #pragma omp critical
                    #endif
                    tmp = kron_table[kron_table.get_kron_tag(phys_i, mpo1.tag_number(b1,b2), mpo2.tag_number(b2,b3))];
                    tmp *= (p1.scale * p2.scale);
                    kron_sums[b3] += tmp;
                }

            }

            for (typename op_scale_map::iterator it = uniform_ops.begin(); it != uniform_ops.end(); ++it) {
                b3 = it->first;
                prempo.push_back(boost::make_tuple(b1, b3, it->second.first, it->second.second));
            }

            for (typename std::set<index_type>::iterator it = non_uniform.begin(); it != non_uniform.end(); ++it) {
                b3 = *it;
                op_t & tmp = kron_sums[b3];

                std::pair<tag_type, value_type> scaled_tag = kron_table.checked_register(tmp);
                prempo.push_back(boost::make_tuple(b1, b3, scaled_tag.first, scaled_tag.second));
            }
        } 

        using boost::tuples::get;
        MPOTensor<MPSMatrix, SymmGroup> mpo_big_tag(mpo1.row_dim(), mpo2.col_dim(), prempo, op_table);

        return mpo_big_tag;

    }

    else {
        prempo_t prempo;

        typename MPOTensor<MPSMatrix, SymmGroup>::op_table_ptr op_table(new OPTable<MPSMatrix, SymmGroup>);
        OPTable<MPSMatrix, SymmGroup> & kron_table = *op_table;

        index_type b1, b2, b3;
        for (b1=0; b1 < mpo1.row_dim(); ++b1) {

            op_map kron_sums;

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
                    kron_sums[b3] += product;
                }
            }

            for (typename op_map::iterator it = kron_sums.begin(); it != kron_sums.end(); ++it) {
                b3 = it->first;
                tag_type new_tag = kron_table.register_op(it->second);
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
    // For now until above function is parallel
    parallel_for(locale::compact(L_ts), locale p = 0; p < L_ts; ++p)
        mpo_out[p] = make_twosite_mpo<MPOMatrix, MPSMatrix>(mpo_orig[p], mpo_orig[p+1], site_dim);
        
    /* Diagnosis */
    boost::shared_ptr<OPTable<MPSMatrix, SymmGroup> > op_table = mpo_orig[0].get_operator_table();
    maquis::cout << "op_table total, site, prod, krons: " << (*op_table).total_size() << ", " << (*op_table).get_num_site_terms() << ", " 
                 << (*op_table).get_num_products() << ", " << (*op_table).get_num_kron_products() << std::endl;
    //maquis::cout << "duplicates in prod/kron_table: " << (*op_table).prod_duplicates() << ", " << (*op_table).kron_duplicates() << std::endl;

}

#endif
