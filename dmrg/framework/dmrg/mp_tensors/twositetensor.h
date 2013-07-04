/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef TWOSITETENSOR_H
#define TWOSITETENSOR_H

#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/mpotensor.h"
#include "dmrg/block_matrix/indexing.h"
#include "dmrg/block_matrix/multi_index.h"
#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"

#include <iostream>
#include <algorithm>

enum TwoSiteStorageLayout {TSRightPaired, TSLeftPaired, TSBothPaired};

template<class Matrix, class SymmGroup>
class TwoSiteTensor
{
public:
    typedef std::size_t size_type;
    typedef typename MultiIndex<SymmGroup>::index_id index_id;
    typedef typename MultiIndex<SymmGroup>::set_id set_id;
    
    TwoSiteTensor(MPSTensor<Matrix, SymmGroup> const & mps1,
                  MPSTensor<Matrix, SymmGroup> const & mps2);

    TwoSiteTensor(MPSTensor<Matrix, SymmGroup> const & twin_mps);

    Index<SymmGroup> const & site_dim() const;
    Index<SymmGroup> const & row_dim() const;
    Index<SymmGroup> const & col_dim() const;
    
    block_matrix<Matrix, SymmGroup> & data();
    block_matrix<Matrix, SymmGroup> const & data() const;
    
    template<class Matrix_, class SymmGroup_>
    friend std::ostream& operator<<(std::ostream&, TwoSiteTensor<Matrix_, SymmGroup_> const &);

    TwoSiteTensor<Matrix, SymmGroup> & operator << ( MPSTensor<Matrix, SymmGroup> const & rhs);
    
    friend struct contraction;
    friend struct compression;
    friend struct multigrid;

    void make_left_paired() const;
    void make_both_paired() const;
    void make_right_paired() const;
    
    MPSTensor<Matrix, SymmGroup> make_mps() const;
    std::pair<MPSTensor<Matrix, SymmGroup>,
              MPSTensor<Matrix, SymmGroup> > split_mps_l2r(std::size_t Mmax, double cutoff, Logger * iter_log = NULL) const;
    std::pair<MPSTensor<Matrix, SymmGroup>,
              MPSTensor<Matrix, SymmGroup> > split_mps_r2l(std::size_t Mmax, double cutoff, Logger * iter_log = NULL) const;
    
    void swap_with(TwoSiteTensor & b);

    friend void swap(TwoSiteTensor & a, TwoSiteTensor & b)
    {
        a.swap_with(b);
    }
    
#ifdef HAVE_ALPS_HDF5
    void load(alps::hdf5::archive & ar);
    void save(alps::hdf5::archive & ar) const;
#endif
    
private:
    MultiIndex<SymmGroup> midx;
    set_id left_paired;
    set_id right_paired;
    set_id both_paired;

    Index<SymmGroup> phys_i, phys_i_orig, left_i, right_i;
    mutable block_matrix<Matrix, SymmGroup> data_;
    mutable TwoSiteStorageLayout cur_storage;
    Indicator cur_normalization;
};

#include "twositetensor.hpp"

#include "dmrg/mp_tensors/mpotensor.h"
template<class MPOMatrix, class MPSMatrix, class SymmGroup>
MPOTensor<MPSMatrix, SymmGroup> make_twosite_mpo(MPOTensor<MPOMatrix, SymmGroup> const & mpo1,
                                                 MPOTensor<MPOMatrix, SymmGroup> const & mpo2, Index<SymmGroup> const & phys_i)
{
    assert(mpo1.col_dim() == mpo2.row_dim());

    // Fast traversal if tags are available
    if (mpo1.tag_ready() && mpo2.tag_ready()) {

        typedef typename MPOTensor<MPOMatrix, SymmGroup>::index_type index_type;
        typedef typename MPOTensor<MPOMatrix, SymmGroup>::row_proxy row_proxy;
        typedef typename MPOTensor<MPOMatrix, SymmGroup>::col_proxy col_proxy;
        typedef typename OPTagTable<MPSMatrix, SymmGroup>::op_tag_t op_tag_t;
        typedef typename OPTagTable<MPSMatrix, SymmGroup>::op_t op_t;
        typedef typename MPSMatrix::value_type value_type;
        typedef std::map<index_type, op_t> op_map;
        typedef std::map<index_type, std::pair<op_tag_t, value_type> > op_scale_map;

        typedef std::vector<boost::tuple<index_type, index_type, op_tag_t, value_type> > prempo_t;
        prempo_t prempo;

        // TODO: This is ugly, improve interface... (maybe make this a MPO member or friend function)
        typename MPOTensor<MPSMatrix, SymmGroup>::tag_table_ptr op_table = mpo1.get_tag_table();

        // use global tag table for all mpotensors
        OPTagTable<MPSMatrix, SymmGroup> & kron_table = *op_table;

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
                    op_tag_t tmptag;

                    std::pair<typename OPTagTable<MPSMatrix, SymmGroup>::op_t const &, typename MPSMatrix::value_type>
                        p1 = mpo1.at(b1,b2), p2 = mpo2.at(b2,b3);

                    if (op_table->is_site_op(mpo1.tag_number(b1,b2)) || op_table->is_site_op(mpo2.tag_number(b2,b3))) {
                        non_uniform.insert(b3);
                    } else {
                        if (non_uniform.count(b3) == 0) {

                            #ifdef MAQUIS_OPENMP
                            #pragma omp critical
                            #endif
                            tmptag = kron_table.get_kron_tag(phys_i, mpo1.tag_number(b1,b2), mpo2.tag_number(b2,b3));

                            if (uniform_ops.count(b3) == 0) {
                                uniform_ops[b3].first = tmptag;
                                uniform_ops[b3].second = (p1.second * p2.second);

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

                    #ifdef MAQUIS_OPENMP
                    #pragma omp critical
                    #endif
                    tmp = kron_table[kron_table.get_kron_tag(phys_i, mpo1.tag_number(b1,b2), mpo2.tag_number(b2,b3))];
                    tmp *= (p1.second * p2.second);
                    kron_sums[b3] += tmp;
                }

            }

            for (typename op_scale_map::iterator it = uniform_ops.begin(); it != uniform_ops.end(); ++it) {
                b3 = it->first;
                //mpo_big(b1, b3) = it->second.second * kron_table[it->second.first];

                prempo.push_back(boost::make_tuple(b1, b3, it->second.first, it->second.second));
            }

            for (typename std::set<index_type>::iterator it = non_uniform.begin(); it != non_uniform.end(); ++it) {
                b3 = *it;
                op_t & tmp = kron_sums[b3];
                //mpo_big(b1, b3) = tmp;

                std::pair<op_tag_t, value_type> scaled_tag = kron_table.checked_register(tmp);
                prempo.push_back(boost::make_tuple(b1, b3, scaled_tag.first, scaled_tag.second));
            }
        } 

        using boost::tuples::get;
        MPOTensor<MPSMatrix, SymmGroup> mpo_big_tag(mpo1.row_dim(), mpo2.col_dim(), prempo, op_table);

        /* Traditional block_matrix filling */
        /*
        for (std::size_t i = 0; i < prempo.size(); ++i){
            index_type i1 = get<0>(prempo[i]);
            index_type i2 = get<1>(prempo[i]);
            op_tag_t tag = get<2>(prempo[i]);
            value_type scale = get<3>(prempo[i]);
            mpo_big_tag(i1, i2) = scale * kron_table[tag];
        }
        */
        return mpo_big_tag;

    }

    else {
        maquis::cout << "MPOTensor not tag ready\n";

        MPOTensor<MPSMatrix, SymmGroup> mpo_big(mpo1.row_dim(), mpo2.col_dim());
        // TODO: use OpenMP, thread-safe reduction needed!
        std::size_t b1, b2, b3;
        for( b1=0; b1 < mpo1.row_dim(); ++b1)
            for( b2=0; b2 < mpo1.col_dim(); ++b2)
                for( b3=0; b3 < mpo2.col_dim(); ++b3)
                {
                    if (! (mpo1.has(b1, b2) && mpo2.has(b2, b3)) )
                        continue;
                    block_matrix<MPSMatrix, SymmGroup> tmp;
                    op_kron(phys_i, mpo1(b1,b2), mpo2(b2,b3), tmp);
                    if (mpo_big.has(b1,b3))
                        mpo_big(b1,b3) += tmp;
                    else
                        mpo_big(b1,b3) = tmp;
                }

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
    boost::shared_ptr<OPTagTable<MPSMatrix, SymmGroup> > op_table = mpo_orig[0].get_tag_table();
    maquis::cout << "op_table total, site, prod, krons: " << (*op_table).total_size() << ", " << (*op_table).get_num_site_terms() << ", " 
                 << (*op_table).get_num_products() << ", " << (*op_table).get_num_kron_products() << std::endl;
    maquis::cout << "duplicates in prod/kron_table: " << (*op_table).prod_duplicates() << ", " << (*op_table).kron_duplicates() << std::endl;

}

#endif
