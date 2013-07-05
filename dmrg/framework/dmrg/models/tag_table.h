/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2013-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MAQUIS_DMRG_MODELS_TAG_TABLE_H
#define MAQUIS_DMRG_MODELS_TAG_TABLE_H

#include <boost/shared_ptr.hpp>

#include <vector>
#include <map>
#include <utility>
#include <stdexcept>

#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"

namespace tag_detail {

    typedef unsigned op_tag_t;

    struct pair_cmp
    {
        bool operator()(std::pair<op_tag_t, op_tag_t> const & i,
                        std::pair<op_tag_t, op_tag_t> const & j) const
        {
            if (i.first < j.first)
                return true;
            else if (i.first > j.first)
                return false;
            else
                return i.second < j.second;
        }
    };
}

template <class Matrix, class SymmGroup>
class OPTagTable : public std::vector<block_matrix<Matrix, SymmGroup> >
{
public:
    typedef tag_detail::op_tag_t op_tag_t;
    typedef block_matrix<Matrix, SymmGroup> op_t;
private:
    typedef typename Matrix::value_type value_type;
    typedef std::map<std::pair<op_tag_t, op_tag_t>, std::pair<op_tag_t, value_type>, tag_detail::pair_cmp> pair_map_t;
    typedef typename pair_map_t::const_iterator pair_map_it_t;
public:

    op_tag_t
    register_op(const op_t & op_);

    std::pair<op_tag_t, value_type>
    checked_register(op_t & sample);

    op_tag_t
    register_site_op(const op_t & op_);

    /* WARNING: not thread safe! */
    std::pair<op_tag_t, value_type>
    get_product_tag(const op_tag_t t1, const op_tag_t t2);

    op_tag_t
    get_kron_tag(Index<SymmGroup> const & phys_i, const op_tag_t t1, const op_tag_t t2);

    // TODO: fix const_element_iterator bug in alps::numeric::matrix to restore const-correctness here
    /* Check if two operators are equal modulo a scale factor*/
    static std::pair<bool, value_type>
    equal(op_t & reference, op_t & sample);

    //static bool full_equal(op_t & op1, op_t & op2);

    /* Diagnostics *************************************/
    op_tag_t kron_duplicates() const { return duplicates_(kron_tags); }
    op_tag_t prod_duplicates() const { return duplicates_(product_tags); }

    op_tag_t get_num_products() const;
    op_tag_t get_num_kron_products() const;
    op_tag_t get_num_site_terms() const { return site_terms.size(); }
    op_tag_t total_size() const { return this->size(); }

    /***************************************************/

    bool is_site_op(op_tag_t tag_) const { return site_terms.count(tag_) > 0; }

private:

    // slow N^2 algorithm, use hashes to get NlogN if necessary    
    template <class Map>
    op_tag_t duplicates_(Map & sample);

    pair_map_t product_tags;
    pair_map_t kron_tags;
    // Keep track of site_terms, they are non-uniformly scaled
    std::set<op_tag_t> site_terms;
};

#include "dmrg/models/tag_table.hpp"

#endif
