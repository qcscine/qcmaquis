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

    typedef unsigned tag_type;

    struct pair_cmp
    {
        bool operator()(std::pair<tag_type, tag_type> const & i,
                        std::pair<tag_type, tag_type> const & j) const
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
class OPTable : public std::vector<block_matrix<Matrix, SymmGroup> >
{
public:
    typedef tag_detail::tag_type tag_type;
    typedef block_matrix<Matrix, SymmGroup> op_t;

private:
    typedef typename Matrix::value_type value_type;
    typedef std::map<std::pair<tag_type, tag_type>, std::pair<tag_type, value_type>, tag_detail::pair_cmp> pair_map_t;
    typedef typename pair_map_t::const_iterator pair_map_it_t;

public:
    tag_type register_op(const op_t & op_);

    std::pair<tag_type, value_type> checked_register(op_t & sample);

    tag_type register_site_op(const op_t & op_);

    /* WARNING: not thread safe! */
    std::pair<tag_type, value_type> get_product_tag(const tag_type t1, const tag_type t2);

    tag_type get_kron_tag(Index<SymmGroup> const & phys_i, const tag_type t1, const tag_type t2);


    /* Diagnostics *************************************/
    tag_type kron_duplicates() const { return duplicates_(kron_tags); }
    tag_type prod_duplicates() const { return duplicates_(product_tags); }

    tag_type get_num_products() const;
    tag_type get_num_kron_products() const;
    tag_type get_num_site_terms() const { return site_terms.size(); }
    tag_type total_size() const { return this->size(); }

    /***************************************************/

    bool is_site_op(tag_type tag_) const { return site_terms.count(tag_) > 0; }

private:

    // TODO: fix const_element_iterator bug in alps::numeric::matrix to restore const-correctness here
    /* Check if two operators are equal modulo a scale factor*/
    static std::pair<bool, value_type> equal(op_t & reference, op_t & sample);

    //static bool full_equal(op_t & op1, op_t & op2);

    // slow N^2 algorithm, use hashes to get NlogN if necessary    
    template <class Map>
    tag_type duplicates_(Map & sample);

    pair_map_t product_tags;
    pair_map_t kron_tags;
    // Keep track of site_terms, they are non-uniformly scaled
    std::set<tag_type> site_terms;
};

#include "dmrg/models/tag_table.hpp"

#endif
