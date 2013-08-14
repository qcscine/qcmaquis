/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2013-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MAQUIS_DMRG_MODELS_OP_HANDLER_H
#define MAQUIS_DMRG_MODELS_OP_HANDLER_H

#include <vector>
#include <map>
#include <utility>
#include <stdexcept>

#include <boost/shared_ptr.hpp>

#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"

#include "dmrg/models/tag_detail.h"

template <class Matrix, class SymmGroup>
class OPTable : public std::vector<block_matrix<Matrix, SymmGroup> >
{
public:
    typedef tag_detail::tag_type tag_type;
    typedef block_matrix<Matrix, SymmGroup> op_t;

private:
    typedef typename Matrix::value_type mvalue_type;

public:
    tag_type register_op(op_t const & op_);
    std::pair<tag_type, mvalue_type> checked_register(op_t & sample);
};

template <class Matrix, class SymmGroup>
class TagHandler
{
public:
    typedef typename OPTable<Matrix, SymmGroup>::tag_type tag_type;
    typedef typename OPTable<Matrix, SymmGroup>::op_t op_t;

protected:
    typedef typename Matrix::value_type value_type;
    typedef std::map<std::pair<tag_type, tag_type>, std::pair<tag_type, value_type>, tag_detail::pair_cmp> pair_map_t;
    typedef typename pair_map_t::const_iterator pair_map_it_t;

public:
    TagHandler() :
        operator_table(new OPTable<Matrix, SymmGroup>())
        { }

    TagHandler(boost::shared_ptr<OPTable<Matrix, SymmGroup> > tbl_) :
        operator_table(tbl_)
        { }
    
    tag_type register_op(const op_t & op_, tag_detail::operator_kind kind);

    std::pair<tag_type, value_type> checked_register(op_t & sample, tag_detail::operator_kind kind) {
        std::pair<tag_type, value_type> ret = operator_table->checked_register(sample);
        if (sign_table.size() < operator_table->size())
            sign_table.push_back(kind);

        assert(sign_table.size() == operator_table->size());

        return ret;
    }

/*
    std::pair<tag_type, value_type> checked_register(op_t & sample, tag_detail::operator_kind kind) {
        std::pair<bool, value_type> cmp_result;
        typename std::vector<op_t>::iterator it_pt = operator_table->begin();
        for (; it_pt != operator_table->end(); ++it_pt) {
            cmp_result = equal(*it_pt, sample);
            if (cmp_result.first)
                break;
        }

        std::pair<tag_type, value_type> ret;
        if (it_pt == operator_table->end()) {
            ret = std::make_pair(operator_table->register_op(sample), 1.0);
            sign_table.push_back(kind);
        } else
            ret = std::make_pair(it_pt - operator_table->begin(), cmp_result.second);

        assert(sign_table.size() == operator_table->size());
        return ret;
    }
*/

    typename OPTable<Matrix, SymmGroup>::value_type & get_op(tag_type i) { return (*operator_table)[i]; }
    typename OPTable<Matrix, SymmGroup>::value_type const & get_op(tag_type i) const { return (*operator_table)[i]; }

    bool is_fermionic(tag_type query_tag) { return sign_table[query_tag]; }

    /* WARNING: not thread safe! */
    std::pair<tag_type, value_type> get_product_tag(const tag_type t1, const tag_type t2);


    /* Diagnostics *************************************/
    tag_type prod_duplicates() const { return duplicates_(product_tags); }

    tag_type get_num_products() const;
    tag_type total_size() const { return operator_table->size(); }
    /***************************************************/

    boost::shared_ptr<OPTable<Matrix, SymmGroup> > get_operator_table() { return operator_table; }

private:
    boost::shared_ptr<OPTable<Matrix, SymmGroup> > operator_table;     

    template <class Map>
    tag_type duplicates_(Map const & sample);

    std::vector<tag_detail::operator_kind> sign_table;
    pair_map_t product_tags;
};

template <class Matrix, class SymmGroup>
class KronHandler : public TagHandler<Matrix, SymmGroup>
{
    typedef TagHandler<Matrix, SymmGroup> base;
    typedef typename OPTable<Matrix, SymmGroup>::tag_type tag_type;
    using typename base::op_t;

public:

    KronHandler(boost::shared_ptr<OPTable<Matrix, SymmGroup> > tbl_)
    :   base(tbl_)
      , kronecker_table(new OPTable<Matrix, SymmGroup>())
    {
        for (typename OPTable<Matrix, SymmGroup>::iterator it = tbl_->begin();
             it != tbl_->end(); ++it)
            uniform.push_back(tag_detail::is_uniform(*it));
    }

    tag_type get_kron_tag(Index<SymmGroup> const & phys_i, tag_type t1, tag_type t2);

    typename OPTable<Matrix, SymmGroup>::value_type & get_op(tag_type i) { return (*kronecker_table)[i]; }
    typename OPTable<Matrix, SymmGroup>::value_type const & get_op(tag_type i) const { return (*kronecker_table)[i]; }

    bool is_uniform(tag_type t) { return uniform[t]; }

    boost::shared_ptr<OPTable<Matrix, SymmGroup> > get_kronecker_table() { return kronecker_table; }

    /* Diagnostics *************************************/
    tag_type get_num_kron_products() const;

private:
    boost::shared_ptr<OPTable<Matrix, SymmGroup> > kronecker_table;
    typename base::pair_map_t kron_tags;
    std::vector<bool> uniform;
};


#include "dmrg/models/op_handler.hpp"

#endif
