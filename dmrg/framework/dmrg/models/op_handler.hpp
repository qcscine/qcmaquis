/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2013-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *
 *
 *****************************************************************************/

#ifndef MAQUIS_DMRG_MODELS_OP_HANDLER_HPP
#define MAQUIS_DMRG_MODELS_OP_HANDLER_HPP

template <class Matrix, class SymmGroup>
typename OPTable<Matrix, SymmGroup>::tag_type
OPTable<Matrix, SymmGroup>::register_op(op_t const & op_)
{
    tag_type ret = this->size();
    this->push_back(op_);
    return ret;
}

template <class Matrix, class SymmGroup>
std::pair<typename OPTable<Matrix, SymmGroup>::tag_type, typename OPTable<Matrix, SymmGroup>::mvalue_type>
OPTable<Matrix, SymmGroup>::checked_register(op_t & sample)
{
    std::pair<bool, mvalue_type> cmp_result;
    typename std::vector<op_t>::iterator it_pt = this->begin();
    for (; it_pt != this->end(); ++it_pt) { 
        cmp_result = tag_detail::equal(*it_pt, sample);
        if (cmp_result.first)
            break;
    }

    if (it_pt == this->end()) {
        return std::make_pair(this->register_op(sample), 1.0);
    } else
        return std::make_pair(it_pt - this->begin(), cmp_result.second);
    
}

template <class Matrix, class SymmGroup>
typename OPTable<Matrix, SymmGroup>::tag_type TagHandler<Matrix, SymmGroup>::register_op(const op_t & op_, tag_detail::operator_kind kind) {
    sign_table.push_back(kind);
    return operator_table->register_op(op_);
}

/*
template <class Matrix, class SymmGroup>
std::pair<typename OPTable<Matrix, SymmGroup>::tag_type, typename OPTable<Matrix, SymmGroup>::value_type> TagHandler<Matrix, SymmGroup>::checked_register(op_t & sample) {

    std::pair<bool, value_type> cmp_result;
    typename std::vector<op_t>::iterator it_pt = this->begin();
    for (; it_pt != this->end(); ++it_pt) { 
        cmp_result = equal(*it_pt, sample);
        if (cmp_result.first)
            break;
    }

    if (it_pt == this->end()) {
        return std::make_pair(this->register_op(sample), 1.0);
    } else
        return std::make_pair(it_pt - this->begin(), cmp_result.second);
    
}
*/

template <class Matrix, class SymmGroup>
std::pair<typename OPTable<Matrix, SymmGroup>::tag_type,
          typename TagHandler<Matrix,SymmGroup>::value_type>
TagHandler<Matrix, SymmGroup>::get_product_tag(const typename
                                                     OPTable<Matrix, SymmGroup>::tag_type t1,
                                                     const typename
                                                     OPTable<Matrix, SymmGroup>::tag_type t2)
{
    assert( t1 < operator_table->size() && t2 < operator_table->size() );

    // return tag of product, if already there
    try {
        return product_tags.at(std::make_pair(t1, t2));
    }

    // compute and register the product, then return the new tag
    catch(const std::out_of_range& e) {

        op_t product;
        op_t& op1 = (*operator_table)[t1];
        op_t& op2 = (*operator_table)[t2];

        gemm(op1, op2, product);
        tag_detail::operator_kind prod_kind = tag_detail::bosonic;
        if (sign_table[t1] != sign_table[t2])
            prod_kind = tag_detail::fermionic;

        std::pair<tag_type, value_type> ret = this->checked_register(product, prod_kind);
        product_tags[std::make_pair(t1, t2)] = ret;
        return ret;
    }
}

// * Diagnostics *************************************
template <class Matrix, class SymmGroup>
template <class Map>
typename OPTable<Matrix, SymmGroup>::tag_type TagHandler<Matrix, SymmGroup>::duplicates_(Map const & sample) {
    typedef typename Map::const_iterator it_t;

    std::vector<tag_type> unique_ops;
    for (it_t it_s = sample.begin(); it_s != sample.end(); ++it_s) 
    {
        bool unique = true;
        for (typename std::vector<tag_type>::iterator it_unique = unique_ops.begin(); it_unique != unique_ops.end(); ++it_unique)
            if (equal((*operator_table)[(*it_s).second.first], (*operator_table)[*it_unique]).first)
            {
                unique = false;
                break;
            }

        if (unique)
            unique_ops.push_back((*it_s).second.first);
    }

    return sample.size() - unique_ops.size();
}

template <class Matrix, class SymmGroup>
typename OPTable<Matrix, SymmGroup>::tag_type TagHandler<Matrix, SymmGroup>::get_num_products() const {
    std::set<tag_type> utags;
    for (pair_map_it_t it = product_tags.begin(); it != product_tags.end(); ++it)
        utags.insert(it->second.first);

    return utags.size();
}
// ***************************************************

template <class Matrix, class SymmGroup>
typename OPTable<Matrix, SymmGroup>::tag_type KronHandler<Matrix, SymmGroup>::get_kron_tag
            (Index<SymmGroup> const & phys_i,
            typename OPTable<Matrix, SymmGroup>::tag_type t1,
            typename OPTable<Matrix, SymmGroup>::tag_type t2)
{
    assert( t1 < base::get_operator_table()->size() && t2 < base::get_operator_table()->size() );

    // return tag of kronecker product, if already there
    try {
        return kron_tags.at(std::make_pair(t1, t2)).first;
    }
    // compute and register the product, then return the new tag
    catch(const std::out_of_range& e) {

        op_t product;
        op_t& op1 = (*base::get_operator_table())[t1];
        op_t& op2 = (*base::get_operator_table())[t2];

        op_kron(phys_i, op1, op2, product);
        
        tag_type ret = kronecker_table->register_op(product);
        kron_tags[std::make_pair(t1, t2)] = std::make_pair(ret, 1.0);

        return ret;
    }
}

template <class Matrix, class SymmGroup>
typename OPTable<Matrix, SymmGroup>::tag_type KronHandler<Matrix, SymmGroup>::get_num_kron_products() const {
    return kronecker_table->size();
}

#endif
