/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2013-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MAQUIS_DMRG_MODELS_TAG_TABLE_HPP
#define MAQUIS_DMRG_MODELS_TAG_TABLE_HPP

template <class Matrix, class SymmGroup>
typename OPTable<Matrix, SymmGroup>::tag_type OPTable<Matrix, SymmGroup>::register_op(const op_t & op_) {
    tag_type ret = this->size();
    this->push_back(op_);
    return ret;
}

template <class Matrix, class SymmGroup>
std::pair<typename OPTable<Matrix, SymmGroup>::tag_type, typename OPTable<Matrix, SymmGroup>::value_type> OPTable<Matrix, SymmGroup>::checked_register(op_t & sample) {

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

template <class Matrix, class SymmGroup>
typename OPTable<Matrix, SymmGroup>::tag_type OPTable<Matrix, SymmGroup>::register_site_op(const op_t & op_) {
    tag_type ret = this->register_op(op_);
    site_terms.insert(ret);
    return ret;
}

template <class Matrix, class SymmGroup>
std::pair<typename OPTable<Matrix, SymmGroup>::tag_type,
          typename OPTable<Matrix, SymmGroup>::value_type> OPTable<Matrix, SymmGroup>::get_product_tag
                    (const typename OPTable<Matrix, SymmGroup>::tag_type t1,
                     const typename OPTable<Matrix, SymmGroup>::tag_type t2) {
    assert( t1 < this->size() && t2 < this->size() );

    // return tag of product, if already there
    try {
        return product_tags.at(std::make_pair(t1, t2));
    }

    // compute and register the product, then return the new tag
    catch(const std::out_of_range& e) {

        op_t product;
        op_t& op1 = (*this)[t1];
        op_t& op2 = (*this)[t2];

        gemm(op1, op2, product);
        std::pair<tag_type, value_type> ret = this->checked_register(product);
        product_tags[std::make_pair(t1, t2)] = ret;
        return ret;
    }
}

template <class Matrix, class SymmGroup>
typename OPTable<Matrix, SymmGroup>::tag_type OPTable<Matrix, SymmGroup>::get_kron_tag
            (Index<SymmGroup> const & phys_i,
            const typename OPTable<Matrix, SymmGroup>::tag_type t1,
            const typename OPTable<Matrix, SymmGroup>::tag_type t2) {
    assert( t1 < this->size() && t2 < this->size() );

    // return tag of kronecker product, if already there
    try {
        return kron_tags.at(std::make_pair(t1, t2)).first;
    }
    // compute and register the product, then return the new tag
    catch(const std::out_of_range& e) {

        assert(t1 < this->size() && t2 < this->size());

        op_t product;
        op_t& op1 = (*this)[t1];
        op_t& op2 = (*this)[t2];

        op_kron(phys_i, op1, op2, product);
        tag_type ret = this->register_op(product);
        kron_tags[std::make_pair(t1, t2)] = std::make_pair(ret, 1.0);

        return ret;
    }
}

template <class Matrix, class SymmGroup>
std::pair<bool, typename OPTable<Matrix, SymmGroup>::value_type> OPTable<Matrix, SymmGroup>::equal(op_t & reference, op_t & sample) {
    if (reference.left_basis() != sample.left_basis() || reference.right_basis() != sample.right_basis())
        return std::make_pair(false, 0.);

    typename Matrix::value_type invscale1, invscale2;
    
    // determine scale of matrices
    for (typename Matrix::const_element_iterator it = reference[0].elements().first;
            it != reference[0].elements().second; ++it)
        if (std::abs(*it) > 1.e-20) {
            invscale1 = 1./(*it);
            break;
        }
    for (typename Matrix::const_element_iterator it = sample[0].elements().first;
            it != sample[0].elements().second; ++it)
        if (std::abs(*it) > 1.e-20) {
            invscale2 = 1./(*it);
            break;
        }

    // Check all blocks for equality modulo scale factor
    for (typename Matrix::size_type b=0; b < reference.n_blocks(); ++b)
    {
        typename Matrix::const_element_iterator it1 = reference[b].elements().first, 
                                                it2 = sample[b].elements().first;
        for ( ; it1 != reference[b].elements().second; ++it1, ++it2)
            if (std::abs(*it1 * invscale1 - *it2 * invscale2) > 1e-20)
                return std::make_pair(false, 0.);
    }

    return std::make_pair(true, invscale1/invscale2);
}

/*
template <class Matrix, class SymmGroup>
bool OPTable<Matrix, SymmGroup>::full_equal(op_t & op1, op_t & op2) {
    if (op1.left_basis() != op2.left_basis() || op1.right_basis() != op2.right_basis())
        return false;

    // Check all blocks for equality
    for (typename Matrix::size_type b=0; b < op1.n_blocks(); ++b)
    {
        if (!std::equal(op1[b].elements().first, op1[b].elements().second, op2[b].elements().first))
            return false;
    }

    return true;
}
*/


// * Diagnostics *************************************
template <class Matrix, class SymmGroup>
template <class Map>
typename OPTable<Matrix, SymmGroup>::tag_type OPTable<Matrix, SymmGroup>::duplicates_(Map & sample) {
    typedef typename Map::const_iterator it_t;

    std::vector<tag_type> unique_ops;
    for (it_t it_s = sample.begin(); it_s != sample.end(); ++it_s) 
    {
        bool unique = true;
        for (typename std::vector<tag_type>::iterator it_unique = unique_ops.begin(); it_unique != unique_ops.end(); ++it_unique)
            if (equal((*this)[(*it_s).second.first], (*this)[*it_unique]).first)
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
typename OPTable<Matrix, SymmGroup>::tag_type OPTable<Matrix, SymmGroup>::get_num_products() const {
    std::set<tag_type> utags;
    for (pair_map_it_t it = product_tags.begin(); it != product_tags.end(); ++it)
        utags.insert(it->second.first);

    return utags.size();
}

template <class Matrix, class SymmGroup>
typename OPTable<Matrix, SymmGroup>::tag_type OPTable<Matrix, SymmGroup>::get_num_kron_products() const {
    std::set<tag_type> utags;
    for (pair_map_it_t it = kron_tags.begin(); it != kron_tags.end(); ++it)
        utags.insert(it->second.first);

    return utags.size();
}

#endif
