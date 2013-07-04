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
    typedef std::map<std::pair<op_tag_t, op_tag_t>, op_tag_t, tag_detail::pair_cmp> pair_map_t;
    typedef typename pair_map_t::const_iterator pair_map_it_t;
    typedef typename Matrix::value_type value_type;
public:

    op_tag_t register_op(const op_t & op_) {
        op_tag_t ret = this->size();
        this->push_back(op_);
        return ret;
    }

    std::pair<op_tag_t, value_type> checked_register(op_t & sample) {

        std::pair<bool, value_type> cmp_result;
        typename std::vector<op_t>::iterator it_pt = this->begin();
        for (; it_pt != this->end(); ++it_pt) { 
            cmp_result = equal(*it_pt, sample);
            if (cmp_result.first)
                break;
        }

        if (it_pt == this->end()) {
            return std::make_pair(this->register_op(sample), 1.0);
            maquis::cout << "checked_register: had to add new operator\n";
        } else
            return std::make_pair(it_pt - this->begin(), cmp_result.second);
        
    }

    op_tag_t register_site_op(const op_t & op_) {
        op_tag_t ret = this->register_op(op_);
        site_terms.insert(ret);
        return ret;
    }

    /* WARNING: not thread safe! */
    op_tag_t get_prod_tag(const op_tag_t t1, const op_tag_t t2) {
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

            // Check table if we have the product already
            typename std::vector<op_t>::iterator it_pt = this->begin();
            for (; it_pt != this->end(); ++it_pt) 
                if (full_equal(*it_pt, product))
                    break;

            op_tag_t ret;
            if (it_pt == this->end()) {
                ret = this->register_op(product);
                product_tags[std::make_pair(t1, t2)] = ret;
            } else {
                ret = it_pt - this->begin();
                product_tags[std::make_pair(t1, t2)] = ret;
            }

            return ret;
        }
    }

    /* WARNING: not thread safe! */
    op_tag_t get_kron_tag(Index<SymmGroup> const & phys_i, const op_tag_t t1, const op_tag_t t2) {
        assert( t1 < this->size() && t2 < this->size() );

        // return tag of kronecker product, if already there
        try {
            return kron_tags.at(std::make_pair(t1, t2));
        }
        // compute and register the product, then return the new tag
        catch(const std::out_of_range& e) {

            assert(t1 < this->size() && t2 < this->size());

            op_t product;
            op_t& op1 = (*this)[t1];
            op_t& op2 = (*this)[t2];

            op_kron(phys_i, op1, op2, product);
            op_tag_t ret = this->register_op(product);
            kron_tags[std::make_pair(t1, t2)] = ret;

            return ret;
        }
    }

    /* Check if two operators are equal modulo a scale factor*/
    // TODO: fix const_element_iterator bug in alps::numeric::matrix to restore const-correctness here
    static std::pair<bool, value_type> equal(op_t & reference, op_t & sample) {
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

    static bool full_equal(op_t & op1, op_t & op2) {
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

    /* Diagnostics *************************************/

    // slow N^2 algorithm, use hashes to get NlogN if necessary    
    template <class Map>
    op_tag_t duplicates_(Map & sample) {
        typedef typename Map::const_iterator it_t;

        std::vector<op_tag_t> unique_ops;
        for (it_t it_s = sample.begin(); it_s != sample.end(); ++it_s) 
        {
            bool unique = true;
            for (typename std::vector<op_tag_t>::iterator it_unique = unique_ops.begin(); it_unique != unique_ops.end(); ++it_unique)
                if (equal((*this)[(*it_s).second], (*this)[*it_unique]).first)
                {
                    unique = false;
                    break;
                }

            if (unique)
                unique_ops.push_back((*it_s).second);
        }

        return sample.size() - unique_ops.size();
    }

    op_tag_t kron_duplicates() { return duplicates_(kron_tags); }
    op_tag_t prod_duplicates() { return duplicates_(product_tags); }


    op_tag_t get_num_products() const {
        std::set<op_tag_t> utags;
        for (pair_map_it_t it = product_tags.begin(); it != product_tags.end(); ++it)
            utags.insert(it->second);

        return utags.size();
    }
    op_tag_t get_num_kron_products() const {
        std::set<op_tag_t> utags;
        for (pair_map_it_t it = kron_tags.begin(); it != kron_tags.end(); ++it)
            utags.insert(it->second);

        return utags.size();
    }
    op_tag_t get_num_site_terms() const { return site_terms.size(); }
    op_tag_t total_size() const { return this->size(); }

    /***************************************************/

    bool is_site_op(op_tag_t tag_) const { return site_terms.count(tag_) > 0; }

private:

    pair_map_t product_tags;
    pair_map_t kron_tags;
    // Keep track of site_terms, they are non-uniformly scaled
    std::set<op_tag_t> site_terms;
};

#endif
