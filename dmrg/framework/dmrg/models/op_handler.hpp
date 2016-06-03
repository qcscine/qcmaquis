/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *
 * 
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 * 
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#ifndef MAQUIS_DMRG_MODELS_OP_HANDLER_HPP
#define MAQUIS_DMRG_MODELS_OP_HANDLER_HPP

// **************************************************************************
// OPTable implementation
// **************************************************************************

template <class Matrix, class SymmGroup>
TagHandler<Matrix, SymmGroup>::TagHandler(TagHandler const & rhs)
    : operator_table(new OPTable<Matrix, SymmGroup>(*rhs.operator_table))
    , sign_table(rhs.sign_table)
    , product_tags(rhs.product_tags)
{
}

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
OPTable<Matrix, SymmGroup>::checked_register(op_t const& sample)
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

// **************************************************************************
// TagHandler implementation
// **************************************************************************

// simple const query
template <class Matrix, class SymmGroup>
typename OPTable<Matrix, SymmGroup>::tag_type TagHandler<Matrix, SymmGroup>::size() const
{
    return operator_table->size();
}

template <class Matrix, class SymmGroup>
boost::shared_ptr<OPTable<Matrix, SymmGroup> > TagHandler<Matrix, SymmGroup>::get_operator_table() const
{
    return operator_table;
}

template <class Matrix, class SymmGroup>
bool TagHandler<Matrix, SymmGroup>::is_fermionic(typename OPTable<Matrix, SymmGroup>::tag_type query_tag) const
{
    assert(query_tag < sign_table.size());
    return sign_table[query_tag];
}

template <class Matrix, class SymmGroup>
typename OPTable<Matrix, SymmGroup>::tag_type TagHandler<Matrix, SymmGroup>::
herm_conj(typename OPTable<Matrix, SymmGroup>::tag_type query_tag) const
{
    assert(query_tag < hermitian.size());
    return hermitian[query_tag];
}

// register new operators
template <class Matrix, class SymmGroup>
typename OPTable<Matrix, SymmGroup>::tag_type TagHandler<Matrix, SymmGroup>::
register_op(const op_t & op_, tag_detail::operator_kind kind)
{
    sign_table.push_back(kind);
    tag_type ret = operator_table->register_op(op_);
    hermitian.push_back(ret);
    assert(sign_table.size() == operator_table->size());
    assert(hermitian.size() == operator_table->size());
    assert(ret < operator_table->size());
    return ret;
}

template <class Matrix, class SymmGroup>
std::pair<typename OPTable<Matrix, SymmGroup>::tag_type,
          typename TagHandler<Matrix,SymmGroup>::value_type> TagHandler<Matrix, SymmGroup>::
checked_register(typename OPTable<Matrix, SymmGroup>::op_t const& sample, tag_detail::operator_kind kind)
{
    std::pair<tag_type, value_type> ret = operator_table->checked_register(sample);
    if (sign_table.size() < operator_table->size())
    {
        sign_table.push_back(kind);
        hermitian.push_back(ret.first);
    }
    
    assert(sign_table.size() == operator_table->size());
    assert(hermitian.size() == operator_table->size());
    assert(ret.first < operator_table->size());
    
    return ret;
}

template <class Matrix, class SymmGroup>
void TagHandler<Matrix, SymmGroup>::hermitian_pair(typename OPTable<Matrix, SymmGroup>::tag_type pair_tag1,
                                                   typename OPTable<Matrix, SymmGroup>::tag_type pair_tag2)
{
    assert(std::max(pair_tag1, pair_tag2) < hermitian.size());
    assert(pair_tag1 != pair_tag2);
    
    if (hermitian[pair_tag1] == pair_tag2 && hermitian[pair_tag2] == pair_tag1) return;
    assert(hermitian[pair_tag1] == pair_tag1 && hermitian[pair_tag2] == pair_tag2);
    std::swap(hermitian[pair_tag1], hermitian[pair_tag2]);
}

// access operators
template <class Matrix, class SymmGroup>
typename OPTable<Matrix, SymmGroup>::value_type & TagHandler<Matrix, SymmGroup>::get_op(tag_type i) { return (*operator_table)[i]; }

template <class Matrix, class SymmGroup>
typename OPTable<Matrix, SymmGroup>::value_type const & TagHandler<Matrix, SymmGroup>::get_op(tag_type i) const { return (*operator_table)[i]; }

template <class Matrix, class SymmGroup>
std::vector<typename OPTable<Matrix, SymmGroup>::value_type> TagHandler<Matrix, SymmGroup>::get_ops(std::vector<tag_type> const & tags) const
{
    std::vector<typename OPTable<Matrix, SymmGroup>::value_type> ret(tags.size());
    for (int k = 0; k < tags.size(); ++k)
        ret[k] = (*operator_table)[tags[k]];

    return ret;
}

// compute products
template <class Matrix, class SymmGroup>
std::pair<typename OPTable<Matrix, SymmGroup>::tag_type,
          typename TagHandler<Matrix,SymmGroup>::value_type> TagHandler<Matrix, SymmGroup>::
get_product_tag(const typename OPTable<Matrix, SymmGroup>::tag_type t1,
                const typename OPTable<Matrix, SymmGroup>::tag_type t2)
{
    assert( t1 < operator_table->size() && t2 < operator_table->size() );
    assert( operator_table->size() == sign_table.size());

    // return tag of product, if already there
    try {
#if defined(__xlC__) || defined(__FCC_VERSION)
        if (product_tags.count(std::make_pair(t1, t2)) == 0)
            throw std::out_of_range("");

        return product_tags[std::make_pair(t1, t2)];
#else
        return product_tags.at(std::make_pair(t1, t2));
#endif
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

        // set the product spin descriptor
        product.spin() = couple(get_op(t2).spin(), get_op(t1).spin());

        std::pair<tag_type, value_type> ret = this->checked_register(product, prod_kind);
        product_tags[std::make_pair(t1, t2)] = ret;
        assert( operator_table->size() == sign_table.size());
        assert( ret.first < operator_table->size() );
        return ret;
    }
}

template <class Matrix, class SymmGroup>
std::pair<std::vector<typename OPTable<Matrix, SymmGroup>::tag_type>, std::vector<typename TagHandler<Matrix,SymmGroup>::value_type> >
TagHandler<Matrix, SymmGroup>::
get_product_tags(std::vector<typename OPTable<Matrix, SymmGroup>::tag_type> const & ops1,
                 std::vector<typename OPTable<Matrix, SymmGroup>::tag_type> const & ops2)
{
    assert(ops1.size() == ops2.size());
    std::pair<std::vector<tag_type>, std::vector<value_type> >ret;
    for (typename SymmGroup::subcharge sc=0; sc < ops1.size(); ++sc) {
        std::pair<tag_type, value_type> ptag = this->get_product_tag(ops1[sc], ops2[sc]);
        ret.first.push_back(ptag.first);
        ret.second.push_back(ptag.second);
    }

    return ret;
}

// * Diagnostics *************************************
template <class Matrix, class SymmGroup>
template <class Map>
typename OPTable<Matrix, SymmGroup>::tag_type TagHandler<Matrix, SymmGroup>::
duplicates_(Map const & sample)
{
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


// **************************************************************************
// KronHandler implementation
// **************************************************************************

template <class Matrix, class SymmGroup>
typename OPTable<Matrix, SymmGroup>::tag_type KronHandler<Matrix, SymmGroup>::
get_kron_tag(Index<SymmGroup> const & phys_i1,
             Index<SymmGroup> const & phys_i2,
             typename OPTable<Matrix, SymmGroup>::tag_type t1,
             typename OPTable<Matrix, SymmGroup>::tag_type t2,
             SpinDescriptor<typename symm_traits::SymmType<SymmGroup>::type> lspin,
             SpinDescriptor<typename symm_traits::SymmType<SymmGroup>::type> mspin,
             SpinDescriptor<typename symm_traits::SymmType<SymmGroup>::type> rspin)
{
    assert( t1 < base::get_operator_table()->size() && t2 < base::get_operator_table()->size() );

    // return tag of kronecker product, if already there
    try {
#if defined(__xlC__) || defined(__FCC_VERSION)
        if (kron_tags.count(std::make_pair(t1, t2)) == 0)
            throw std::out_of_range("");

        return kron_tags[std::make_pair(t1, t2)].first;
#else
        return kron_tags.at(std::make_pair(t1, t2)).first;
#endif
    }
    // compute and register the product, then return the new tag
    catch(const std::out_of_range& e) {

        op_t product;
        op_t const& op1 = (*base::get_operator_table())[t1];
        op_t const& op2 = (*base::get_operator_table())[t2];

        op_kron(phys_i1, phys_i2, op1, op2, product, lspin, mspin, rspin);

        tag_detail::remove_empty_blocks(product);
        
        tag_type ret = kronecker_table->register_op(product);
        kron_tags[std::make_pair(t1, t2)] = std::make_pair(ret, 1.0);

        return ret;
    }
}

template <class Matrix, class SymmGroup>
typename OPTable<Matrix, SymmGroup>::tag_type KronHandler<Matrix, SymmGroup>::get_num_kron_products() const
{
    return kronecker_table->size();
}

#endif
