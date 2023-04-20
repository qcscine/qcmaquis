/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef TAG_HANDLER_H
#define TAG_HANDLER_H

#include <vector>
#include <map>
#include <utility>
#include <stdexcept>

#include <boost/shared_ptr.hpp>
#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"
#include "dmrg/block_matrix/site_operator.h"
#include "dmrg/block_matrix/site_operator_algorithms.h"
#include "dmrg/models/tag_detail.h"

template <class Matrix, class SymmGroup>
class TagHandler
{
public:
    typedef typename OPTable<Matrix, SymmGroup>::tag_type tag_type;
    typedef typename OPTable<Matrix, SymmGroup>::op_t op_t;

protected:
    typedef typename Matrix::value_type value_type;
    typedef std::pair<tag_type, tag_type> tag_pair_t;
    typedef std::map<tag_pair_t, std::pair<tag_type, value_type>, compare_pair<tag_pair_t> > pair_map_t;
    typedef typename pair_map_t::const_iterator pair_map_it_t;

public:
    // constructors
    TagHandler() : operator_table(new OPTable<Matrix, SymmGroup>()) { }
    TagHandler(std::shared_ptr<OPTable<Matrix, SymmGroup> > tbl_) : operator_table(tbl_) { }
    TagHandler(TagHandler const & a);

    // simple const query
    tag_type size() const;
    std::shared_ptr<OPTable<Matrix, SymmGroup> > get_operator_table() const;
    bool is_fermionic (tag_type query_tag) const;
    tag_type herm_conj(tag_type query_tag) const;

    // register new operators
    tag_type register_op(const op_t & op_, tag_detail::operator_kind kind);
    std::pair<tag_type, value_type> checked_register(op_t const& sample, tag_detail::operator_kind kind);
    bool hasRegistered(const op_t& sample);

    void hermitian_pair(tag_type pair_tag1, tag_type pair_tag2);

    // access operators
    typename OPTable<Matrix, SymmGroup>::value_type & get_op(tag_type i);
    typename OPTable<Matrix, SymmGroup>::value_type const & get_op(tag_type i) const;
    std::vector<typename OPTable<Matrix, SymmGroup>::value_type> get_ops(std::vector<tag_type> const & i) const;

    // compute products (WARNING: not thread safe!)
    std::pair<tag_type, value_type> get_product_tag(const tag_type t1, const tag_type t2);
    std::pair<std::vector<tag_type>, std::vector<value_type> > get_product_tags(const std::vector<tag_type> & t1, const std::vector<tag_type> & t2);

    // Diagnostics
    tag_type prod_duplicates() const { return duplicates_(product_tags); }
    tag_type get_num_products() const;
    tag_type total_size() const { return operator_table->size(); }
    bool product_is_null(const tag_type t1, const tag_type t2);

private:
    std::shared_ptr<OPTable<Matrix, SymmGroup> > operator_table;

    template <class Map> tag_type duplicates_(Map const & sample);

    std::vector<tag_detail::operator_kind> sign_table;
    pair_map_t product_tags;

    std::vector<tag_type> hermitian;
};

#include "TagHandler.hpp"

#endif
