/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

#ifndef MAQUIS_DMRG_MODELS_OP_HANDLER_H
#define MAQUIS_DMRG_MODELS_OP_HANDLER_H

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
class OPTable : public std::vector<typename operator_selector<Matrix, SymmGroup>::type>
{
public:
    typedef tag_detail::tag_type tag_type;
    typedef typename operator_selector<Matrix, SymmGroup>::type op_t;

private:
    typedef typename Matrix::value_type mvalue_type;

public:
    tag_type register_op(op_t const & op_);
    std::pair<tag_type, mvalue_type> checked_register(op_t const& sample);
};

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
    bool is_null(const tag_type t1, const tag_type t2);

private:
    std::shared_ptr<OPTable<Matrix, SymmGroup> > operator_table;

    template <class Map> tag_type duplicates_(Map const & sample);

    std::vector<tag_detail::operator_kind> sign_table;
    pair_map_t product_tags;

    std::vector<tag_type> hermitian;
};

template <class Matrix, class SymmGroup>
class KronHandler : public TagHandler<Matrix, SymmGroup>
{
    typedef TagHandler<Matrix, SymmGroup> base;
    typedef typename OPTable<Matrix, SymmGroup>::tag_type tag_type;
    typedef typename base::op_t op_t;

public:

    KronHandler(std::shared_ptr<OPTable<Matrix, SymmGroup> > tbl_)
    :   base(tbl_)
      , kronecker_table(new OPTable<Matrix, SymmGroup>()) { }

    tag_type get_kron_tag(Index<SymmGroup> const & phys_i1, Index<SymmGroup> const & phys_i2, tag_type t1, tag_type t2,
                          SpinDescriptor<typename symm_traits::SymmType<SymmGroup>::type> lspin,
                          SpinDescriptor<typename symm_traits::SymmType<SymmGroup>::type> mspin,
                          SpinDescriptor<typename symm_traits::SymmType<SymmGroup>::type> rspin);

    typename OPTable<Matrix, SymmGroup>::value_type & get_op(tag_type i) { return (*kronecker_table)[i]; }
    typename OPTable<Matrix, SymmGroup>::value_type const & get_op(tag_type i) const { return (*kronecker_table)[i]; }

    std::shared_ptr<OPTable<Matrix, SymmGroup> > get_kronecker_table() { return kronecker_table; }

    /* Diagnostics *************************************/
    tag_type get_num_kron_products() const;

private:
    std::shared_ptr<OPTable<Matrix, SymmGroup> > kronecker_table;
    typename base::pair_map_t kron_tags;
};


#include "dmrg/models/op_handler.hpp"

#endif
