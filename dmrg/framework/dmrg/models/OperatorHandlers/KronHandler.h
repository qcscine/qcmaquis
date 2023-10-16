/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef KRON_HANDLER_H
#define KRON_HANDLER_H

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
#include "TagHandler.h"

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

    tag_type get_num_kron_products() const;

private:

    std::shared_ptr<OPTable<Matrix, SymmGroup> > kronecker_table;
    typename base::pair_map_t kron_tags;
};


#include "KronHandler.hpp"

#endif
