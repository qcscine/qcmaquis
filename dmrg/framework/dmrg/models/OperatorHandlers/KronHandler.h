/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *               2021 by Alberto Baiardi <abaiardi@ethz.ch>
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
