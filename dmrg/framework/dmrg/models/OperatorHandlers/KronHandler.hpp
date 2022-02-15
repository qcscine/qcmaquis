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

#ifndef KRON_HANDLER_HPP
#define KRON_HANDLER_HPP

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
