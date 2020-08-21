/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2016 Institute for Theoretical Physics, ETH Zurich
 *                    Laboratory for Physical Chemistry, ETH Zurich
 *               2016-2016 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

#ifndef MAQUIS_DMRG_MPO_TIMES_MPS_HPP
#define MAQUIS_DMRG_MPO_TIMES_MPS_HPP

#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/mpotensor.h"

template <class MPOMatrix, class MPSMatrix, class SymmGroup>
MPSTensor<MPSMatrix, SymmGroup> mpo_times_mps(MPOTensor<MPOMatrix, SymmGroup> const & mpo,
                                              MPSTensor<MPSMatrix, SymmGroup> const & mps,
                                              typename SymmGroup::charge & in_delta)
{
    using MPOTensor_detail::term_descriptor;
    using boost::tuples::get;

    typedef typename SymmGroup::charge charge;
    typedef typename MPSMatrix::value_type value_type;

    mps.make_right_paired();
    block_matrix<MPSMatrix, SymmGroup> const & data = mps.data();

    Index<SymmGroup> const & right_i = mps.col_dim();
    //maquis::cout << "      mps.site_dim: " << mps.site_dim() << std::endl;
    //maquis::cout << "      mps.col_dim : " << mps.col_dim() << std::endl;
    ProductBasis<SymmGroup> right_pb(mps.site_dim(), mps.col_dim(),
                                     boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                        -boost::lambda::_1, boost::lambda::_2));

    term_descriptor<MPOMatrix, SymmGroup, true> access = mpo.at(0,0);
    typename operator_selector<MPOMatrix, SymmGroup>::type const & W = access.op();

    charge W_delta = SymmGroup::fuse(W.basis().right_charge(0), -W.basis().left_charge(0));
    charge out_delta = SymmGroup::fuse(in_delta, W_delta);

    Index<SymmGroup> new_left_i, new_right_i, new_phys_i;
    for (size_t w_block = 0; w_block < W.basis().size(); ++w_block)
    {
        charge phys_in = W.basis().left_charge(w_block);
        if (! mps.site_dim().has(phys_in) ) continue;
        charge phys_out = W.basis().right_charge(w_block);
        new_phys_i.insert(std::make_pair(phys_out, W.basis().right_size(w_block)));
    }

    for (size_t b = 0; b < data.n_blocks(); ++b)
    {
        charge lc = data.basis().left_charge(b);
        charge rc = data.basis().right_charge(b); 

        for (size_t w_block = 0; w_block < W.basis().size(); ++w_block)
        {
            charge phys_in = W.basis().left_charge(w_block);
            charge phys_out = W.basis().right_charge(w_block);

            // add the operator deltas from previous sites to the left charge
            charge out_l_charge = SymmGroup::fuse(lc, in_delta); // unpaired
            if (! charge_detail::physical<SymmGroup>(out_l_charge)) continue;

            charge in_r_charge = SymmGroup::fuse(rc, phys_in); // unpaired
            if (!right_i.has(in_r_charge)) continue; // do we have phys_in in block b?

            charge out_r_charge = SymmGroup::fuse(out_l_charge, phys_out); // unpaired

            if (!new_left_i.has(out_l_charge)) new_left_i.insert(std::make_pair(out_l_charge, data.basis().left_size(b)));
            if (!new_right_i.has(out_r_charge)) new_right_i.insert(std::make_pair(out_r_charge, right_i.size_of_block(in_r_charge)));
        }
    }

    //maquis::cout << "      new_left_i: " << new_left_i << std::endl;
    //maquis::cout << "      new_right_i: " << new_right_i << std::endl;
    //maquis::cout << "      new_phys_i: " << new_phys_i << std::endl;

    ProductBasis<SymmGroup> out_right_pb(new_phys_i, new_right_i,
                                         boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                                                 -boost::lambda::_1, boost::lambda::_2));
    block_matrix<MPSMatrix, SymmGroup> prod;

    for (size_t b = 0; b < data.n_blocks(); ++b)
    {
        charge lc = data.basis().left_charge(b);
        charge rc = data.basis().right_charge(b); 

        for (size_t w_block = 0; w_block < W.basis().size(); ++w_block)
        {
            charge phys_in = W.basis().left_charge(w_block);
            charge phys_out = W.basis().right_charge(w_block);

            charge out_l_charge = SymmGroup::fuse(lc, in_delta); // unpaired
            if (! charge_detail::physical<SymmGroup>(out_l_charge)) continue;

            charge in_r_charge = SymmGroup::fuse(rc, phys_in); // unpaired
            if (!mps.site_dim().has(phys_in)) continue; // do we have phys_in in block b ?
            if (!right_i.has(in_r_charge)) continue; // do we have phys_in in block b?

            charge out_r_charge = SymmGroup::fuse(out_l_charge, phys_out); // unpaired

            // source     -> data[b](·, in_right_offset + 1:rsize)  
            // destination -> prod[o](·, out_right_offset + 1:rsize)
            size_t in_right_offset  = right_pb(phys_in,  in_r_charge); 
            size_t out_right_offset = out_right_pb(phys_out, out_r_charge); 
            size_t l_size = data.basis().left_size(b);
            size_t r_size = right_i.size_of_block(in_r_charge);

            MPSMatrix const & iblock = data[b];

            size_t o = prod.find_block(out_l_charge, out_l_charge); // out_l_charge = out_r_charge right-paired
            if (o == prod.n_blocks())
                o = prod.insert_block(MPSMatrix(l_size, out_right_pb.size(-phys_out, out_r_charge)), out_l_charge, out_l_charge);

            MPSMatrix & oblock = prod[o];

            value_type alfa = access.scale() * W[w_block](0,0);
            //maquis::cout << " access.scale()  ... " << alfa << std::endl;
            //maquis::cout << " W[w_block](0,0) ... " << W[w_block](0,0) << "for block " << w_block << std::endl;
            //maquis::cout << " alfa            ... " << access.scale() << std::endl;
            for(size_t rr = 0; rr < r_size; ++rr)


                maquis::dmrg::detail::iterator_axpy(&iblock(0, in_right_offset + rr),
                                                    &iblock(0, in_right_offset + rr) + l_size,
                                                    &oblock(0, out_right_offset + rr),
                                                    alfa);
        }
    } 
    std::swap(in_delta, out_delta);

    MPSTensor<MPSMatrix, SymmGroup> ret;
    ret.make_right_paired();
    ret.left_i = new_left_i;
    ret.right_i = new_right_i;
    ret.phys_i = new_phys_i;
    swap(ret.data(), prod);

    return ret;
}

#endif
