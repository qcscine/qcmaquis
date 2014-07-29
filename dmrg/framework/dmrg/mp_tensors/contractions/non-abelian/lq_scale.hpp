/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2014-2014 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

#ifndef CONTRACTIONS_SU2_LQ_SCALE_HPP
#define CONTRACTIONS_SU2_LQ_SCALE_HPP


#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/mp_tensors/contractions/non-abelian/gemm.hpp"

namespace SU2 {

    template<class Matrix, class SymmGroup>
    void lq_prescale(Index<SymmGroup> const & physical_i,
                     Index<SymmGroup> const & left_i,
                     Index<SymmGroup> const & right_i,
                     block_matrix<Matrix, SymmGroup> & m)
    {
        typedef std::size_t size_t;
        typedef typename SymmGroup::charge charge;

        ProductBasis<SymmGroup> in_right(physical_i, right_i,
                                         boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                             -boost::lambda::_1, boost::lambda::_2));
        for (size_t block = 0; block < m.n_blocks(); ++block)
        {
            size_t l = left_i.position(m.left_basis_charge(block));
            if(l == left_i.size()) continue;
            charge in_l_charge = left_i[l].first;
            charge in_r_charge = m.right_basis_charge(block);

            for (size_t s = 0; s < physical_i.size(); ++s)
            {
                size_t r = right_i.position(SymmGroup::fuse(in_r_charge, physical_i[s].first));
                if(r == right_i.size()) continue;

                charge out_r_charge = right_i[r].first;

                size_t in_right_offset = in_right(physical_i[s].first, out_r_charge);
                size_t r_length = right_i[r].second;

                Matrix & in_block = m[block];

                typename Matrix::value_type norm_factor = sqrt( (out_r_charge[1]+1.0)/(in_l_charge[1]+1.0) );

                for(size_t c = in_right_offset; c < in_right_offset + r_length; ++c)
                    std::for_each(in_block.col(c).first, in_block.col(c).second, boost::lambda::_1 *= norm_factor);
            }
        }
    }

    template<class Matrix, class SymmGroup>
    void lq_postscale(Index<SymmGroup> const & physical_i,
                     Index<SymmGroup> const & left_i,
                     Index<SymmGroup> const & right_i,
                     block_matrix<Matrix, SymmGroup> & m)
    {
        typedef std::size_t size_t;
        typedef typename SymmGroup::charge charge;
        maquis::cout << "postscale\n";

        ProductBasis<SymmGroup> in_right(physical_i, right_i,
                                         boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                             -boost::lambda::_1, boost::lambda::_2));
        for (size_t block = 0; block < m.n_blocks(); ++block)
        {
            size_t l = left_i.position(m.left_basis_charge(block));
            if(l == left_i.size()) continue;
            charge in_l_charge = left_i[l].first;
            charge in_r_charge = m.right_basis_charge(block);

            for (size_t s = 0; s < physical_i.size(); ++s)
            {
                size_t r = right_i.position(SymmGroup::fuse(in_r_charge, physical_i[s].first));
                if(r == right_i.size()) continue;

                charge out_r_charge = right_i[r].first;

                size_t in_right_offset = in_right(physical_i[s].first, out_r_charge);
                size_t r_length = right_i[r].second;

                Matrix & in_block = m[block];

                typename Matrix::value_type norm_factor = sqrt( (in_l_charge[1]+1.0)/(out_r_charge[1]+1.0) );

                for(size_t c = in_right_offset; c < in_right_offset + r_length; ++c)
                    std::for_each(in_block.col(c).first, in_block.col(c).second, boost::lambda::_1 *= norm_factor);
            }
        }
    }

} // namespace SU2

#endif
