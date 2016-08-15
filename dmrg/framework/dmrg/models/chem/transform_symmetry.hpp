/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2015 Laboratory of Physical Chemistry, ETH Zurich
 *               2015-2015 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

#ifndef TRANSFORM_SYMMETRY_HPP
#define TRANSFORM_SYMMETRY_HPP

#include <boost/mpl/if.hpp>

#include "dmrg/block_matrix/symmetry/gsl_coupling.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/models/chem/util.h"

namespace transform_detail
{

    template<class SymmIn, class SymmOut>
    std::vector<typename SymmOut::charge> transform_charge(typename SymmIn::charge cin)
    {
        typedef typename SymmIn::subcharge subcharge;

        subcharge sz_min = -SymmIn::spin(cin);
        subcharge sz_max =  SymmIn::spin(cin);
        subcharge N      =  SymmIn::particleNumber(cin);

        std::vector<typename SymmOut::charge> ret;

        // sz = 2 * S_z
        for (subcharge sz = sz_min; sz <= sz_max; sz += 2)
        {
            typename SymmOut::charge tcharge;
            tcharge[0] = (N + sz)/2;
            tcharge[1] = (N - sz)/2;
            tcharge = PGCharge<SymmOut>()(tcharge, getPG<SymmIn>()(cin));
            ret.push_back(tcharge);
        }

        return ret;
    }
}

template<class Matrix, class SymmIn, class SymmOut>
void transform_site(MPSTensor<Matrix, SymmIn> const & mps_in,
                    MPSTensor<Matrix, SymmOut> & mps_out)
{
    typedef std::size_t size_t;
    typedef typename SymmIn::charge charge;
    typedef typename SymmOut::charge out_charge;

    Index<SymmIn> const & physical_i = mps_in.site_dim();
    Index<SymmIn> const & left_i = mps_in.row_dim();
    Index<SymmIn> const & right_i = mps_in.col_dim();

    Index<SymmOut> const & physical_i_out = mps_out.site_dim();
    Index<SymmOut> const & left_i_out = mps_out.row_dim();
    Index<SymmOut> const & right_i_out = mps_out.col_dim();

    block_matrix<Matrix, SymmIn> const & m1 = mps_in.data();
    block_matrix<Matrix, SymmOut> & m2 = mps_out.data();

    ProductBasis<SymmIn> in_left_pb(physical_i, left_i);

    // data for the layout of the output MPS/block_matrix
    // each 2u1 sector (out_charge) contains a SU2 Index to describe the SU2 blocks within
    // the larger 2u1 block
    typedef std::map<out_charge, Index<SymmIn> > subsector_map_t;
    subsector_map_t left_subblocks, right_subblocks;
    Index<SymmOut> new_left_i, new_right_i;

    for (int pass = 0; pass < 2; ++pass)
    {
        // allocate blocks in output 2u1 block_matrix
        if (pass == 1)
        {
            for (typename subsector_map_t::iterator it = left_subblocks.begin(); it != left_subblocks.end(); ++it)
                new_left_i.insert(std::make_pair(it->first, (it->second).sum_of_sizes()));

            for (typename subsector_map_t::iterator it = right_subblocks.begin(); it != right_subblocks.end(); ++it)
                new_right_i.insert(std::make_pair(it->first, (it->second).sum_of_sizes()));

            mps_out = MPSTensor<Matrix, SymmOut>(physical_i_out, new_left_i, new_right_i, false, 0.);
        }    

        ProductBasis<SymmOut> out_left_pb(physical_i_out, new_left_i);

    for (size_t block = 0; block < m1.n_blocks(); ++block)
    {
        size_t r = right_i.position(m1.basis().right_charge(block));
        if(r == right_i.size()) continue;
        charge in_r_charge = right_i[r].first;
        charge in_l_charge_paired = m1.basis().left_charge(block);

        for (size_t s = 0; s < physical_i.size(); ++s)
        {
            size_t l = left_i.position(SymmIn::fuse(m1.basis().left_charge(block), -physical_i[s].first));
            if(l == left_i.size()) continue;

            charge in_l_charge = left_i[l].first;

            // transform one SU2 charge to corresponding 2U1 charges
            std::vector<out_charge> l_sectors = transform_detail::transform_charge<SymmIn, SymmOut>(in_l_charge);
            std::vector<out_charge> r_sectors = transform_detail::transform_charge<SymmIn, SymmOut>(in_r_charge);

            // form pairs from matching right and left sectors
            std::vector<std::pair<out_charge, out_charge> > sectors;
            for (typename std::vector<out_charge>::const_iterator it1 = l_sectors.begin(); it1 != l_sectors.end(); ++it1)
                for (typename std::vector<out_charge>::const_iterator it2 = r_sectors.begin(); it2 != r_sectors.end(); ++it2)
                    for (typename Index<SymmOut>::const_iterator itp = physical_i_out.begin(); itp != physical_i_out.end(); ++itp)
                        if (SymmOut::fuse(*it1, itp->first) == *it2)
                            sectors.push_back(std::make_pair(*it1, *it2));

            // record positions of the SU2 blocks within the larger 2U1 blocks
            if (pass == 0)
            {
                // insert the source sector into a non-paired target symmetry block_matrix
                for (typename std::vector<std::pair<out_charge, out_charge> >::const_iterator it = sectors.begin(); it != sectors.end(); ++it)
                {
                    out_charge leftc = it->first, rightc = it->second, physc = SymmOut::fuse(leftc, -rightc);

                    if ( !left_i_out.has(leftc) || !right_i_out.has(rightc) )
                        continue;

                    if (!left_subblocks[leftc].has(in_l_charge))
                        left_subblocks[leftc].insert(left_i[l]);

                    if (!right_subblocks[rightc].has(in_r_charge))
                        right_subblocks[rightc].insert(right_i[r]);

                    assert(left_subblocks[leftc].size_of_block(in_l_charge) == left_i[l].second);
                    assert(right_subblocks[rightc].size_of_block(in_r_charge) == right_i[r].second);
                }
            }
            // transfer the blocks
            else
            {
                std::size_t in_left_offset = in_left_pb(physical_i[s].first, left_i[l].first);
                std::size_t ldim = left_i[l].second;
                Matrix const & iblock = m1[block];
                Matrix source_block(ldim, right_i[r].second);

                // extract source block
                for (std::size_t ci = 0; ci < num_cols(iblock); ++ci)
                    std::copy(iblock.col(ci).first + in_left_offset, iblock.col(ci).first + in_left_offset + ldim, source_block.col(ci).first);

                for (typename std::vector<std::pair<out_charge, out_charge> >::const_iterator it = sectors.begin(); it != sectors.end(); ++it)
                {
                    out_charge leftc = it->first, rightc = it->second, physc = SymmOut::fuse(-leftc, rightc);

                    if (!m2.has_block(rightc, rightc))
                        continue;
                    Matrix & current_block = m2(rightc, rightc); // left_paired

                    std::size_t  out_left_offset_2u1 = out_left_pb(physc, leftc);
                    std::size_t  out_left_offset_su2 = left_subblocks[leftc].position(std::make_pair(in_l_charge, 0));
                    std::size_t out_right_offset_su2 = right_subblocks[rightc].position(std::make_pair(in_r_charge, 0));

                    int l1 = SymmIn::spin(in_l_charge), l2 = std::abs(SymmIn::spin(physical_i[s].first)), l3 = SymmIn::spin(in_r_charge);
                    int m1 = leftc[0] - leftc[1], m2 = physc[0] - physc[1], m3 = rightc[0] - rightc[1];
                    double clebsch_gordan = pow(-1.0,(l1-l2+m3)/2)*sqrt(l3+1.0)*gsl_sf_coupling_3j(l1,l2,l3,m1,m2,-m3);

                    for (std::size_t ci = 0; ci < num_cols(source_block); ++ci)
                        std::transform(source_block.col(ci).first, source_block.col(ci).second,
                                       current_block.col(ci + out_right_offset_su2).first + out_left_offset_2u1 + out_left_offset_su2,
                                       boost::lambda::_1*clebsch_gordan);
                }
            }
        } // SU2 input physical_i
    } // m1 block
    } // pass
}

template <class Matrix, class SymmGroup, class = void>
struct transform_mps
{
    typedef typename boost::mpl::if_<symm_traits::HasPG<SymmGroup>, TwoU1PG, TwoU1>::type SymmOut;

    void operator()(MPS<Matrix, SymmGroup> const & mps_in, MPS<Matrix, SymmOut> & mps_out)
    {}
};

template <class Matrix, class SymmGroup>
struct transform_mps<Matrix, SymmGroup, typename boost::enable_if<symm_traits::HasSU2<SymmGroup> >::type>
{
    typedef typename boost::mpl::if_<symm_traits::HasPG<SymmGroup>, TwoU1PG, TwoU1>::type SymmOut;

    MPS<Matrix, SymmOut> operator()(MPS<Matrix, SymmGroup> mps_in, int Nup, int Ndown)
    {
        BaseParameters parms = chem_detail::set_2u1_parameters(mps_in.size(), Nup, Ndown);
        parms.set("init_bond_dimension", 1000);

        // determine the irreps per site
        std::string site_types;
        for (Lattice::pos_t p = 0; p < mps_in.size(); ++p)
            for (std::size_t i = 0; i < mps_in[p].site_dim().size(); ++i)
            {
                if (SymmGroup::particleNumber(mps_in[p].site_dim()[i].first) % 2 != 0)
                {
                    site_types += boost::lexical_cast<std::string>(getPG<SymmGroup>()(mps_in[p].site_dim()[i].first)) + ",";
                    break;
                }
                if (i == mps_in[p].site_dim().size())
                    site_types += "0,";
            }
        parms.set("site_types", site_types);

        Lattice lat(parms);
        Model<Matrix, SymmOut> model(lat, parms);

        typename SymmOut::charge initc;
        initc[0] = Nup;
        initc[1] = Ndown;
        initc = PGCharge<SymmOut>()(initc, getPG<SymmGroup>()(mps_in[mps_in.size()-1].col_dim()[0].first));

        std::vector<Index<SymmOut> > site_bases;
        for (int i = 0; i <= lat.maximum_vertex_type(); ++i)
            site_bases.push_back(model.phys_dim(i));

        const_mps_init<Matrix, SymmOut> mpsinit(parms, site_bases, initc, parms["site_types"]);
        MPS<Matrix, SymmOut> mps_out(mps_in.size(), mpsinit);

        // clean the input MPS, ensure consistent indices across bonds
        for (Lattice::pos_t p = 0; p < mps_out.length()-1; ++p)
        {
            mps_in[p].make_left_paired();
            mps_in[p+1].make_right_paired();
            block_matrix<Matrix, SymmGroup> bm1 = mps_in[p].data(), bm2 = mps_in[p+1].data();

            for (size_t b = 0; b < bm1.n_blocks(); ++b)
            {
                typename SymmGroup::charge c = bm1.basis().right_charge(b);
                if (! bm2.basis().has(c, c))
                    bm1.remove_block(b--);
            }
            for (size_t b = 0; b < bm2.n_blocks(); ++b)
            {
                typename SymmGroup::charge c = bm2.basis().left_charge(b);
                if (! bm1.basis().has(c, c))
                    bm2.remove_block(b--);
            }

            assert(bm1.right_basis() == bm2.left_basis());

            mps_in[p].replace_left_paired(bm1);
            mps_in[p+1].replace_right_paired(bm2);
        }

        for (Lattice::pos_t p = 0; p < mps_out.length(); ++p)
        {
            mps_in[p].make_left_paired();
            transform_site(mps_in[p], mps_out[p]);
        }

        return mps_out;
    }
};

#endif
