

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

#ifndef CI_ENCODE_HPP
#define CI_ENCODE_HPP

#include <iostream>

#include <vector>
#include <boost/lexical_cast.hpp>

#include "dmrg/utils/DmrgParameters.h"
#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mps_initializers.h"

/*
    Small function to read determinants from a input file
*/
template<class Matrix, class SymmGroup>
std::vector<std::vector<typename SymmGroup::charge> >
parse_config(std::string file, std::vector<Index<SymmGroup> > const & site_dims)
{
    std::ifstream config_file;
    config_file.open(file.c_str());

    std::vector<std::vector<typename SymmGroup::charge> > configs;

    for (std::string line; std::getline(config_file, line); ) {
        std::vector<std::string> det_coeff;
        boost::split(det_coeff, line, boost::is_any_of(" "));

        std::string det = det_coeff[0];

        if (det.size() != site_dims.size())
            throw std::runtime_error("The determinant length doesn't match the mps length\n");

        std::vector<typename SymmGroup::charge> tmp;
        for (std::size_t i = 0; i < det.size(); ++i) {
            int occ = boost::lexical_cast<int>(det[i]);
            switch(occ) {
                case 4:
                    tmp.push_back(site_dims[i][0].first); // doubly occ
                    break;
                case 3:
                    tmp.push_back(site_dims[i][1].first); // up
                    break;
                case 2:
                    tmp.push_back(site_dims[i][2].first); // down 
                    break;
                case 1:
                    tmp.push_back(site_dims[i][3].first); // empty
                    break;
            }
        }
        configs.push_back(tmp);
    }
    return configs;
}


template <class Matrix, class SymmGroup>
typename Matrix::value_type extract_coefficient(MPS<Matrix, SymmGroup> const & mps, std::vector<typename SymmGroup::charge> const & det)
{
    typedef typename SymmGroup::charge charge;

    if (mps.length() != det.size())
        throw std::runtime_error("extract_coefficient: length of mps != length of basis state\n");

    charge total = std::accumulate(det.begin(), det.end(), SymmGroup::IdentityCharge);
    charge target = mps[mps.length()-1].col_dim()[0].first;
    if (total != target) {
        std::stringstream ss;
        std::copy(det.begin(), det.end(), std::ostream_iterator<charge>(ss, " "));
        ss << " (Has: " << total << ", should be: " << target << ")\n";
        throw std::runtime_error("Determinant has wrong number of up/down electrons: " + ss.str());
    }

    mps[0].make_left_paired();

    block_matrix<Matrix, SymmGroup> const & block0 = mps[0].data();

    charge sector = det[0];
    std::size_t b = block0.left_basis().position(sector);
    if (b == block0.left_basis().size())
        return 0.0;

    Matrix coeff = block0[b];

    for(unsigned p = 1; p < mps.length(); ++p)
    {
        mps[p].make_left_paired();

        charge site_charge = det[p];
        if (! mps[p].site_dim().has(site_charge))
            return 0.0;
        charge left_input = sector;
        sector = SymmGroup::fuse(sector, site_charge);

        block_matrix<Matrix, SymmGroup> const & data = mps[p].data();

        std::size_t b = data.left_basis().position(sector);
        if(b == data.left_basis().size())
            return 0.0;

        Matrix const & sector_matrix = data[b];

        // determine voffset - row offset for site_charge in sector_matrix
        // if permformance is an issue, do not recalculate the ProductBasis on every call
        ProductBasis<SymmGroup> left_pb(mps[p].site_dim(), mps[p].row_dim());
        std::size_t voffset = left_pb(site_charge, left_input);

        // determine vlen and set hlen, the size of the matrix subblock
        unsigned vlen = mps[p-1].col_dim().size_of_block(SymmGroup::fuse(sector, -site_charge));
        unsigned hlen = data.right_basis().size_of_block(sector);

        // extract the subblock
        Matrix site_block(vlen, hlen);
        for (unsigned rowcnt = 0; rowcnt < vlen; ++rowcnt)
            std::copy(sector_matrix.row(voffset+rowcnt).first,
                      sector_matrix.row(voffset+rowcnt).second, site_block.row(rowcnt).first);

        if (coeff.num_cols() != site_block.num_rows())
        {
            maquis::cout << "\nsite charge: " << site_charge << std::endl;
            maquis::cout << "current sector: " << sector << std::endl;
            maquis::cout << "site " << p << ", selected " << vlen << "x" << hlen << " block, offset " << voffset
                         << ", sector_matrix " << sector_matrix.num_rows() << "x" << sector_matrix.num_cols() << std::endl;
            maquis::cout << "  ---> sector_matrix " << std::endl; 
            maquis::cout << sector_matrix << std::endl;
            maquis::cout << "  --->  vlen " << vlen << " hlen " << hlen << std::endl;

            maquis::cout << "coeff" << std::endl << coeff;
            throw std::runtime_error("coeff.num_cols() != vlen\n");
        }

        // multiply subblock with previous subblock
        Matrix tmp(coeff.num_rows(), site_block.num_cols());

        gemm(coeff, site_block, tmp);
        coeff=tmp;
    }

    #ifndef NDEBUG
    if (coeff.num_rows() != 1 || coeff.num_cols() != 1)
        throw std::runtime_error("final coefficient is not a scalar: "
                        + boost::lexical_cast<std::string>(coeff.num_rows()) + "x"
                        + boost::lexical_cast<std::string>(coeff.num_cols()) + "\n");
    #endif

    return coeff(0,0);
}

template <class Matrix, class SymmGroup>
void set_coefficient(MPS<Matrix, SymmGroup> & mps, std::vector<typename SymmGroup::charge> const & det, typename Matrix::value_type coefficient)
{
    typedef typename SymmGroup::charge charge;

    if (mps.length() != det.size())
        throw std::runtime_error("extract_coefficient: length of mps != length of basis state\n");

    charge total = std::accumulate(det.begin(), det.end(), SymmGroup::IdentityCharge),
           target = mps[mps.length()-1].col_dim()[0].first;
    if (total != target) {
        std::stringstream ss;
        std::copy(det.begin(), det.end(), std::ostream_iterator<charge>(ss, " "));
        ss << " (Has: " << total << ", should be: " << target << ")\n";
        throw std::runtime_error("Determinant has wrong number of up/down electrons: " + ss.str());
    }


    Index<SymmGroup> const & phys = mps[0].site_dim();
    std::size_t L = mps.length();

    charge left_sector = det[0];
    std::size_t left_offset = 0;

    mps[0].make_left_paired();
    mps[0].data()[mps[0].data().left_basis().position(left_sector)](0,0) = 1;

    for(std::size_t p=1; p < (L+1)/2; ++p)
    {
        mps[p].make_left_paired();

        charge site_charge = det[p];
        charge input = left_sector;
        left_sector = SymmGroup::fuse(left_sector, site_charge);

        maquis::cout << "coding site " << p << ", site_charge " << site_charge << ", left_sector " << left_sector << std::endl; 

        block_matrix<Matrix, SymmGroup> & data = mps[p].data();
        Matrix & sector_matrix = data[data.left_basis().position(left_sector)];

        // determine voffset - row offset for site_charge in sector_matrix
        std::size_t voffset = 0;
        std::size_t sci = 0;
        while (phys[sci].first > site_charge)
        {
            charge offset_charge = SymmGroup::fuse(left_sector, -phys[sci].first);
            voffset += mps[p-1].col_dim().size_of_block(offset_charge, true);
            sci++;
        }

        // determine vlen and set hlen
        std::size_t vlen = mps[p-1].col_dim().size_of_block(SymmGroup::fuse(left_sector, -site_charge));
        std::size_t hlen = data.right_basis().size_of_block(left_sector);
        maquis::cout << "  site " << p << ", selected " << vlen << "x" << hlen << " block, offset " << voffset
                     << ", sector_matrix " << sector_matrix.num_rows() << "x" << sector_matrix.num_cols() << std::endl;

        // activate subblock
        left_offset += voffset;
        std::size_t col_continue = std::min(left_offset, sector_matrix.num_cols()-1);

        if (left_offset >= sector_matrix.num_rows()) {
            throw std::runtime_error("left_offset too big\n");
        }

        if (p < (L+1)/2-1) {
            maquis::cout << "  activating site " << p << " at " << left_offset << ", " << col_continue << std::endl;
            sector_matrix(left_offset, col_continue) = 1;
            left_offset = col_continue;
        }
        else {
            maquis::cout << "  final left_offset: " << left_offset << std::endl;
        }
    }

    charge right_sector = SymmGroup::fuse(-det[L-1], mps[L-1].col_dim()[0].first);
    std::size_t right_offset = 0;

    mps[L-1].make_right_paired();
    mps[L-1].data()[mps[L-1].data().left_basis().position(right_sector)](0,0) = 1;
    
    for(std::size_t p=L-2; p >= (L+1)/2; --p)
    {
        mps[p].make_right_paired();

        charge site_charge = det[p];
        charge input = right_sector;
        right_sector = SymmGroup::fuse(right_sector, -site_charge);

        maquis::cout << "coding site " << p << ", site_charge " << site_charge << ", right_sector " << right_sector << std::endl; 

        block_matrix<Matrix, SymmGroup> & data = mps[p].data();
        Matrix & sector_matrix = data[data.right_basis().position(right_sector)];

        // determine hoffset - row offset for site_charge in sector_matrix
        std::size_t hoffset = 0;
        std::size_t sci = 0;
        while (phys[sci].first > site_charge)
        {
            charge offset_charge = SymmGroup::fuse(right_sector, phys[sci].first);
            //maquis::cout << "adding offset for phys " << phys[sci].first << ", " << offset_charge << std::endl;
            hoffset += mps[p+1].row_dim().size_of_block(offset_charge, true);
            sci++;
        }

        maquis::cout << "  site " << p << ", offset " << hoffset
                     << ", sector_matrix " << sector_matrix.num_rows() << "x" << sector_matrix.num_cols() << std::endl;

        // activate subblock
        right_offset += hoffset;
        std::size_t row_continue = std::min(right_offset, sector_matrix.num_rows()-1);

        if (right_offset >= sector_matrix.num_cols()) {
            throw std::runtime_error("right_offset too big\n");
        }
        if (row_continue >= sector_matrix.num_rows()) {
            throw std::runtime_error("cannot continue with right_offset: too big\n");
        }
        maquis::cout << "  activating site " << p << " at " << row_continue << ", " << right_offset << std::endl;
        sector_matrix(row_continue, right_offset) = 1;
        right_offset = row_continue;
    }

    if (left_sector != right_sector)
        throw std::runtime_error("code-sector right/left mismatch\n");
    
    // set coefficient
    block_matrix<Matrix, SymmGroup> & center = mps[(L+1)/2 - 1].data();
    Matrix & code_block = center[center.left_basis().position(left_sector)];

    if (left_offset >= code_block.num_rows()) {
        maquis::cout << "code_block " << left_sector << ": " << code_block.num_rows() << "x" << code_block.num_cols() << std::endl;
        throw std::runtime_error("center left_offset too big: " + boost::lexical_cast<std::string>(left_offset) + "\n");
    }
    if (right_offset >= code_block.num_cols())
        throw std::runtime_error("center right_offset too big\n");

    //code_block[code_block.left_basis().position(left_sector)](left_offset, right_offset) = coefficient;
    code_block(left_offset, right_offset) = coefficient;

    maquis::cout << std::endl;
}

#endif


