/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2013-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *
 *
 *****************************************************************************/

#ifndef CI_ENCODE_HPP
#define CI_ENCODE_HPP

#include <iostream>

#include <vector>
#include <boost/lexical_cast.hpp>

#include "dmrg/utils/DmrgParameters2.h"
#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mps_initializers.h"


template <class Matrix, class SymmGroup>
typename Matrix::value_type extract_coefficient(MPS<Matrix, SymmGroup> const & mps, std::vector<typename SymmGroup::charge> const & det)
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
    block_matrix<Matrix, SymmGroup> const & block0 = mps[0].data();

    charge sector = det[0];
    Matrix coeff = block0[block0.left_basis().position(sector)];

    mps[0].make_left_paired();
    for(unsigned p=1; p<mps.length(); ++p)
    {
        mps[p].make_left_paired();

        charge site_charge = det[p];
        charge input = sector;
        sector = SymmGroup::fuse(sector, site_charge);

        //maquis::cout << "site " << p << ", site_charge " << site_charge << ", sector " << sector << std::endl; 

        block_matrix<Matrix, SymmGroup> const & data = mps[p].data();
        Matrix const & sector_matrix = data[data.left_basis().position(sector)];

        // determine voffset - row offset for site_charge in sector_matrix
        unsigned voffset = 0;
        unsigned sci = 0;
        while (phys[sci].first > site_charge)
        {
            charge offset_charge = SymmGroup::fuse(sector, -phys[sci].first);
            voffset += mps[p-1].col_dim().size_of_block(offset_charge, true);
            sci++;
        }

        // determine vlen and set hlen
        unsigned vlen = mps[p-1].col_dim().size_of_block(SymmGroup::fuse(sector, -site_charge));
        unsigned hlen = data.right_basis().size_of_block(sector);

        //maquis::cout << "site " << p << ", selected " << vlen << "x" << hlen << " block, offset " << voffset
        //             << ", sector_matrix " << sector_matrix.num_rows() << "x" << sector_matrix.num_cols() << std::endl;

        // extract the subblock
        Matrix site_block(vlen, hlen);
        for (unsigned rowcnt = 0; rowcnt < vlen; ++rowcnt)
            std::copy(sector_matrix.row(voffset+rowcnt).first,
                      sector_matrix.row(voffset+rowcnt).second, site_block.row(rowcnt).first);

        if (coeff.num_cols() != site_block.num_rows())
            throw std::runtime_error("coeff.num_cols() != site_block.num_rows()\n");

        // multiply subblock with previous subblock
        Matrix tmp(coeff.num_rows(), site_block.num_cols());
        gemm(coeff, site_block, tmp);
        coeff=tmp;
        
        //maquis::cout << std::endl;
    }

    if (coeff.num_rows() != 1 || coeff.num_cols() != 1)
        throw std::runtime_error("final coefficient is not a scalar: "
                        + boost::lexical_cast<std::string>(coeff.num_rows()) + "x"
                        + boost::lexical_cast<std::string>(coeff.num_cols()) + "\n");

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
