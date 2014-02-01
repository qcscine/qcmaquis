/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
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

#ifndef TWOU1PG_UTIL_H
#define TWOU1PG_UTIL_H

#include <fstream>
#include <sstream>
#include <algorithm>
#include <boost/tokenizer.hpp>

#include "dmrg/block_matrix/indexing.h"

#include "dmrg/models/op_handler.h"


// TODO: rename to PGIndexConverter
template <class SymmGroup>
class PGDecorator
{
public:
    Index<SymmGroup> operator()(Index<SymmGroup> const & rhs, int irr)
    {
       return rhs;
    }
};

template < >
class  PGDecorator<TwoU1PG>
{
public:
    typedef TwoU1PG::subcharge subcharge;
    Index<TwoU1PG> operator()(Index<TwoU1PG> rhs, subcharge irr)
    {
        for(Index<TwoU1PG>::iterator it = rhs.begin(); it != rhs.end(); ++it)
            if ( (it->first[0] + it->first[1]) % 2 == 0)
                it->first[2] = 0;
            else
                it->first[2] = irr;

        return rhs;
    }
};

template <class SymmGroup>
class PGCharge
{
public:
    void operator()(typename SymmGroup::charge & rhs, int irr)
    { }
};

template < >
class  PGCharge<TwoU1PG>
{
public:
    typedef typename TwoU1PG::subcharge subcharge;
    void operator()(typename TwoU1PG::charge & rhs, subcharge irr)
    {
        rhs[2] = irr;
    }
};


template <class SymmGroup> inline
std::vector<int>
parse_symm(int L, BaseParameters& model)
{
    return std::vector<int>(L, 0);
}

// TODO: This function moved to lattice, remove as soon as possible
template < > inline
std::vector<typename TwoU1PG::subcharge>
parse_symm<TwoU1PG>(int L, BaseParameters& model)
{
    typedef typename TwoU1PG::subcharge subcharge;

    // TODO: pos_t type consistency
    std::vector<int> order(L);
    if (!model.is_set("orbital_order"))
        for (int p = 0; p < L; ++p)
            order[p] = p;
    else
        order = model["orbital_order"].as<std::vector<int> >();

    std::vector<subcharge> irreps(L, 0);
    if (model.is_set("integral_file")) {
        std::ifstream orb_file;
        orb_file.open(model["integral_file"].c_str());
        orb_file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');

        std::string line;
        std::getline(orb_file, line);

        std::vector<std::string> split_line;
        boost::split(split_line, line, boost::is_any_of("="));
        std::vector<subcharge> symm_vec;

        std::replace(split_line[1].begin(), split_line[1].end(), ',', ' ');
        std::istringstream iss(split_line[1]);
        subcharge number;
        while( iss >> number )
            symm_vec.push_back(number-1);

        assert( L == symm_vec.size() );
        for (int p = 0; p < L; ++p)
            irreps[p] = symm_vec[order[p]-1];

        maquis::cout << "Symmetry string (reordered): ";
        std::copy(irreps.begin(), irreps.end(), std::ostream_iterator<subcharge>(maquis::cout, ","));
        maquis::cout << std::endl;

        orb_file.close();
    }
    else
        throw std::runtime_error("\"integral_file\" in model input file is not set\n");

    return irreps;
}

#endif
