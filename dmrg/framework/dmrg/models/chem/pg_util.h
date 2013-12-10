/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2013-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
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
        order = model["orbital_order"].template as<std::vector<int> >();

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
        std::copy(irreps.begin(), irreps.end(), std::ostream_iterator<subcharge>(std::cout, ","));
        maquis::cout << std::endl;

        orb_file.close();
    }
    else
        throw std::runtime_error("\"integral_file\" in model input file is not set\n");

    return irreps;
}

#endif
