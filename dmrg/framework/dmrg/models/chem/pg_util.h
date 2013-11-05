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


template <class SymmGroup>
class PGDecorator
{
public:
    typedef int irrep_t;
    Index<SymmGroup> operator()(Index<SymmGroup> const & rhs, irrep_t irr)
    {
       return rhs;
    }
};

template < >
class  PGDecorator<TwoU1PG>
{
public:
    typedef int irrep_t;
    Index<TwoU1PG> operator()(Index<TwoU1PG> rhs, irrep_t irr)
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
    typedef typename PGDecorator<SymmGroup>::irrep_t irrep_t;
    void operator()(typename SymmGroup::charge & rhs, irrep_t irr)
    { }
};

template < >
class  PGCharge<TwoU1PG>
{
public:
    typedef PGDecorator<TwoU1PG>::irrep_t irrep_t;
    void operator()(TwoU1PG::charge & rhs, irrep_t irr)
    {
        rhs[2] = irr;
    }
};

template <class SymmGroup>
class PGChargeTransfer
{
public:
    typedef typename PGDecorator<SymmGroup>::irrep_t irrep_t;
    void operator()(typename SymmGroup::charge const & ref, typename SymmGroup::charge & target)
    { }
};

template < >
class  PGChargeTransfer<TwoU1PG>
{
public:
    typedef PGDecorator<TwoU1PG>::irrep_t irrep_t;
    void operator()(TwoU1PG::charge const & ref, TwoU1PG::charge & target)
    {
        if ( (target[0] + target[1]) % 2 )
            target[2] = ref[2];
    }
};

template <class SymmGroup>
class PGSensible
{
public:
    typedef typename PGDecorator<SymmGroup>::irrep_t irrep_t;
    bool operator()(typename SymmGroup::charge & rhs)
    {
        return true;
    }
};

template < >
class  PGSensible<TwoU1PG>
{
public:
    typedef PGDecorator<TwoU1PG>::irrep_t irrep_t;
    bool operator()(TwoU1PG::charge & rhs)
    {
        irrep_t irr = rhs[2];
        if (irr > 7 || irr < 0)
            return false;
        else
            return true;
    }
};

template <class SymmGroup> inline 
std::vector<typename PGDecorator<SymmGroup>::irrep_t>
parse_symm(int L, BaseParameters& model)
{
    return std::vector<typename PGDecorator<SymmGroup>::irrep_t>(L, 0);
}

template < > inline
std::vector<PGDecorator<TwoU1PG>::irrep_t>
parse_symm<TwoU1PG>(int L, BaseParameters& model)
{
    typedef PGDecorator<TwoU1PG>::irrep_t irrep_t;

    // TODO: pos_t type consistency
    std::vector<int> order(L);
    if (!model.is_set("orbital_order"))
        for (int p = 0; p < L; ++p)
            order[p] = p;
    else
        order = model["orbital_order"];

    std::vector<irrep_t> irreps(L, 0);
    if (model.is_set("integral_file")) {
        std::ifstream orb_file;
        orb_file.open(model["integral_file"].c_str());
        orb_file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');

        std::string line;
        std::getline(orb_file, line);

        std::vector<std::string> split_line;
        boost::split(split_line, line, boost::is_any_of("="));
        std::vector<irrep_t> symm_vec;

        std::replace(split_line[1].begin(), split_line[1].end(), ',', ' ');
        std::istringstream iss(split_line[1]);
        irrep_t number;
        while( iss >> number )
            symm_vec.push_back(number-1);

        assert( L == symm_vec.size() );
        for (int p = 0; p < L; ++p)
            irreps[p] = symm_vec[order[p]-1];

        maquis::cout << "Symmetry string (reordered): ";
        std::copy(irreps.begin(), irreps.end(), std::ostream_iterator<irrep_t>(maquis::cout, ","));
        maquis::cout << std::endl;

        orb_file.close();
    }
    else
        throw std::runtime_error("need an integral (FCIDUMP) file to load orbital symmetries\n");

    return irreps;
}

#endif
