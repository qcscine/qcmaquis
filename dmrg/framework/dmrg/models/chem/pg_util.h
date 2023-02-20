/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef TWOU1PG_UTIL_H
#define TWOU1PG_UTIL_H

#include <fstream>
#include <sstream>
#include <algorithm>
#include <boost/tokenizer.hpp>

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/block_matrix/dual_index.h"

// function objects to set point group symmetry in Indices depending on the occupation
// (do not set irrep if empty or doubly occupied)

// TODO: rename to PGIndexConverter
template <class SymmGroup, class = void>
class PGDecorator
{
public:
    DualIndex<SymmGroup> operator()(DualIndex<SymmGroup> const & rhs, int irr)
    {
       return rhs;
    }
};

template <class SymmGroup>
class  PGDecorator<SymmGroup, typename std::enable_if<symm_traits::HasPG<SymmGroup>::value >::type>
{
public:
    typedef typename SymmGroup::subcharge subcharge;
    typedef typename SymmGroup::charge charge;

    DualIndex<SymmGroup> operator()(DualIndex<SymmGroup> rhs, subcharge irr)
    {
        for(typename DualIndex<SymmGroup>::iterator it = rhs.begin(); it != rhs.end(); ++it)
        {
            modify(it->lc, irr);
            modify(it->rc, irr);
        }
        return rhs;
    }

    static void modify(charge & rhs, subcharge irr)
    {
        if ( (SymmGroup::particleNumber(rhs)) % 2 == 0)
            SymmGroup::irrep(rhs) = 0;
        else
            SymmGroup::irrep(rhs) = irr;
    }
};

template <class SymmGroup, class = void>
class PGCharge
{
public:
    typename SymmGroup::charge operator()(typename SymmGroup::charge rhs, int irr)
    {
        return rhs;
    }
};

template <class SymmGroup>
class  PGCharge<SymmGroup, typename std::enable_if<symm_traits::HasPG<SymmGroup>::value >::type>
{
public:
    typedef typename SymmGroup::subcharge subcharge;
    typename SymmGroup::charge operator()(typename SymmGroup::charge rhs, subcharge irr)
    {
        SymmGroup::irrep(rhs) = irr;
        return rhs;
    }
};

template <class SymmGroup, class = void>
class getPG
{
public:
    typename SymmGroup::subcharge operator()(typename SymmGroup::charge rhs)
    {
        return 0;
    }
};

template <class SymmGroup>
class  getPG<SymmGroup, typename std::enable_if<symm_traits::HasPG<SymmGroup>::value >::type>
{
public:
    typedef typename SymmGroup::subcharge subcharge;
    typename SymmGroup::subcharge operator()(typename SymmGroup::charge rhs)
    {
        return SymmGroup::irrep(rhs);
    }
};

#endif
