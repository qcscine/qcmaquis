/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
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
#include "dmrg/block_matrix/dual_index.h"

#include "dmrg/models/op_handler.h"

// function objects to set point group symmetry in Indices depending on the occupation
// (do not set irrep if empty or doubly occupied)

// TODO: rename to PGIndexConverter
template <class SymmGroup>
class PGDecorator
{
public:
    DualIndex<SymmGroup> operator()(DualIndex<SymmGroup> const & rhs, int irr)
    {
       return rhs;
    }
};

template < >
class  PGDecorator<TwoU1PG>
{
public:
    typedef TwoU1PG::subcharge subcharge;
    DualIndex<TwoU1PG> operator()(DualIndex<TwoU1PG> rhs, subcharge irr)
    {
        for(DualIndex<TwoU1PG>::iterator it = rhs.begin(); it != rhs.end(); ++it)
        {
            if ( (it->lc[0] + it->lc[1]) % 2 == 0)
                it->lc[2] = 0;
            else
                it->lc[2] = irr;

            if ( (it->rc[0] + it->rc[1]) % 2 == 0)
                it->rc[2] = 0;
            else
                it->rc[2] = irr;
        }

        return rhs;
    }
};

template < >
class  PGDecorator<TwoU1LPG>
{
public:
    typedef TwoU1LPG::subcharge subcharge;
    DualIndex<TwoU1LPG> operator()(DualIndex<TwoU1LPG> rhs, subcharge irr)
    {
        for(DualIndex<TwoU1LPG>::iterator it = rhs.begin(); it != rhs.end(); ++it)
		{
            if ( (it->lc[0] + it->lc[1]) % 2 == 0)
                it->lc[2] = 0;
            else
                it->lc[2] = irr;

            if ( (it->rc[0] + it->rc[1]) % 2 == 0)
                it->rc[2] = 0;
            else
                it->rc[2] = irr;
		}
        return rhs;
    }
};

template < >
class  PGDecorator<U1LPG>
{
public:
    typedef U1LPG::subcharge subcharge;
    DualIndex<U1LPG> operator()(DualIndex<U1LPG> rhs, subcharge irr)
    {
        for(DualIndex<U1LPG>::iterator it = rhs.begin(); it != rhs.end(); ++it)
		{
            if ( it->lc[0] % 2 == 0)
                it->lc[1] = 0;
            else
                it->lc[1] = irr;

            if ( it->rc[0] % 2 == 0)
                it->rc[1] = 0;
            else
                it->rc[1] = irr;
		}
        return rhs;
    }
};

template < >
class  PGDecorator<SU2U1PG>
{
public:
    typedef SU2U1PG::subcharge subcharge;
    DualIndex<SU2U1PG> operator()(DualIndex<SU2U1PG> rhs, subcharge irr)
    {
        for(DualIndex<SU2U1PG>::iterator it = rhs.begin(); it != rhs.end(); ++it)
        {
            if ( (it->lc[0]) % 2 == 0)
                it->lc[2] = 0;
            else
                it->lc[2] = irr;

            if ( (it->rc[0]) % 2 == 0)
                it->rc[2] = 0;
            else
                it->rc[2] = irr;
        }
        return rhs;
    }
};

//////////////////////////////////////////////////

template <class SymmGroup, class PGTag>
class PGCharge_
{
public:
    typename SymmGroup::charge operator()(typename SymmGroup::charge rhs, int irr)
    { 
        return rhs;
    }
};

template <class SymmGroup>
class  PGCharge_<SymmGroup, symm_traits::PGat1>
{
public:
    typedef typename SymmGroup::subcharge subcharge;
    typename SymmGroup::charge operator()(typename SymmGroup::charge rhs, subcharge irr)
    {
        rhs[1] = irr;
        return rhs;
    }
};

template <class SymmGroup>
class  PGCharge_<SymmGroup, symm_traits::PGat2>
{
public:
    typedef typename SymmGroup::subcharge subcharge;
    typename SymmGroup::charge operator()(typename SymmGroup::charge rhs, subcharge irr)
    {
        rhs[2] = irr;
        return rhs;
    }
};

template <class SymmGroup>
class PGCharge
{
public:
    typename SymmGroup::charge operator()(typename SymmGroup::charge rhs, int irr)
    { 
        return PGCharge_<SymmGroup, typename symm_traits::PGType<SymmGroup>::type>()(rhs, irr);
    }
};

#endif
