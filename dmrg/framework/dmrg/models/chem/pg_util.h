/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2013-2015 by Sebastian Keller <sebkelle@phys.ethz.ch>
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
class  PGDecorator<SymmGroup, typename boost::enable_if<symm_traits::HasPG<SymmGroup> >::type>
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
class  PGCharge<SymmGroup, typename boost::enable_if<symm_traits::HasPG<SymmGroup> >::type>
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
class  getPG<SymmGroup, typename boost::enable_if<symm_traits::HasPG<SymmGroup> >::type>
{
public:
    typedef typename SymmGroup::subcharge subcharge;
    typename SymmGroup::subcharge operator()(typename SymmGroup::charge rhs)
    {
        return SymmGroup::irrep(rhs);
    }
};

#endif
