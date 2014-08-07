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
#include "dmrg/block_matrix/dual_index.h"

#include "dmrg/models/op_handler.h"


// TODO: rename to PGIndexConverter
template <class SymmGroup>
class PGDecorator
{
public:
    PGDecorator(BaseParameters & parms) {}
    PGDecorator(bool su2_) {}

    DualIndex<SymmGroup> operator()(DualIndex<SymmGroup> const & rhs, int irr)
    {
       return rhs;
    }
};

template < >
class  PGDecorator<TwoU1PG>
{
public:

    PGDecorator(BaseParameters & parms) : su2(false)
    {
        if (parms["MODEL"] == "quantum_chemistry_SU2")
            su2 = true;
    }

    PGDecorator(bool su2_) : su2(su2_) {}

    typedef TwoU1PG::subcharge subcharge;
    DualIndex<TwoU1PG> operator()(DualIndex<TwoU1PG> rhs, subcharge irr)
    {
        if(su2)
            for(DualIndex<TwoU1PG>::iterator it = rhs.begin(); it != rhs.end(); ++it)
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

        else
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

private:
    bool su2;
};

template <class SymmGroup>
class PGCharge
{
public:
    typename SymmGroup::charge operator()(typename SymmGroup::charge rhs, int irr)
    { 
        return rhs;
    }
};

template < >
class  PGCharge<TwoU1PG>
{
public:
    typedef TwoU1PG::subcharge subcharge;
    TwoU1PG::charge operator()(TwoU1PG::charge rhs, subcharge irr)
    {
        rhs[2] = irr;
        return rhs;
    }
};

#endif
