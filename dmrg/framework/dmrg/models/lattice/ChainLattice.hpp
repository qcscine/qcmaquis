/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2012-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 *                            Sebastian Keller <sebkelle@phys.ethz.ch>
 *               2020- by Robin Feldmann <robinfe@phys.chem.ethz.ch>
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

#ifndef MAQUIS_DMRG_lattice_hpp
#define MAQUIS_DMRG_lattice_hpp

#include "dmrg/models/lattice/lattice.h"
#include <sstream>
#include <vector>
#include <set>
#include <boost/lexical_cast.hpp>
#include <boost/lambda/lambda.hpp>
#include <numeric>
#include "dmrg/utils/BaseParameters.h"

class ChainLattice : public lattice_impl
{
public:
    typedef lattice_impl::pos_t pos_t;

    ChainLattice (BaseParameters & parms, bool pbc_=false)
    : L(parms["L"])
    , a(parms["a"])
    , pbc(pbc_)
    { }

    ChainLattice (int L_, bool pbc_=false, double a_=1.)
    : L(L_)
    , a(a_)
    , pbc(pbc_)
    { }

    std::vector<pos_t> forward(pos_t i) const
    {
        std::vector<pos_t> ret;
        if (i < L-1)
            ret.push_back(i+1);
        if (pbc && i == L-1)
            ret.push_back(0);
        return ret;
    }
    std::vector<pos_t> all(pos_t i) const
    {
        std::vector<pos_t> ret;
        if (i < L-1)
            ret.push_back(i+1);
        if (i > 0)
            ret.push_back(i-1);
        if (pbc && i == L-1)
            ret.push_back(0);
        if (pbc && i == 0)
            ret.push_back(L-1);
        return ret;
    }

    boost::any get_prop_(std::string const & property, std::vector<pos_t> const & pos) const
    {
        if (property == "label" && pos.size() == 1)
            return boost::any( site_label(pos[0]) );
        else if (property == "label" && pos.size() == 2)
            return boost::any( bond_label(pos[0], pos[1]) );
        else if (property == "type" && pos.size() == 1)
            return boost::any( 0 );
        else if (property == "type" && pos.size() == 2)
            return boost::any( 0 );
        else if (property == "x" && pos.size() == 1)
            return boost::any( a * pos[0] );
        else if (property == "at_open_boundary" && pos.size() == 1)
            return boost::any( (!pbc) && (pos[0]==0 || pos[0]==L-1) );
        else if (property == "at_open_left_boundary" && pos.size() == 1)
            return boost::any( (!pbc) && pos[0]==0 );
        else if (property == "at_open_right_boundary" && pos.size() == 1)
            return boost::any( (!pbc) && pos[0]==L-1 );
        else if (property == "wraps_pbc" && pos.size() == 2)
            return boost::any( (pos[0] < pos[1]) );
        else if (property == "NumTypes")
            return boost::any( 1 );
        else if (property == "ParticleType" && pos.size() == 1)
            return boost::any( 0 );
        else {
            std::ostringstream ss;
            ss << "No property '" << property << "' with " << pos.size() << " points implemented.";
            throw std::runtime_error(ss.str());
            return boost::any();
        }
    }

    pos_t size() const
    {
        return L;
    }

    int maximum_vertex_type() const
    {
        return 0;
    }

private:

    std::string site_label (int i) const
    {
        return "( " + boost::lexical_cast<std::string>(a * i) + " )";
    }

    std::string bond_label (int i, int j) const
    {
        return (  "( " + boost::lexical_cast<std::string>(a * i) + " )"
                + " -- "
                + "( " + boost::lexical_cast<std::string>(a * j) + " )");
    }

private:
    int L;
    double a;
    bool pbc;

};

#endif
