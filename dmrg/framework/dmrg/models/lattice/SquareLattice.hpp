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

#ifndef MAQUIS_DMRG_SQUARE_LATTICE
#define MAQUIS_DMRG_SQUARE_LATTICE

#include "dmrg/models/lattice/lattice.h"
#include <sstream>
#include <vector>
#include <set>
#include <boost/lexical_cast.hpp>
#include <boost/lambda/lambda.hpp>
#include <numeric>
#include "dmrg/utils/BaseParameters.h"

class SquareLattice : public lattice_impl
{
public:
    SquareLattice(BaseParameters & parms)
    : L_(parms["L"])
    , W_(parms["W"])
    , a(parms["a"])
    { }

    /*
     0 4  8 12
     1 5  9 13
     2 6 10 14
     3 7 11 15
     */
    std::vector<int> forward(int p) const
    {
        std::vector<int> ret;
        if (p+1 < L_*W_ && (p+1) % W_ != 0)
            ret.push_back(p+1);
        if (p+W_ < L_*W_)
            ret.push_back(p+W_);

        //        maquis::cout << p << " -> ";
        //        std::copy(ret.begin(), ret.end(), std::ostream_iterator<int>(maquis::cout, " "));
        //        maquis::cout << std::endl;

        return ret;
    }

    std::vector<int> all(int p) const
    {
        std::vector<int> ret = forward(p);
        if (p >= 1 && p % W_ != 0)
            ret.push_back(p-1);
        if (p >= W_)
            ret.push_back(p-W_);

        return ret;
    }

    int size() const { return L_*W_; }

    int maximum_vertex_type() const
    {
        return 0;
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
            return boost::any( x(pos[0]) );
        else if (property == "y" && pos.size() == 1)
            return boost::any( y(pos[0]) );
        else if (property == "wraps_pbc" && pos.size() == 2)
            return boost::any( false );
        else if (property == "NumTypes")
            return boost::any( 1 );
        else if (property == "ParticleType")
            return boost::any( 0 );
        else {
            std::ostringstream ss;
            ss << "No property '" << property << "' with " << pos.size() << " points implemented.";
            throw std::runtime_error(ss.str());
            return boost::any();
        }
    }

private:

    double x (int i) const
    { return a * int(i/W_); }
    double y (int i) const
    { return a * (i%W_); }

    std::string site_label (int i) const
    {
        return "( " + ( boost::lexical_cast<std::string>(x(i))
                       + "," + boost::lexical_cast<std::string>(y(i)) ) + " )";
    }

    std::string bond_label (int i, int j) const
    {
        return (  "( " + ( boost::lexical_cast<std::string>(x(i))
                          + "," + boost::lexical_cast<std::string>(y(i)) ) + " )"
                + " -- "
                + "( " + ( boost::lexical_cast<std::string>(x(j))
                          + "," + boost::lexical_cast<std::string>(y(j)) ) + " )" );
    }


    int L_, W_;
    double a;
};

#endif
