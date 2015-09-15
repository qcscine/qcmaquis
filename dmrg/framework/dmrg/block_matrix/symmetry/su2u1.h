/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2014-2014 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

#ifndef SYMMETRY_SU2U1_H
#define SYMMETRY_SU2U1_H

#include "dmrg/block_matrix/symmetry/nu1_tpl.h"
#include "dmrg/block_matrix/symmetry/nu1pg.h"

// SU2 x U1 Symmetry
// ---
// reusing TwoU1 charges

template<class S = int>
class SU2U1_template
{
public:
    typedef S subcharge;
    typedef NU1Charge<2, S> charge;
    typedef std::vector<charge> charge_v;

    static const charge IdentityCharge;
    static const bool finite = false;

    static subcharge particleNumber(charge rhs) { return rhs[0]; }

    static subcharge & spin(charge & rhs) { return rhs[1]; }
    static subcharge const & spin(charge const & rhs) { return rhs[1]; }

    static charge fuse(charge a, charge b)
    {
        return a+b;
    }

    template<int R> static charge fuse(boost::array<charge, R> const & v)
    {
        charge ret = v[0];
        for (int i = 1; i < R; ++i)
            ret = fuse(ret, v[i]);
        return ret;
    }
};

template<class S> const typename SU2U1_template<S>::charge SU2U1_template<S>::IdentityCharge = typename SU2U1_template<S>::charge();

typedef SU2U1_template<> SU2U1;


// SU2 x U1 x PG Symmetry
// ---
// reusing TwoU1PG charges

template<class S = int>
class SU2U1PG_template
{
public:
    typedef S subcharge;
    typedef NU1ChargePG<2, S> charge;
    typedef std::vector<charge> charge_v;

    static const charge IdentityCharge;
    static const bool finite = false;

    static subcharge particleNumber(charge rhs) { return rhs[0]; }

    static subcharge & spin(charge & rhs) { return rhs[1]; }
    static subcharge const & spin(charge const & rhs) { return rhs[1]; }

    static subcharge & irrep(charge & rhs) { return rhs[2]; }
    static subcharge const & irrep(charge const & rhs) { return rhs[2]; }

    static charge fuse(charge a, charge b)
    {
        return a+b;
    }

    template<int R> static charge fuse(boost::array<charge, R> const & v)
    {
        charge ret = v[0];
        for (int i = 1; i < R; ++i)
            ret = fuse(ret, v[i]);
        return ret;
    }
};

template<class S> const typename SU2U1PG_template<S>::charge SU2U1PG_template<S>::IdentityCharge = typename SU2U1PG_template<S>::charge();

typedef SU2U1PG_template<> SU2U1PG;

#endif
