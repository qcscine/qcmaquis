/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

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

    static subcharge adjoin(subcharge I)
    {
        return I;
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
