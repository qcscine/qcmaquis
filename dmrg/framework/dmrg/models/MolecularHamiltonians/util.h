/*****************************************************************************
 *
 * QCMaquis DMRG Project
 *
 * Copyright (C) 2015 Laboratory for Physical Chemistry, ETH Zurich
 *               2012-2015 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *               2022- by Alberto Baiardi <abaiardi@ethz.ch>
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

#ifndef QC_CHEM_UTIL_H
#define QC_CHEM_UTIL_H

#include <string>

#include "dmrg/utils/BaseParameters.h"
#include "dmrg/models/model.h"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/utils/align.h"

namespace chem {
namespace detail {

template <class SymmGroup>
struct qn_helper
{
    typename SymmGroup::charge total_qn(BaseParameters & parms)
    {
        typename SymmGroup::charge ret(0);
        ret[0] = parms["u1_total_charge1"];
        ret[1] = parms["u1_total_charge2"];
        return ret;
    }
};

template <>
struct qn_helper<TwoU1PG>
{
    TwoU1PG::charge total_qn(BaseParameters & parms)
    {
        TwoU1PG::charge ret(0);
        ret[0] = parms["u1_total_charge1"];
        ret[1] = parms["u1_total_charge2"];
        ret[2] = parms["irrep"];
        return ret;
    }
};

template <>
struct qn_helper<U1DG>
{
    U1DG::charge total_qn(BaseParameters & parms)
    {
        U1DG::charge ret(0);
        ret[0] = parms["nelec"];
        ret[1] = parms["irrep"];
        return ret;
    }
};

template <>
struct qn_helper<SU2U1>
{
    SU2U1::charge total_qn(BaseParameters & parms)
    {
        SU2U1::charge ret(0);
        ret[0] = parms["nelec"];
        ret[1] = parms["spin"];
        return ret;
    }
};

template <>
struct qn_helper<SU2U1PG>
{
    SU2U1PG::charge total_qn(BaseParameters & parms)
    {
        SU2U1PG::charge ret(0);
        ret[0] = parms["nelec"];
        ret[1] = parms["spin"];
        ret[2] = parms["irrep"];
        return ret;
    }
};

template <class T>
void append(std::vector<T> & target, std::vector<T> const & source) {
    std::copy(source.begin(), source.end(), std::back_inserter(target));
}

template<class T>
typename std::enable_if<boost::is_complex<T>::value,T>::type cconj(T a) {
    return std::conj(a);
}

template<class T>
typename std::enable_if<!boost::is_complex<T>::value,T>::type cconj(T a) {
    return a;
}

// Create 2U1 parameter set from a given L, number of spin-up and spin-down electrons (Nup,Ndown)
// copy over parameters from existing parameters if provided
inline BaseParameters set_2u1_parameters(int L, int Nup, int Ndown,
                                         const BaseParameters& existing_pars = BaseParameters())
{
    BaseParameters ret(existing_pars);
    ret.set("lattice_library", "coded");
    ret.set("LATTICE", "orbitals");
    ret.set("model_library", "coded");
    ret.set("MODEL", "quantum_chemistry");
    ret.set("L", L);
    ret.set("u1_total_charge1", Nup);
    ret.set("u1_total_charge2", Ndown);
    return ret;
}

template <class Matrix, class SymmGroup>
inline std::string infer_site_types(MPS<Matrix, SymmGroup> const & mps)
{
    // determine the irreps per site
    std::string site_types;
    for (Lattice::pos_t p = 0; p < mps.size(); ++p)
        for (std::size_t i = 0; i < mps[p].site_dim().size(); ++i)
        {
            if (SymmGroup::particleNumber(mps[p].site_dim()[i].first) % 2 != 0)
            {
                site_types += boost::lexical_cast<std::string>(getPG<SymmGroup>()(mps[p].site_dim()[i].first)) + ",";
                break;
            }
            if (i == mps[p].site_dim().size() -1)
                site_types += "0,";
        }
    return site_types;
}

template <class SymmGroup>
inline typename SymmGroup::charge make_2u1_initc(int Nup, int Ndown, int irrep)
{
    typename SymmGroup::charge ret;
    ret[0] = Nup;
    ret[1] = Ndown;
    ret = PGCharge<SymmGroup>()(ret, irrep);
    return ret;
}

template <class Matrix, class SymmGroup>
inline std::vector<Index<SymmGroup> > make_2u1_site_basis(int L, int Nup, int Ndown,
                                                          std::string site_types)
{
    BaseParameters parms = set_2u1_parameters(L, Nup, Ndown);
    parms.set("site_types", site_types);
    Lattice lat(parms);
    Model<Matrix, SymmGroup> model(lat, parms);
    std::vector<Index<SymmGroup> > site_bases;
    for (int i = 0; i <= lat.maximum_vertex_type(); ++i)
        site_bases.push_back(model.phys_dim(i));
    return site_bases;
}

} // chem
} // detail

#endif
