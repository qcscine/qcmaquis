/*****************************************************************************
 *
 * QCMaquis DMRG Project
 *
 * Copyright (C) 2015 Laboratory for Physical Chemistry, ETH Zurich
 *               2012-2015 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *
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
#include "dmrg/models/lattice.h"

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

    class IndexTuple : public NU1Charge<4>
    {
    public:
        IndexTuple() {}
        IndexTuple(int i, int j, int k, int l) {
            (*this)[0] = i; (*this)[1] = j; (*this)[2] = k; (*this)[3] = l;
        }
    };

	template <class SymmGroup>
    inline IndexTuple align(int i, int j, int k, int l) {
        if (i<j) std::swap(i,j);
        if (k<l) std::swap(k,l);
        if (i<k) { std::swap(i,k); std::swap(j,l); }
        if (i==k && j<l) { std::swap(j,l); }
        return IndexTuple(i,j,k,l);
    }

	template <>
    inline IndexTuple align<U1DG>(int i, int j, int k, int l) {
        return IndexTuple(i,j,k,l);
    }

	template <class SymmGroup>
    inline IndexTuple align(IndexTuple const & rhs) {
        return align<SymmGroup>(rhs[0], rhs[1], rhs[2], rhs[3]);
    }

    inline int sign(IndexTuple const & idx)
    {
        int inv_count=0, n=4;
        for(int c1 = 0; c1 < n - 1; c1++)
            for(int c2 = c1+1; c2 < n; c2++)
                if(idx[c1] > idx[c2]) inv_count++;

        return 1 - 2 * (inv_count % 2);
    }

    inline std::ostream& operator<<(std::ostream & os, IndexTuple const & c) {
        os << "<";
        for (int i = 0; i < 4; ++i) {
            os << c[i];
            if (i+1 < 4)
                os << ",";
        }
        os << ">";
        return os;
    }

    class SixTuple : public NU1Charge<6>
    {
    public:
        SixTuple() {}
        SixTuple(int t0, int t1, int t2, int t3, int t4, int t5) {
            (*this)[0] = t0;
            (*this)[1] = t1;
            (*this)[2] = t2;
            (*this)[3] = t3;
            (*this)[4] = t4;
            (*this)[5] = t5;
        }
    };

    class EightTuple : public NU1Charge<8>
    {
    public:
        EightTuple() {}
        EightTuple(IndexTuple const & a, IndexTuple const & b) {
            for (int i=0; i<4; i++) { (*this)[i] = a[i]; (*this)[i+4] = b[i]; }
        }
    };

    template <class T>
    void append(std::vector<T> & target, std::vector<T> const & source)
    {
        std::copy(source.begin(), source.end(), std::back_inserter(target));
    }

    template<class T>
    typename boost::enable_if<boost::is_complex<T>,T>::type cconj(T a)
    {
        return std::conj(a);
    }
    template<class T>
    typename boost::disable_if<boost::is_complex<T>,T>::type cconj(T a)
    {
        return a;
    }

    inline
    BaseParameters set_2u1_parameters(int L, int Nup, int Ndown)
    {
        BaseParameters ret;

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
    inline
    std::string infer_site_types(MPS<Matrix, SymmGroup> const & mps)
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
    inline
    typename SymmGroup::charge make_2u1_initc(int Nup, int Ndown, int irrep)
    {
        typename SymmGroup::charge ret;
        ret[0] = Nup;
        ret[1] = Ndown;
        ret = PGCharge<SymmGroup>()(ret, irrep);

        return ret;
    }

    template <class Matrix, class SymmGroup>
    inline
    std::vector<Index<SymmGroup> > make_2u1_site_basis(int L, int Nup, int Ndown, std::string site_types)
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
}
}

#endif
