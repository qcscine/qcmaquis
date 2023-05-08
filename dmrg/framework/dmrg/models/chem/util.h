/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

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

    class IndexTuple : public NU1Charge<4>
    {
    public:
        typedef std::array<int, 4> array_type;
        typedef NU1Charge<4> base;
        using base::base;

        IndexTuple() {}
        IndexTuple(int i, int j, int k, int l) {
            (*this)[0] = i; (*this)[1] = j; (*this)[2] = k; (*this)[3] = l;
        }
    };

	template <class SymmGroup=TrivialGroup>
    inline IndexTuple align(int i, int j, int k, int l) {
        return IndexTuple(maquis::detail::align<false>({i,j,k,l}));
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
    typename std::enable_if<boost::is_complex<T>::value,T>::type cconj(T a)
    {
        return std::conj(a);
    }
    template<class T>
    typename std::enable_if<!boost::is_complex<T>::value,T>::type cconj(T a)
    {
        return a;
    }

    // Create 2U1 parameter set from a given L, number of spin-up and spin-down electrons (Nup,Ndown)
    // copy over parameters from existing parameters if provided
    inline
    BaseParameters set_2u1_parameters(int L, int Nup, int Ndown, const BaseParameters& existing_pars = BaseParameters())
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
    inline std::vector<Index<SymmGroup> > make_2u1_site_basis(int L, int Nup, int Ndown, std::string site_types)
    {
        BaseParameters parms = set_2u1_parameters(L, Nup, Ndown);
        parms.set("site_types", site_types);
        Lattice lat(parms);
        Model<Matrix, SymmGroup> model(lat, parms);
        std::vector<Index<SymmGroup> > site_bases;
        for (int iType = 0; iType < lat.getMaxType(); iType++)
            site_bases.push_back(model.phys_dim(iType));
        return site_bases;
    }
}
}

#endif
