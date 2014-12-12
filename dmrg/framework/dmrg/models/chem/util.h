/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2012-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

namespace chem_detail {

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
            ret[2] = parms["irrep_charge"];
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

    inline IndexTuple align(int i, int j, int k, int l) {
        if (i<j) std::swap(i,j);
        if (k<l) std::swap(k,l);
        if (i<k) { std::swap(i,k); std::swap(j,l); }
        if (i==k && j<l) { std::swap(j,l); }
        return IndexTuple(i,j,k,l);
    }
    
    inline IndexTuple align(IndexTuple const & rhs) {
        return align(rhs[0], rhs[1], rhs[2], rhs[3]);
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

    class TermTuple : public NU1Charge<8>
    {
    public:
        TermTuple() {}
        TermTuple(IndexTuple const & a, IndexTuple const & b) {
            for (int i=0; i<4; i++) { (*this)[i] = a[i]; (*this)[i+4] = b[i]; }
        }
    };
}

#endif
