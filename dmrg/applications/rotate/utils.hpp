/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2016 Laboratory of Physical Chemistry, ETH Zurich
 *               2016 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

#ifndef MPS_ROTATE_UTILS_HPP
#define MPS_ROTATE_UTILS_HPP

namespace debug
{
    template <class Matrix, class SymmGroup>
    void mps_print(MPS<Matrix, SymmGroup> const & mps, std::string title)
    {
        for (Lattice::pos_t p = 0; p < mps.length(); ++p)
        {
            maquis::cout << title << p << ":" << std::endl;        
            maquis::cout << mps[p];
        }
    }
    
    template <class Matrix, class SymmGroup>
    void mps_print(MPSTensor<Matrix, SymmGroup> const & mps, std::string title)
    {
        maquis::cout << title << ":" << std::endl;        
        maquis::cout << mps;
    }

    template <class Matrix>
    void mps_print_ci_(MPS<Matrix, TwoU1PG> const & mps, std::string detfile)
    {
        typedef Lattice::pos_t pos_t;
        typedef typename TwoU1PG::charge charge;

        // extract physical basis for every site from MPS
        pos_t L = mps.length();
        TwoU1PG::subcharge Nup = mps[L-1].col_dim()[0].first[0];
        TwoU1PG::subcharge Ndown = mps[L-1].col_dim()[0].first[1];
        std::string site_types = chem_detail::infer_site_types(mps);

        // extract physical basis for every site from MPS
        BaseParameters parms;
        parms.set("site_types", site_types);
        std::vector<TwoU1PG::subcharge> irreps = parms["site_types"];
        std::vector<Index<TwoU1PG> > per_site, phys_dims = chem_detail::make_2u1_site_basis<Matrix, TwoU1PG>(L, Nup, Ndown, site_types);
        for (pos_t q = 0; q < L; ++q)
            per_site.push_back(phys_dims[irreps[q]]);

        // load the determinants
        std::vector<std::vector<charge> > determinants = parse_config<Matrix, TwoU1PG>(detfile, per_site);
        // printout the determinants
        for (pos_t q = 0; q < determinants.size(); ++q){
           for (pos_t p = 0; p < L; ++p)
               std::cout << determinants[q][p];

           std::cout << std::endl;
        }

        // compute the CI coefficients for all determinants in the input
        for (typename std::vector< std::vector<charge> >::iterator it = determinants.begin(); it != determinants.end(); ++it)
            maquis::cout << "CI coefficient of det " << it - determinants.begin() + 1 << ": " << extract_coefficient(mps, *it) << std::endl;

        maquis::cout << std::endl;
    }

    template <class Matrix, class SymmGroup>
    void mps_print_ci(MPS<Matrix, SymmGroup> const & mps_in, std::string detfile)
    {
    }

    template <>
    void mps_print_ci(MPS<alps::numeric::matrix<double>, TwoU1PG> const & mps_in, std::string detfile)
    {
        mps_print_ci_(mps_in, detfile);
    }

    template <>
    void mps_print_ci(MPS<alps::numeric::matrix<double>, SU2U1PG> const & mps_in, std::string detfile)
    {
        typedef alps::numeric::matrix<double> Matrix;

        // the output MPS
        size_t L = mps_in.size();
        size_t N = SU2U1PG::particleNumber(mps_in[L-1].col_dim()[0].first);
        size_t S = SU2U1PG::spin(mps_in[L-1].col_dim()[0].first);
        MPS<Matrix, TwoU1PG> mps = transform_mps<Matrix, SU2U1PG>()(mps_in, (N+S)/2, (N-S)/2);

        mps_print_ci_(mps, detfile);
    }
} // namespace debug

#endif
