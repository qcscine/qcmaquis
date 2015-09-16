/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2015 Laboratory of Physical Chemistry, ETH Zurich
 *               2015-2015 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

#ifndef MEASURE_TRANSFORM_HPP
#define MEASURE_TRANSFORM_HPP

#include <boost/mpl/if.hpp>

#include "dmrg/models/chem/transform_symmetry.hpp"
#include "dmrg/models/measurements.h"

namespace transform_detail {

    template <class SymmGroup, class = void>
    struct insert_symmetry
    {
        void operator()(DmrgParameters parms, DmrgParameters & parms_out)
        {
            parms_out.set("symmetry", "2u1");
        }
    };

    template <class SymmGroup>
    struct insert_symmetry<SymmGroup, typename boost::enable_if<symm_traits::HasPG<SymmGroup> >::type>
    {
        void operator()(DmrgParameters parms, DmrgParameters & parms_out)
        {
            parms_out.set("symmetry", "2u1pg");
            parms_out.set("irrep", parms["irrep"]);
        }
    };
}

template <class Matrix, class SymmGroup, class = void>
struct measure_transform
{
    void operator()(std::string rfile, std::string result_path, DmrgParameters parms, Lattice lat, MPS<Matrix, SymmGroup> const & mps)
    {}
};

template <class Matrix, class SymmGroup>
struct measure_transform<Matrix, SymmGroup, typename boost::enable_if<symm_traits::HasSU2<SymmGroup> >::type>
{
    void operator()(std::string rfile, std::string result_path, DmrgParameters parms, Lattice lat, MPS<Matrix, SymmGroup> const & mps)
    {
        typedef typename boost::mpl::if_<symm_traits::HasPG<SymmGroup>, TwoU1PG, TwoU1>::type SymmOut;

        int N = parms["nelec"];
        int TwoS = parms["spin"];
        int Nup = (N + TwoS) / 2;
        int Ndown = (N - TwoS) / 2;

        DmrgParameters parms_tmp;
        parms_tmp.set("model_library", "coded");
        parms_tmp.set("MODEL", "quantum_chemistry");
        parms_tmp.set("integral_file", parms["integral_file"]);

        parms_tmp.set("u1_total_charge1", Nup);
        parms_tmp.set("u1_total_charge2", Ndown);
        parms_tmp.set("L", parms["L"]);
        parms_tmp.set("init_state", "const");
        parms_tmp.set("MEASURE[ChemEntropy]", 1);

        transform_detail::insert_symmetry<SymmGroup>()(parms, parms_tmp);

        Model<Matrix, SymmOut> model_tmp(lat, parms);

        MPS<Matrix, SymmOut> mps_tmp(lat.size(), *(model_tmp.initializer(lat, parms_tmp)));
        transform_mps<Matrix, SymmGroup>()(mps, mps_tmp);

        typename Model<Matrix, SymmOut>::measurements_type transformed_measurements = model_tmp.measurements();
        std::for_each(transformed_measurements.begin(), transformed_measurements.end(),
                      measure_and_save<Matrix, SymmOut>(rfile, result_path, mps_tmp));
    };
};

#endif
