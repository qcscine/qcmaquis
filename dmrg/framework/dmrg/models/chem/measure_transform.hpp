/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2015 Laboratory of Physical Chemistry, ETH Zurich
 *               2015-2015 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *               2020 by Leon Freitag <lefreita@ethz.ch>
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


template <class Matrix, class SymmGroup, class = void>
struct measure_transform
{
    typedef typename Model<Matrix, SymmGroup>::results_map_type results_map_type;
    typedef typename boost::mpl::if_<symm_traits::HasPG<SymmGroup>, TwoU1PG, TwoU1>::type SymmOut;

    void operator()(std::string rfile, std::string result_path, Lattice lat, MPS<Matrix, SymmGroup> const & mps,
        BaseParameters const & measurement_parms = BaseParameters()) {}

    results_map_type meas_out(Lattice lat, MPS<Matrix, SymmGroup> const & mps,
        BaseParameters const & measurement_parms = BaseParameters(), const std::string & rfile="",
        const std::string & result_path="") { return results_map_type(); };
};

template <class Matrix, class SymmGroup>
struct measure_transform<Matrix, SymmGroup, typename boost::enable_if<symm_traits::HasSU2<SymmGroup> >::type>
{
    typedef typename Model<Matrix, SymmGroup>::results_map_type results_map_type;
    typedef typename boost::mpl::if_<symm_traits::HasPG<SymmGroup>, TwoU1PG, TwoU1>::type SymmOut;

    // Measure and output into file rfile
    void operator()(std::string rfile, std::string result_path, Lattice lat, MPS<Matrix, SymmGroup> const & mps,
    BaseParameters const & measurement_parms = BaseParameters() )
    {

        MPS<Matrix, SymmOut> mps_tmp;
        typename Model<Matrix, SymmOut>::measurements_type transformed_measurements;

        // --- disabled for now ---
        // For compatibility with older versions we need MEASURE[ChemEntropy]
        // in transformed measurements here -- at least I guess so
        // if (measurement_parms.empty())
        // {
        //     BaseParameters parms_chementropy;
        //     parms_chementropy.set("MEASURE[ChemEntropy]", 1);
        //     std::tie(mps_tmp, transformed_measurements) = prepare_measurements(lat, mps, parms_chementropy);
        // }
        // else
        {
            std::tie(mps_tmp, transformed_measurements) = prepare_measurements(lat, mps, measurement_parms);
        }

        std::for_each(transformed_measurements.begin(), transformed_measurements.end(),
                      measure_and_save<Matrix, SymmOut>(rfile, result_path, mps_tmp));
    };

    // Measure and return a map with results
    results_map_type meas_out(Lattice lat, MPS<Matrix, SymmGroup> const & mps,
        BaseParameters const & measurement_parms = BaseParameters(), const std::string & rfile="", const std::string & result_path="")
    {
        results_map_type ret;
        MPS<Matrix, SymmOut> mps_tmp;
        typename Model<Matrix, SymmOut>::measurements_type transformed_measurements;

        // Unlike above, no need for setting ChemEntropy here
        std::tie(mps_tmp, transformed_measurements) = prepare_measurements(lat, mps, measurement_parms);

        for(auto&& meas: transformed_measurements)
            ret[meas.name()] = measure_and_save<Matrix,SymmOut>(rfile, result_path, mps_tmp).meas_out(meas);

        return ret;
    }

    private:
        // Prepare measurements and return the transformed MPS and the measurements

        std::pair<MPS<Matrix, SymmOut>, typename Model<Matrix, SymmOut>::measurements_type>
        prepare_measurements(Lattice lat, MPS<Matrix, SymmGroup> const & mps,
         const BaseParameters& measurement_parms = BaseParameters())
        {
            typedef typename boost::mpl::if_<symm_traits::HasPG<SymmGroup>, TwoU1PG, TwoU1>::type SymmOut;

            int N = SymmGroup::particleNumber(mps[mps.size()-1].col_dim()[0].first);
            int TwoS = SymmGroup::spin(mps[mps.size()-1].col_dim()[0].first);
            int Nup = (N + TwoS) / 2;
            int Ndown = (N - TwoS) / 2;

            BaseParameters parms_tmp = chem::detail::set_2u1_parameters(mps.size(), Nup, Ndown, measurement_parms);
            //parms_tmp.set("MEASURE[ChemEntropy]", 1);

            Model<Matrix, SymmOut> model_tmp(lat, parms_tmp);

            return std::make_pair(transform_mps<Matrix, SymmGroup>()(mps, Nup, Ndown), model_tmp.measurements());
        }

    // transformed measurements for MS-MPS
    // Not ready yet
    /*
    void operator()(std::vector<std::string> rfiles, std::string result_path, Lattice lat, const std::vector<MPS<Matrix,SymmGroup > > & mps_vec)
    {
        typedef typename boost::mpl::if_<symm_traits::HasPG<SymmGroup>, TwoU1PG, TwoU1>::type SymmOut;

        assert(mps_vec.length() > 0);
        const MPS<Matrix,SymmGroup> & mps = mps_vec[0];

        int N = SymmGroup::particleNumber(mps[mps.size()-1].col_dim()[0].first);
        int TwoS = SymmGroup::spin(mps[mps.size()-1].col_dim()[0].first);
        int Nup = (N + TwoS) / 2;
        int Ndown = (N - TwoS) / 2;

        BaseParameters parms_tmp = chem::detail::set_2u1_parameters(mps.size(), Nup, Ndown);
        parms_tmp.set("MEASURE[ChemEntropy]", 1);

        Model<Matrix, SymmOut> model_tmp(lat, parms_tmp);

        std::vector<MPS<Matrix,SymmGroup> > mps_out;

        std::transform(mps_vec.begin(), mps_vec.end(), mps_out.begin(), [&](const MPS<Matrix,SymmGroup> & m) { return transform_mps<Matrix, SymmGroup>()(m, Nup, Ndown); });
        // now put the equivalent of measure_and_save for multi-state MPS here
        ...
    }
    */
};

#endif
