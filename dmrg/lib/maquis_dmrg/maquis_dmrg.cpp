/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
*               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
*               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
*               2014-2014    Sebastian Keller <sebkelle@phys.ethz.ch>
*               2018-2019    Leon Freitag <lefreita@ethz.ch>
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

#include "maquis_dmrg.h"
#include <complex>

namespace maquis
{
    template <class V, bool Relativistic>
    DMRGInterface<V,Relativistic>::DMRGInterface(DmrgParameters & parms_)
        : parms(parms_), sim(::dmrg::symmetry_factory<simulation_traits<V> >(parms_, parms_)) {};

    template <class V, bool Relativistic>
    void DMRGInterface<V,Relativistic>::optimize()
    {
        try
        {
            sim->run();
        }
        catch (std::exception & e)
        {
            maquis::cerr << "Exception thrown!" << std::endl;
            maquis::cerr << e.what() << std::endl;
            exit(1);
        }
    }

    template <class V, bool Relativistic>
    void DMRGInterface<V,Relativistic>::run_measure()
    {
        try
        {
            sim->run_measure();
        }
        catch (std::exception & e)
        {
            maquis::cerr << "Exception thrown!" << std::endl;
            maquis::cerr << e.what() << std::endl;
            exit(1);
        }
    }

    template <class V, bool Relativistic>
    V DMRGInterface<V,Relativistic>::energy()
    {
        return sim->get_energy();
    }

    template <class V, bool Relativistic>
    results_collector& DMRGInterface<V,Relativistic>::get_iteration_results()
    {
        return sim->get_iteration_results();
    }

    template <class V, bool Relativistic>
    int DMRGInterface<V,Relativistic>::get_last_sweep()
    {
        return sim->get_last_sweep();
    }

    template <class V, bool Relativistic>
    void DMRGInterface<V,Relativistic>::measure()
    {
        measurements_ = sim->measure_out();
    }

    template <class V, bool Relativistic>
    void DMRGInterface<V,Relativistic>::update_integrals(const integral_map<V,Relativistic> & integrals)
    {
        sim->update_integrals(integrals);
    }

    template <class V, bool Relativistic>
    const typename DMRGInterface<V,Relativistic>::results_map_type& DMRGInterface<V,Relativistic>::measurements()
    {
        if (measurements_.empty())
            measure();
        return measurements_;
    };

    // TODO: This does not work for 2U1/2U1PG symmetry because "oneptdm" measurement is not recognised by the model!
    // Fix the model to recognise it!
    template <class V, bool Relativistic>
    const typename DMRGInterface<V,Relativistic>::meas_with_results_type& DMRGInterface<V,Relativistic>::onerdm()
    {
        return measurements().at("oneptdm");
    }

    template <class V, bool Relativistic>
    const typename DMRGInterface<V,Relativistic>::meas_with_results_type& DMRGInterface<V,Relativistic>::twordm()
    {
        return measurements().at("twoptdm");
    }

    //template class DMRGInterface<double, static const bool Relativistic = false>;
   // template class DMRGInterface<std::complex<double>, static const bool Relativistic = false>;
}
