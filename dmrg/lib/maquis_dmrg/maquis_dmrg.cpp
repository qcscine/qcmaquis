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

#include "dmrg/sim/symmetry_factory.h"
#include "dmrg/sim/matrix_types.h"
#include "interface_sim.h"

namespace maquis
{
    template<class V>
    struct simulation_traits {
        typedef std::shared_ptr<abstract_interface_sim<tmatrix<V> > > shared_ptr;
        template <class SymmGroup> struct F {
            typedef interface_sim<tmatrix<V>, SymmGroup> type;
        };
    };

    template <class V>
    struct DMRGInterface<V>::Impl
    {
        typedef typename simulation_traits<V>::shared_ptr sim_ptr;
        sim_ptr sim;

        Impl(sim_ptr sim_) : sim(sim_) {};
        ~Impl() = default;
    };

    template <class V>
    DMRGInterface<V>::DMRGInterface(DmrgParameters & parms_)
        : parms(parms_), impl_(new Impl(::dmrg::symmetry_factory<simulation_traits<V> >(parms_, parms_))) {};

    template <class V>
    void DMRGInterface<V>::optimize()
    {
        try
        {
            impl_->sim->run();
        }
        catch (std::exception & e)
        {
            maquis::cerr << "Exception thrown!" << std::endl;
            maquis::cerr << e.what() << std::endl;
            exit(1);
        }
    }

    template <class V>
    void DMRGInterface<V>::run_measure()
    {
        try
        {
            impl_->sim->run_measure();
        }
        catch (std::exception & e)
        {
            maquis::cerr << "Exception thrown!" << std::endl;
            maquis::cerr << e.what() << std::endl;
            exit(1);
        }
    }

    template <class V>
    V DMRGInterface<V>::energy()
    {
        return impl_->sim->get_energy();
    }

    template <class V>
    results_collector& DMRGInterface<V>::get_iteration_results()
    {
        return impl_->sim->get_iteration_results();
    }

    template <class V>
    int DMRGInterface<V>::get_last_sweep()
    {
        return impl_->sim->get_last_sweep();
    }

    template <class V>
    void DMRGInterface<V>::measure()
    {
        measurements_ = impl_->sim->measure_out();
    }

    template <class V>
    void DMRGInterface<V>::update_integrals(const integral_map<V> & integrals)
    {
        impl_->sim->update_integrals(integrals);
    }

    template <class V>
    const typename DMRGInterface<V>::results_map_type& DMRGInterface<V>::measurements()
    {
        if (measurements_.empty())
            measure();
        // This is probably not going to work if we call optimize() several times
        // TODO: handle also these cases!
        return measurements_;
    };

    // TODO: This does not work for 2U1/2U1PG symmetry because "oneptdm" measurement is not recognised by the model!
    // Fix the model to recognise it!
    template <class V>
    const typename DMRGInterface<V>::meas_with_results_type& DMRGInterface<V>::onerdm()
    {
        return measurements().at("oneptdm");
    }

    template <class V>
    const typename DMRGInterface<V>::meas_with_results_type& DMRGInterface<V>::twordm()
    {
        return measurements().at("twoptdm");
    }

    template <class V>
    const typename DMRGInterface<V>::meas_with_results_type& DMRGInterface<V>::fourrdm()
    {
        parms.set("MEASURE[4rdm]", 1); // required for 4-RDM measurement
        return measurements().at("fourptdm");
    }

    template <class V>
    DMRGInterface<V>::~DMRGInterface() = default;

    template <class V>
    V DMRGInterface<V>::overlap(const std::string& aux_mps_name)
    {
        return impl_->sim->get_overlap(aux_mps_name);
    }

    template class DMRGInterface<double>;
    template class DMRGInterface<std::complex<double> >;
}