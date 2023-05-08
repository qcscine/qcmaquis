/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#include "maquis_dmrg.h"
#include <complex>

#include "dmrg/sim/symmetry_factory.h"
#include "dmrg/sim/matrix_types.h"
#include "dmrg/sim/interface_sim.h"

namespace maquis
{
    #if defined(HAVE_SU2U1PG)
    typedef SU2U1PG SU2U1grp;
    typedef TwoU1PG TwoU1grp;
    #elif defined(HAVE_SU2U1)
    typedef SU2U1 SU2U1grp;
    typedef TwoU1 TwoU1grp;
    #endif

    template<class ScalarType>
    struct simulation_traits {
        typedef std::shared_ptr<abstract_interface_sim<tmatrix<ScalarType> > > shared_ptr;
        template <class SymmGroup> struct F {
            typedef interface_sim<tmatrix<ScalarType>, SymmGroup> type;
        };
    };

    template <typename ScalarType>
    struct DMRGInterface<ScalarType>::Impl
    {
        typedef typename simulation_traits<ScalarType>::shared_ptr sim_ptr;
        sim_ptr sim;

        Impl(sim_ptr sim_) : sim(sim_) {};
        ~Impl() = default;
    };

    template <typename ScalarType>
    DMRGInterface<ScalarType>::DMRGInterface(DmrgParameters & parms_)
        : parms(parms_), impl_(new Impl(::dmrg::symmetry_factory<simulation_traits<ScalarType> >(parms_, parms_))) {};

    template <typename ScalarType>
    void DMRGInterface<ScalarType>::optimize()
    {
        try {
            impl_->sim->run("optimize");
        }
        catch (std::exception & e) {
            maquis::cerr << "Exception thrown!" << std::endl;
            maquis::cerr << e.what() << std::endl;
            exit(1);
        }
    }

    template <typename ScalarType>
    void DMRGInterface<ScalarType>::evolve()
    {
        try {
            impl_->sim->run("evolve");
        }
        catch (std::exception & e) {
            maquis::cerr << "Exception thrown!" << std::endl;
            maquis::cerr << e.what() << std::endl;
            exit(1);
        }
    }

    template <typename ScalarType>
    void DMRGInterface<ScalarType>::run_measure()
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

    template <typename ScalarType>
    void DMRGInterface<ScalarType>::runInversePowerIteration() {
        try {
            impl_->sim->run("ipi");
        }
        catch (std::exception & e) {
            maquis::cerr << "Exception thrown!" << std::endl;
            maquis::cerr << e.what() << std::endl;
            exit(1);
        }
    }

    template <typename ScalarType>
    void DMRGInterface<ScalarType>::runFEAST() {
        try {
            impl_->sim->run("feast");
        }
        catch (std::exception& e) {
            maquis::cerr << "Exception thrown!" << std::endl;
            maquis::cerr << e.what() << std::endl;
            throw;
            //exit(1);
        }
    }

    template <typename ScalarType>
    ScalarType DMRGInterface<ScalarType>::energy()
    {
        return impl_->sim->get_energy();
    }

    template <typename ScalarType>
    ScalarType DMRGInterface<ScalarType>::energyFEAST(int iState)
    {
        return impl_->sim->getFEASTEnergy(iState);
    }

    template <typename ScalarType>
    ScalarType DMRGInterface<ScalarType>::getCICoefficient(std::string determinantString)
    {
        return impl_->sim->getCICoefficient(determinantString);
    }

    template <typename ScalarType>
    results_collector& DMRGInterface<ScalarType>::get_iteration_results()
    {
        return impl_->sim->get_iteration_results();
    }

    template <typename ScalarType>
    int DMRGInterface<ScalarType>::get_last_sweep()
    {
        return impl_->sim->get_last_sweep();
    }

    template <typename ScalarType>
    void DMRGInterface<ScalarType>::measure()
    {
        measurements_ = impl_->sim->measure_out();
    }

    template <typename ScalarType>
    void DMRGInterface<ScalarType>::update_integrals(const integral_map<ScalarType> & integrals)
    {
        impl_->sim->update_integrals(integrals);
    }

    template <typename ScalarType>
    const typename DMRGInterface<ScalarType>::results_map_type& DMRGInterface<ScalarType>::measurements()
    {
        if (measurements_.empty())
            measure();
        // This is probably not going to work if we call optimize() several times
        // TODO: handle also these cases!
        return measurements_;
    };

    template <typename ScalarType>
    const typename DMRGInterface<ScalarType>::meas_with_results_type& DMRGInterface<ScalarType>::mutinf()
    {
        return measurements().at("mutinf");
    }

    // TODO: This does not work for 2U1/2U1PG symmetry because "oneptdm" measurement is not recognised by the model!
    // Fix the model to recognise it!
    template <typename ScalarType>
    const typename DMRGInterface<ScalarType>::meas_with_results_type& DMRGInterface<ScalarType>::onerdm()
    {
        return measurements().at("oneptdm");
    }

    template <typename ScalarType>
    const typename DMRGInterface<ScalarType>::meas_with_results_type& DMRGInterface<ScalarType>::onespdm()
    {
        return measurements().at("oneptspdm");
    }

    template <typename ScalarType>
    const typename DMRGInterface<ScalarType>::meas_with_results_type& DMRGInterface<ScalarType>::twordm()
    {
        return measurements().at("twoptdm");
    }

    template <typename ScalarType>
    const typename DMRGInterface<ScalarType>::meas_with_results_type& DMRGInterface<ScalarType>::threerdm()
    {
        parms.set("MEASURE[3rdm]", 1); // required for 3-RDM measurement
        return measurements().at("threeptdm");
    }

    template <typename ScalarType>
    const typename DMRGInterface<ScalarType>::meas_with_results_type& DMRGInterface<ScalarType>::fourrdm()
    {
        parms.set("MEASURE[4rdm]", 1); // required for 4-RDM measurement
        return measurements().at("fourptdm");
    }

    template <typename ScalarType>
    const typename DMRGInterface<ScalarType>::meas_with_results_type& DMRGInterface<ScalarType>::getMeasurement(std::string measName)
    {
        if (measurements().find(measName) == measurements().end())
            throw std::runtime_error("Measurement not available!");
        return measurements().at(measName);
    }

    #define measure_and_save_rdm(N) \
        BaseParameters meas_parms = parms.measurements(); \
        parms.erase_measurements(); \
        parms.set("MEASURE[" #N "rdm]", 1); \
        impl_->sim->run_measure(); \
        parms.erase_measurements(); \
        parms << meas_parms

    template <typename ScalarType>
    void DMRGInterface<ScalarType>::measure_and_save_3rdm()
    {
        measure_and_save_rdm(3);
    }

    template <typename ScalarType>
    void DMRGInterface<ScalarType>::measure_and_save_4rdm()
    {
        // Clear all unnecessary measurements before running 4-RDM measurement
        // FIXME: clearing parms here has NO EFFECT on the measurements! This has to be changed in another way!
        // For now the measurements are modified in maquis_cinterface.cpp, but it won't work if DMRGInterface is called directly!
        // Back up measurements
        measure_and_save_rdm(4);
    }

    #undef measure_and_save_rdm

    template <typename ScalarType>
    void DMRGInterface<ScalarType>::measure_and_save_trans3rdm(const std::string & bra_name)
    {
        BaseParameters meas_parms = parms.measurements();
        parms.erase_measurements();
        parms.set("MEASURE[trans3rdm]", bra_name);
        impl_->sim->run_measure();
        parms.erase_measurements();
        parms << meas_parms;
    }

    template <typename ScalarType>
    DMRGInterface<ScalarType>::~DMRGInterface() = default;

    template <typename ScalarType>
    ScalarType DMRGInterface<ScalarType>::overlap(const std::string& aux_mps_name)
    {
        return impl_->sim->get_overlap(aux_mps_name);
    }

    template <typename ScalarType>
    void DMRGInterface<ScalarType>::dump_parameters(const std::string & file)
    {
        std::ofstream fs(file);
        fs << parms;
    }

    // Explicit template instantiation
    template class DMRGInterface<double>;
    template class DMRGInterface<std::complex<double> >;
}
