/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef APP_MULTIGRID_MEAS_SIM_H
#define APP_MULTIGRID_MEAS_SIM_H

#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/stat.h>

#include <boost/shared_ptr.hpp>

#include "dmrg/sim/sim.h"
#include "dmrg/models/continuum/factory.h"

#include "dmrg/optimize/optimize.h"
#include "dmrg/mp_tensors/multigrid.h"

inline BaseParameters compute_initial_parms(BaseParameters parms)
{
    int initial_graining = 0;
    
    std::string chkpfile = boost::trim_right_copy_if(parms["chkpfile"].str(), boost::is_any_of("/ "));
    boost::filesystem::path p(chkpfile);
    if (boost::filesystem::exists(p) && boost::filesystem::exists(p / "mps0.h5")) {
        storage::archive ar(chkpfile+"/props.h5");
        if (ar.is_data("/status/graining") && ar.is_scalar("/status/graining"))
            ar["/status/graining"] >> initial_graining;
    }
    
    parms << parms.iteration_params("graining", initial_graining);
    return parms;
}


template <class Matrix, class SymmGroup>
class mg_meas_sim : public sim<Matrix, SymmGroup> {
    
    typedef sim<Matrix, SymmGroup> base;
    typedef optimizer_base<Matrix, SymmGroup, storage::disk> opt_base_t;
    typedef typename base::status_type status_type;

    enum measure_t {sweep_measure, mg_measure};
    
    using base::mps;
    using base::mpo;
    using base::lat;
    using base::mpoc;
    using base::parms;
    using base::model;
    using base::all_measurements;
    using base::stop_callback;
    using base::init_sweep;
    using base::init_site;
    using base::rfile;
    
public:
    mg_meas_sim(DmrgParameters & parms_)
    : base(compute_initial_parms(parms_))
    , initial_graining(0)
    {
        if (this->restore)
        {
            storage::archive ar(this->chkpfile+"/props.h5");
            ar["/status/graining"] >> initial_graining;
        }
    }
    
    void model_init()
    {
        /// Model initialization
        this->lat = Lattice(this->parms);
        this->model = Model<Matrix, SymmGroup>(this->lat, this->parms);
        this->mpo = make_mpo(this->lat, this->model);
        this->all_measurements = this->model.measurements();
        this->all_measurements << overlap_measurements<Matrix, SymmGroup>(this->parms);
    }

    void run()
    {
        /// Set current status in parms
        parms << parms.iteration_params("graining", initial_graining);
        /// Build current model and load/build MPS
        this->model_init();
        
        this->measure("/spectrum/results/", all_measurements);
        
        double energy = maquis::real(expval(mps, mpo));
        // MD: removed redundant energy calculation
        // maquis::cout << "Energy before: " << maquis::real(expval(mps, mpo)) << std::endl;
        maquis::cout << "Energy: " << maquis::real(expval(mps, mpo)) << std::endl;
        {
            storage::archive ar(rfile, "w");
            ar["/spectrum/results/Energy/mean/value"] << std::vector<double>(1, energy);
        }
        
        if (parms["MEASURE[EnergyVariance]"] > 0) {
            MPO<Matrix, SymmGroup> mpo2 = square_mpo(mpo);
            mpo2.compress(1e-12);
            
            double energy2 = maquis::real(expval(mps, mpo2, true));
            
            maquis::cout << "Energy^2: " << energy2 << std::endl;
            maquis::cout << "Variance: " << energy2 - energy*energy << std::endl;
            
            {
                storage::archive ar(rfile, "w");
                ar["/spectrum/results/Energy^2/mean/value"] << std::vector<double>(1, energy2);
                ar["/spectrum/results/EnergyVariance/mean/value"] << std::vector<double>(1, energy2 - energy*energy);
            }
        }
    }
    
private:
    int initial_graining;
};

#endif
