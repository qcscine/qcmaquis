/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef APP_MULTIGRID_MEAS_SIM_H
#define APP_MULTIGRID_MEAS_SIM_H

#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/stat.h>

#include <boost/shared_ptr.hpp>

#include "dmrg/sim/sim.h"
#include "dmrg/models/continuum/factory.h"

#include "dmrg/mp_tensors/optimize.h"
#include "dmrg/mp_tensors/multigrid.h"

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
    using base::measurements;
    using base::stop_callback;
    using base::init_sweep;
    using base::init_site;
    using base::rfile;
    
public:
    mg_meas_sim(DmrgParameters & parms_, ModelParameters & model_)
    : base(parms_, model_)
    , initial_graining(0)
    {
        assert(parms["lattice_library"] == "continuum");
        assert(parms["model_library"] == "continuum");
        
        if (this->restore)
        {
            storage::archive ar(this->chkpfile+"/props.h5");
            ar["/status/graining"] >> initial_graining;
        }
    }
    
    void run()
    {
        /// Set current status in parms
        parms = parms.get_at_index("graining", initial_graining);
        model = model.get_at_index("graining", initial_graining);
        /// Build current model and load/build MPS
        this->model_init();
        this->mps_init();
        
        this->measure("/spectrum/results/", measurements);
        
        double energy = maquis::real(expval(mps, mpoc));
        // MD: removed redundant energy calculation
        // maquis::cout << "Energy before: " << maquis::real(expval(mps, mpo)) << std::endl;
        maquis::cout << "Energy: " << maquis::real(expval(mps, mpoc)) << std::endl;
        {
            storage::archive ar(rfile, "w");
            ar["/spectrum/results/Energy/mean/value"] << std::vector<double>(1, energy);
        }
        
        if (parms["calc_h2"] > 0) {
            MPO<Matrix, SymmGroup> mpo2 = square_mpo(mpoc);
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
