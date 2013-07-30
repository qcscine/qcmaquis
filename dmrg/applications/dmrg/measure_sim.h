/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2013 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef APP_DMRG_MEASURE_SIM_H
#define APP_DMRG_MEASURE_SIM_H

#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/stat.h>

#include <boost/shared_ptr.hpp>

#include "dmrg/sim/sim.h"
#include "dmrg/mp_tensors/optimize.h"


template <class Matrix, class SymmGroup>
class measure_sim : public sim<Matrix, SymmGroup> {
    
    typedef sim<Matrix, SymmGroup> base;
    typedef optimizer_base<Matrix, SymmGroup, storage::disk> opt_base_t;
    
    using base::mps;
    using base::mpo;
    using base::mpoc;
    using base::parms;
    using base::model;
    using base::measurements;
    using base::stop_callback;
    using base::rfile;
    
public:
    
    measure_sim (DmrgParameters & parms_, ModelParameters & model_)
    : base(parms_, model_)
    { }
    
    void run()
    {
        this->measure("/spectrum/results", measurements);
        
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
};

#endif
