/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
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

template <class Matrix, class SymmGroup>
class mg_meas_sim : public sim<Matrix, SymmGroup> {
    
    typedef sim<Matrix, SymmGroup> base;
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
        parms << parms.iteration_params("graining", initial_graining);
        model << model.iteration_params("graining", initial_graining);
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
        
        if (parms["MEASURE[EnergyVariance]"] > 0) {
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
