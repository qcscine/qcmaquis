/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef APP_DMRG_MEASURE_SIM_H
#define APP_DMRG_MEASURE_SIM_H

#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/stat.h>

#include <boost/shared_ptr.hpp>

#include "dmrg/sim/sim.h"
#include "dmrg/optimize/optimize.h"


template <class Matrix, class SymmGroup>
class measure_sim : public sim<Matrix, SymmGroup> {
    
    typedef sim<Matrix, SymmGroup> base;
    typedef optimizer_base<Matrix, SymmGroup, storage::disk> opt_base_t;
    
    using base::mps;
    using base::mpo;
    using base::parms;
    using base::all_measurements;
    using base::stop_callback;
    using base::rfile;
    
public:
    
    measure_sim (DmrgParameters & parms_)
    : base(parms_)
    { }
    
    void run()
    {
        this->measure("/spectrum/results/", all_measurements);
        
        /// MPO creation
        MPO<Matrix, SymmGroup> mpoc = mpo;
        if (parms["use_compressed"])
            mpoc.compress(1e-12);

        double energy = maquis::real(expval(mps, mpoc));

        if (parms["MEASURE[Energy]"]) {
            maquis::cout << "Energy: " << maquis::real(expval(mps, mpoc)) << std::endl;
            {
                storage::archive ar(rfile, "w");
                ar["/spectrum/results/Energy/mean/value"] << std::vector<double>(1, energy);
            }
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
};

#endif
