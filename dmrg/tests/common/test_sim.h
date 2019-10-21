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

#ifndef APP_DMRG_SIM_H
#define APP_DMRG_SIM_H

#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/stat.h>

#include <boost/shared_ptr.hpp>

#include "dmrg/sim/sim.h"
#include "dmrg/optimize/optimize.h"


template <class Matrix, class SymmGroup>
class test_sim : public sim<Matrix, SymmGroup> {

    typedef sim<Matrix, SymmGroup> base;
    typedef optimizer_base<Matrix, SymmGroup, storage::disk> opt_base_t;
    typedef typename base::status_type status_type;
    typedef typename base::measurements_type measurements_type;

    using base::mps;
    using base::mpo;
    using base::parms;
    using base::all_measurements;
    using base::stop_callback;
    using base::init_sweep;
    using base::init_site;
    using base::base;

public:

    void run()
    {

        //
        // Optimizer initialization
        // ------------------------
        std::shared_ptr<opt_base_t> optimizer;
        if (parms["optimization"] == "singlesite")
        {
            optimizer.reset( new ss_optimize<Matrix, SymmGroup, storage::disk>
                            (mps, mpo, parms, stop_callback, init_site) );
        }
        else if(parms["optimization"] == "twosite")
        {
            optimizer.reset( new ts_optimize<Matrix, SymmGroup, storage::disk>
                            (mps, mpo, parms, stop_callback, init_site) );
        }
        else {
            throw std::runtime_error("Don't know this optimizer");
        }

        //
        // DMRG Sweep optimization
        // -----------------------
        for (int sweep=init_sweep; sweep < parms["nsweeps"]; ++sweep)
        {
            optimizer->sweep(sweep, Both);
        }

    }

    typename Matrix::value_type energy()
    {
        return expval(mps,mpo) + mpo.getCoreEnergy();
    }

};

#endif
