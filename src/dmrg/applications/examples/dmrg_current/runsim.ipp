/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 *                            Bela Bauer <bauerb@comp-phys.org>
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

#include <iostream>
#include <sstream>

#include <boost/shared_ptr.hpp>

#include "dmrg/block_matrix/detail/alps.hpp"

#include "dmrg/models/factory.h"
#include "dmrg/models/measurements.h"
#include "dmrg/mp_tensors/mpo.h"

#include "dmrg/mp_tensors/mps.h"
#include "dmrg/optimize/ss_optimize.hpp"


namespace maquis { namespace dmrg {
    template <>
    void run_sim<symm>(DmrgParameters & parms, ModelParameters & model_parms)
    {
        /// Parsing model
        boost::shared_ptr<Lattice> lattice;
        boost::shared_ptr<Model<matrix, symm> > model;
        model_parser<matrix, symm>(parms["lattice_library"], parms["model_library"], model_parms,
                                   lattice, model);
        
        Hamiltonian<matrix, symm> H = model->H();
        MPO<matrix, symm> mpo = make_mpo(lattice->size(), H);
        
        
        /// Initialize MPS
        MPS<matrix, symm> mps(lattice->size(), parms.get<std::size_t>("init_bond_dimension"),
                              H.get_phys(), model->initc(model_parms),
                              *(model->initializer(parms)));
        
        /// Initialize optimizer
        storage::setup(parms);
        ss_optimize<matrix, symm, storage::disk> optimizer(mps, mpo, parms);
        int sweeps = parms.get<int>("nsweeps");
        
        /// Optimize
        for (int sweep=0; sweep<parms.get<int>("nsweeps"); ++sweep) {
            int exit_site = optimizer.sweep(sweep, Both, -1 /* starting site */, -1 /* runlimit */);
            
            /// Write iteration results
            {
                std::stringstream iteration_path;
                iteration_path << "/simulation/sweep" << sweep;
                
                storage::archive ar(parms["resultfile"].str(), "w");
                
                ar[iteration_path.str() + "/parameters"] << parms;
                ar[iteration_path.str() + "/parameters"] << model_parms;
                
                ar[iteration_path.str() + "/results"] << storage::log;
            }
        }
        mps = optimizer.get_current_mps();
        
        
        /// Measurement on final MPS
        measure_on_mps(mps, *lattice, model->measurements(), parms["resultfile"]);
        
        /// Write parameters
        {
            storage::archive ar(parms["resultfile"].str(), "w");
            
            ar["/parameters"] << parms;
            ar["/parameters"] << model_parms;
        }
        
        /// Finalize worker thread
        storage::disk::sync();
    }
} }
