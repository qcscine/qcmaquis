/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 *                            Bela Bauer <bauerb@comp-phys.org>
 *
 *****************************************************************************/

#include <iostream>
#include <sstream>

#include <boost/shared_ptr.hpp>

#include <alps/hdf5.hpp>


#include "dmrg/block_matrix/detail/alps.hpp"

#include "dmrg/models/factory.h"
#include "dmrg/models/measurements.h"
#include "dmrg/mp_tensors/mpo.h"

#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/ss_optimize.hpp"


namespace maquis { namespace dmrg {
    template <>
    void run_sim<symm>(DmrgParameters & parms, ModelParameters & model_parms)
    {
        /// Parsing model
        boost::shared_ptr<Lattice> lattice;
        boost::shared_ptr<Model<matrix, symm> > model;
        model_parser<matrix, symm>(parms.get<std::string>("lattice_library"), parms.get<std::string>("model_library"), model_parms,
                                   lattice, model);
        
        Hamiltonian<matrix, symm> H = model->H();
        MPO<matrix, symm> mpo = make_mpo(lattice->size(), H);
        
        
        /// Initialize MPS
        MPS<matrix, symm> mps(lattice->size(), parms.get<std::size_t>("init_bond_dimension"),
                              H.get_phys(), model->initc(model_parms),
                              *(model->initializer(parms)));
        
        /// Initialize optimizer
        StreamStorageMaster ssm(parms.get<std::string>("storagedir"));
        ss_optimize<matrix, symm, StreamStorageMaster> optimizer(mps, mpo, parms, ssm);
        int sweeps = parms.get<int>("nsweeps");
        
        /// Optimize
        for (int sweep=0; sweep<parms.get<int>("nsweeps"); ++sweep) {
            Logger iteration_log;
            int exit_site = optimizer.sweep(sweep, iteration_log, Both, -1 /* starting site */, -1 /* runlimit */);
            
            /// Write iteration results
            {
                std::stringstream iteration_path;
                iteration_path << "/simulation/sweep" << sweep;
                
                alps::hdf5::archive h5ar(parms.get<std::string>("resultfile"), "w");
                
                h5ar << alps::make_pvp(iteration_path.str() + "/parameters", parms);
                h5ar << alps::make_pvp(iteration_path.str() + "/parameters", model_parms);
                
                h5ar << alps::make_pvp(iteration_path.str() + "/results", iteration_log);
            }
        }
        mps = optimizer.get_current_mps();
        
        
        /// Measurement on final MPS
        measure_on_mps(mps, *lattice, model->measurements(), parms.get<std::string>("resultfile"));
        
        /// Write parameters
        {
            alps::hdf5::archive h5ar(parms.get<std::string>("resultfile"), "w");
            
            h5ar << alps::make_pvp("/parameters", parms);
            h5ar << alps::make_pvp("/parameters", model_parms);
        }
        
        /// Finalize worker thread
        ssm.sync();
    }
} }
