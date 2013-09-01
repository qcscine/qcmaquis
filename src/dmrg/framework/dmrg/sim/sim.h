/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2013 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef APP_SIM_H
#define APP_SIM_H

#include "dmrg_version.h"

#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/stat.h>

#include <boost/filesystem.hpp>
#include <boost/optional.hpp>

#include "utils/data_collector.hpp"

#include "dmrg/utils/DmrgParameters2.h"

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/twositetensor.h"
#include "dmrg/mp_tensors/contractions.h"
#include "dmrg/mp_tensors/contractions.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/mpo_ops.h"
#include "dmrg/mp_tensors/mps_initializers.h"

#include "dmrg/utils/random.hpp"
#include "dmrg/utils/time_stopper.h"
#include "utils/timings.h"

#include "dmrg/models/factory.h"

template <class Matrix, class SymmGroup>
class sim {
public:
    sim(DmrgParameters const &, ModelParameters const &);
    ~sim();
    
    virtual void run() =0;
    
protected:
    typedef std::map<std::string, int> status_type;
    
    virtual std::string results_archive_path(status_type const&) const;
    
    virtual void model_init(boost::optional<int> opt_sweep=boost::optional<int>());
    virtual void mps_init();
    virtual void measure(std::string archive_path, Measurements<Matrix, SymmGroup> const& meas);
    // TODO: can be made const, now only problem are parameters
    
    virtual void checkpoint_simulation(MPS<Matrix, SymmGroup> const& state, status_type const&);
    
protected:
    DmrgParameters parms;
    ModelParameters model;
        
    int init_sweep, init_site;
    bool restore;
    bool dns;
    std::string chkpfile;
    std::string rfile;
    
    time_stopper stop_callback;
    
    Lattice_ptr lat;
    typename model_traits<Matrix, SymmGroup>::model_ptr phys_model;
    Hamiltonian<Matrix, SymmGroup> H;
    Index<SymmGroup> phys;
    typename SymmGroup::charge initc;
    MPS<Matrix, SymmGroup> mps;
    MPO<Matrix, SymmGroup> mpo, mpoc;
    Measurements<Matrix, SymmGroup> measurements;
};

#include "sim.hpp"
#endif
