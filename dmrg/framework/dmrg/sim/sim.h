/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
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
#include "utils/timings.h"

#include "dmrg/models/factory.h"

template <class Matrix, class SymmGroup>
class sim {
    
    typedef std::vector<MPOTensor<Matrix, SymmGroup> > mpo_t;
    typedef Boundary<Matrix, SymmGroup> boundary_t;
        
public:
    
    sim (DmrgParameters const &, ModelParameters const &, bool fullinit=true);    
    ~sim();
    
    virtual bool run ();
    
    virtual void measure ();
    
protected:
    
    virtual std::string sweep_archive_path ();
    
    virtual void model_init ();
    virtual void mps_init ();
    
    virtual int do_sweep (double=-1) =0;
    virtual void do_sweep_measure ();
    virtual int advance (int=1, double=-1);
    
protected:
    DmrgParameters parms;
    ModelParameters model;
    
    timeval now, then, snow, sthen;
    
    bool dns;
    int sweep;
    int site;
    std::string chkpfile;
    std::string rfile;
    bool restore;
    
    Lattice_ptr lat;
    typename model_traits<Matrix, SymmGroup>::model_ptr phys_model;
    Hamiltonian<Matrix, SymmGroup> H;
    Index<SymmGroup> phys;
    typename SymmGroup::charge initc;
    MPS<Matrix, SymmGroup> mps;
    MPO<Matrix, SymmGroup> mpo, mpoc, ts_cache_mpo;
    Measurements<Matrix, SymmGroup> measurements;
    Measurements<Matrix, SymmGroup> meas_always;
};

#include "sim.hpp"
#endif
