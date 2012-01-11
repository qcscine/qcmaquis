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

#include "types/dense_matrix/dense_matrix.h"
#include "types/dense_matrix/matrix_interface.hpp"
#include "types/dense_matrix/resizable_matrix_interface.hpp"
#include "types/dense_matrix/algorithms.hpp"
#include "types/dense_matrix/matrix_algorithms.hpp"
#include "types/dense_matrix/dense_matrix_blas.hpp"
#include "types/dense_matrix/aligned_allocator.h"

#ifdef USE_MTM
#include "types/mt_matrix/algorithms.hpp"
#endif

#ifdef USE_GPU
#include <cublas.h>
#endif

#include <alps/hdf5.hpp>

#include "dmrg/utils/DmrgParameters2.h"

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/contractions.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/mpo_ops.h"
#include "dmrg/mp_tensors/mps_initializers.h"

#include "dmrg/utils/stream_storage.h"
#include "dmrg/utils/logger.h"
#include "dmrg/utils/random.hpp"
#include "utils/timings.h"

#include "dmrg/models/factory.h"


namespace app {
    
    template <class Matrix, class SymmGroup>
    class sim {
        
        typedef std::vector<MPOTensor<Matrix, SymmGroup> > mpo_t;
        typedef Boundary<Matrix, SymmGroup> boundary_t;
        
    public:
        
        sim (DmrgParameters const &, ModelParameters const &, bool fullinit=true);    
        ~sim();
        
        virtual void run ();
        virtual int do_sweep (Logger&, double=-1) =0;
        virtual void do_sweep_measure (Logger&);
        
        virtual void measure ();
        
    protected:
        
        virtual void model_init ();
        virtual void mps_init ();
        
        
    private:
        
        mps_initializer<Matrix, SymmGroup> * initializer_factory(BaseParameters & params);
        
        
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
        MPO<Matrix, SymmGroup> mpo, mpoc;
        Measurements<Matrix, SymmGroup> measurements;
        Measurements<Matrix, SymmGroup> meas_always;
        
        StreamStorageMaster ssm;
        
    private:
        
    };
    
} //namespace app


#include "sim.hpp"


#endif
