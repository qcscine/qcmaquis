/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef APP_DMRG_TEVOL_MPO_SIM_H
#define APP_DMRG_TEVOL_MPO_SIM_H

#include <cmath>
#include <sys/stat.h>

#include "dmrg/utils/storage.h"
#include "dmrg/sim/te_utils.hpp"
#include "dmrg/mp_tensors/te.h"
#include "dmrg/utils/results_collector.h"

// ******   SIMULATION CLASS   ******
template <class Matrix, class SymmGroup>
class mpo_evolver {
public:
    mpo_evolver(DmrgParameters * parms_, MPS<Matrix, SymmGroup> * mps_,
                Lattice_ptr lat_, Hamiltonian<Matrix, SymmGroup> const* H_,
                int init_sweep=0)
    : parms(parms_)
    , mps(mps_)
    , lat(lat_)
    , H(H_)
    , sweep_(init_sweep)
    , hamils(separate_overlaps(*H))
    {
        maquis::cout << "Using MPO time evolution." << std::endl;
        
        maquis::cout << "Found " << hamils.size() << " non overlapping Hamiltonians." << std::endl;
    }
    
    void prepare_te_terms()
    {
        double dt = (*parms)["dt"];
        typename Matrix::value_type I;
        if (sweep_ < (*parms)["nsweeps_img"])
            I = maquis::traits::real_identity<typename Matrix::value_type>::value;
        else
            I = maquis::traits::imag_identity<typename Matrix::value_type>::value;
        typename Matrix::value_type alpha = -I*dt;
        
        Uterms.resize(hamils.size());
        for (int i=0; i<hamils.size(); ++i)
            Uterms[i] = make_exp_mpo(lat->size(), hamils[i], alpha);
    }
    
    void operator()(int nsteps)
    {
        iteration_results_.clear();
        
        int ns = sweep_ + nsteps;
        for (int i=sweep_; i < ns; ++i) {
            sweep_ = i;
            (*parms).set("sweep", sweep_);
            evolve_time_step();
        }
        assert(sweep_ == ns-1);
    }
    
    int sweep() const
    {
        return sweep_;
    }
    
    results_collector const& iteration_results() const
    {
        return iteration_results_;
    }
    
private:
    void evolve_time_step()
    {
        // TODO: use iteration_results
        for (int which = 0; which < Uterms.size(); ++which)
        {
            time_evolve<Matrix, SymmGroup, storage::nop> evolution(*mps,
                                                                   Uterms[which],
                                                                   (*parms));
            for (int k = 0; k < 5; ++k)
                evolution.sweep(sweep_);
            evolution.finalize();
            *mps = evolution.get_current_mps();
        }
    }

    
private:        
    DmrgParameters * parms;
    MPS<Matrix, SymmGroup> * mps;
    Lattice_ptr lat;
    Hamiltonian<Matrix, SymmGroup> const * H;
    int sweep_;
    
    results_collector iteration_results_;
    
    std::vector<Hamiltonian<Matrix, SymmGroup> > hamils;
    std::vector<MPO<Matrix, SymmGroup> > Uterms;
};

#endif
