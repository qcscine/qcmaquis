/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef APP_DMRG_TEVOL_SIM_H
#define APP_DMRG_TEVOL_SIM_H

#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/stat.h>

#include "dmrg/sim/sim.h"

#include "dmrg/sim/te_utils.hpp"

template <class Matrix, class SymmGroup>
class dmrg_tevol_sim : public sim<Matrix, SymmGroup> {
    
    typedef sim<Matrix, SymmGroup> base;
    
public:
    typedef std::vector<MPOTensor<Matrix, SymmGroup> > mpo_t;
    typedef Boundary<Matrix, SymmGroup> boundary_t;
        
    dmrg_tevol_sim (DmrgParameters const & parms_, ModelParameters const  & model_)
    : base(parms_, model_, false)
    , parms_orig(parms_)
    , model_orig(model_)
    , initialized(false)
    {
        base::parms = parms_orig.get_at_index("t", base::sweep);
        base::model = model_orig.get_at_index("t", base::sweep);
        
        base::model_init();
        base::mps_init();
    }
    
    int advance (int nsteps, double time_limit)
    {
        if (this->sweep == 0 || !initialized)
        {
            int pc = 0, mc = 0;
            this->parms = parms_orig.get_at_index("t", this->sweep, &pc);
            this->model = model_orig.get_at_index("t", this->sweep, &mc);
            if (mc > 0 || pc > 0)
                this->model_init();
            prepare_te_terms();
            initialized = true;
        } else if (this->parms.template get<int>("update_each") > -1 && (this->sweep % this->parms.template get<int>("update_each")) == 0)
        {
            int pc = 0, mc = 0;
            this->parms = parms_orig.get_at_index("t", this->sweep, &pc);
            this->model = model_orig.get_at_index("t", this->sweep, &mc);
            if (mc > 0 || pc > 0) {
                this->model_init();
                prepare_te_terms();
            }
        } else {
            if (this->sweep == this->parms.template get<int>("nsweeps_img"))
                // since this is just a change in the time step, there is
                // no need to split the hamiltonian in non-overlapping terms.
                prepare_te_terms(false);
        }
        
        // apply time evolution operators
        evolve_ntime_steps(nsteps);
        
        double energy = maquis::real(expval(this->mps, this->mpo));
        maquis::cout << "Energy " << energy << std::endl;
        storage::log << std::make_pair("Energy", energy);
        
        return -1; // no early exit
    }

    int do_sweep (double time_limit = -1)
    {
        throw std::runtime_error("do_sweep not implemented for time evolution.");
        return -1;
    }
    
protected:
    virtual void prepare_te_terms(bool split_hamil=true) =0;
    virtual void evolve_time_step() =0;
    
    virtual void evolve_ntime_steps(int nsteps)
    {
        int ns = base::sweep + nsteps;
        for (; base::sweep < ns; ++base::sweep)
        {
            this->parms.set("sweep", base::sweep);
            evolve_time_step();
        }
        // base::sweep = ns!
        --base::sweep;
    }
    
    DmrgParameters parms_orig;
    ModelParameters model_orig;
    bool initialized;
};

#endif
