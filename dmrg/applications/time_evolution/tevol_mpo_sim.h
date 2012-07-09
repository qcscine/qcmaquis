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

#include "tevol_sim.h"

#include "dmrg/sim/te_utils.hpp"
#include "dmrg/mp_tensors/te.h"


// ******   SIMULATION CLASS   ******
template <class Matrix, class SymmGroup>
class dmrg_tevol_mpo_sim : public dmrg_tevol_sim<Matrix, SymmGroup> {
    
    typedef dmrg_tevol_sim<Matrix, SymmGroup> base;
    
public:
    typedef typename base::mpo_t mpo_t;
    typedef typename base::boundary_t boundary_t;
    
    dmrg_tevol_mpo_sim(DmrgParameters const & parms_, ModelParameters const  & model_)
    : base(parms_, model_)
    {
        maquis::cout << "Using MPO time evolution." << std::endl;
    }
    
protected:
    void prepare_te_terms(bool split_hamil=true)
    {
        if (split_hamil) {
            hamils = separate_overlaps(this->H);
            maquis::cout << "Found " << hamils.size() << " non overlapping Hamiltonians." << std::endl;
        }

        typename Matrix::value_type I;
        if (this->sweep < this->parms.template get<int>("nsweeps_img"))
            I = maquis::traits::real_identity<typename Matrix::value_type>::value;
        else
            I = maquis::traits::imag_identity<typename Matrix::value_type>::value;
        typename Matrix::value_type alpha = -I*this->parms.template get<double>("dt");
        
        Uterms.resize(hamils.size());
        for (int i=0; i<hamils.size(); ++i)
            Uterms[i] = make_exp_mpo(this->lat->size(), hamils[i], alpha);
    }
    
    void evolve_time_step(Logger & iteration_log)
    {
        for (int which = 0; which < Uterms.size(); ++which)
        {
            time_evolve<Matrix, SymmGroup, NoopStorageMaster> evolution(this->mps,
                                                                        Uterms[which],
                                                                        this->parms, nossm);
            for (int k = 0; k < 5; ++k)
                evolution.sweep(this->sweep, iteration_log);
            evolution.finalize();
            this->mps = evolution.get_current_mps();
        }
    }

        
private:        
    std::vector<Hamiltonian<Matrix, SymmGroup> > hamils;
    std::vector<MPO<Matrix, SymmGroup> > Uterms;
    
    NoopStorageMaster nossm;
};

#endif
