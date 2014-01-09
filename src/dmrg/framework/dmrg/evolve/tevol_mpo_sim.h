/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef APP_DMRG_TEVOL_MPO_SIM_H
#define APP_DMRG_TEVOL_MPO_SIM_H

#include <cmath>

#include "dmrg/utils/storage.h"
#include "dmrg/evolve/te_utils.hpp"
#include "dmrg/mp_tensors/mpo_contractor_ss.h"
#include "dmrg/utils/results_collector.h"

// ******   SIMULATION CLASS   ******
template <class Matrix, class SymmGroup>
class mpo_evolver {
    typedef term_descriptor<typename Matrix::value_type> term_t;
public:
    mpo_evolver(DmrgParameters * parms_, MPS<Matrix, SymmGroup> * mps_,
                Lattice const& lattice_, Model<Matrix, SymmGroup> const& model_,
                int init_sweep=0)
    : parms(parms_)
    , mps(mps_)
    , lattice(lattice_) // shallow copy
    , model(model_) // shallow copy
    , sweep_(init_sweep)
    , hamils(separate_hamil_terms(model.hamiltonian_terms()))
    {
        maquis::cout << "Using MPO time evolution." << std::endl;
        
        maquis::cout << "Found " << hamils.size() << " non overlapping Hamiltonians." << std::endl;
        
        /// compute the time evolution gates
        prepare_te_terms();
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
            Uterms[i] = make_exp_mpo(lattice, model, hamils[i], alpha);
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
        for (int which = 0; which < Uterms.size(); ++which)
        {
            int maxiter = 6; double tol = 1e-6;
            mpo_contractor_ss<Matrix, SymmGroup, storage::nop> evolution(*mps, Uterms[which], (*parms));
            for (int k = 0; k < maxiter; ++k) {
                std::pair<double,double> eps = evolution.sweep(sweep_);
                double rel_error = std::abs( (eps.first-eps.second) / eps.second );
                if (rel_error < tol)
                    break;
            }
            evolution.finalize();
            *mps = evolution.get_current_mps();
        }
    }

    
private:        
    DmrgParameters * parms;
    MPS<Matrix, SymmGroup> * mps;
    Lattice lattice;
    Model<Matrix,SymmGroup> model;
    int sweep_;
    
    results_collector iteration_results_;
    
    std::vector<std::vector<term_t> > hamils;
    std::vector<MPO<Matrix, SymmGroup> > Uterms;
};

#endif
