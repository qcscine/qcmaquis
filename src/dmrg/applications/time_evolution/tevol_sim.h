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

using std::cerr;
using std::cout;
using std::endl;

#include "dmrg/sim/sim.h"

#include "dmrg/sim/te_utils.hpp"
#include "dmrg/mp_tensors/evolve.h"
#include "dmrg/mp_tensors/te.h"

#include "utils/types.h"
//#include "dmrg/utils/noop_storage.h"

using namespace app;


template <class Matrix, class SymmGroup>
class dmrg_tevol_sim : public sim<Matrix, SymmGroup> {
    
    typedef sim<Matrix, SymmGroup> base;
    
    typedef std::vector<MPOTensor<Matrix, SymmGroup> > mpo_t;
    typedef Boundary<Matrix, SymmGroup> boundary_t;
    
    enum TEvol_t {te_uknown, te_nn, te_mpo};
    
public:
    
    dmrg_tevol_sim (DmrgParameters & parms_, ModelParameters & model_)
    : base(parms_, model_)
    , te_type( parse_te_type() )
    , split_H( separate_overlaps(base::H) )
    
    {
        std::cout << split_H.size() << " non overlapping Hamiltonians" << std::endl;
        
        if (te_type == te_nn)
            Unn = getUnn(base::parms.template get<double>("dt"),
                         base::sweep < base::parms.template get<int>("nsweeps_img"));
        else if (te_type == te_mpo)
            Umpo = getUmpo(base::parms.template get<double>("dt"),
                           base::sweep < base::parms.template get<int>("nsweeps_img"));
    }
    
    
    bool do_sweep (Logger& iteration_log)
    {        
        
        if (te_type == te_nn)
            nn_time_evolve(iteration_log);
        else if (te_type == te_mpo)
            mpo_time_evolve(iteration_log);
        
        double energy = expval(base::mps, base::mpo);
        std::cout << "Energy " << energy << std::endl;
        iteration_log << make_log("Energy", energy);
        
        return false;
    }
    
private:
    
    TEvol_t parse_te_type ()
    {
        if (base::parms.template get<std::string>("te_type") == "nn")
            return te_nn;
        else if (base::parms.template get<std::string>("te_type") == "mpo")
            return te_mpo;
        else {
            throw std::runtime_error("Don't know this type of te");
            return te_uknown;
        }
    }
    
    /* nearest-neighbors time-evolution */
    
    std::vector<std::map<std::size_t, block_matrix<Matrix, grp> > >
    getUnn (double dt, bool img)
    {
        typename Matrix::value_type I;
        if (img)
            I = utils::real_identity<typename Matrix::value_type>::value;
        else
            I = utils::imag_identity<typename Matrix::value_type>::value;
        typename Matrix::value_type alpha = -I*dt;
        
        std::vector<std::map<std::size_t, block_matrix<Matrix, grp> > > expH(split_H.size());
        for (int i=0; i<split_H.size(); ++i)
            expH[i] = make_exp_nn(split_H[i], alpha);
        return expH;
    }
    
    void nn_time_evolve(Logger& iteration_log)
    {
        for (int i=0; i < Unn.size(); ++i)
            base::mps = evolve(base::mps, Unn[i],
                               base::parms.template get<std::size_t>("max_bond_dimension"),
                               base::parms.template get<double>("truncation_final"));
    }
    
    
    /* mpo time-evolution */
    
    std::vector<MPO<Matrix, grp> >
    getUmpo(double dt, bool img)
    {
        typename Matrix::value_type I;
        if (img)
            I = utils::real_identity<typename Matrix::value_type>::value;
        else
            I = utils::imag_identity<typename Matrix::value_type>::value;
        typename Matrix::value_type alpha = -I*dt;
        
        std::vector<MPO<Matrix, SymmGroup> > expMPO(split_H.size(), MPO<Matrix, SymmGroup>(base::lat->size()));
        for (int i=0; i<split_H.size(); ++i)
            expMPO[i] = make_exp_mpo(base::lat->size(), split_H[i], alpha);
        return expMPO;
    }
    
    
    void mpo_time_evolve(Logger& iteration_log)
    {
        for (int which = 0; which < Umpo.size(); ++which)
        {
            time_evolve<Matrix, SymmGroup, NoopStorageMaster> evolution(base::mps,
                                                                        Umpo[which],
                                                                        base::parms, nossm);
            for (int k = 0; k < 5; ++k)
                evolution.sweep(base::sweep, iteration_log);
            evolution.finalize();
            base::mps = evolution.get_current_mps();
        }
        
    }
    
    
    TEvol_t te_type;
    std::vector<Hamiltonian<Matrix, SymmGroup> > split_H;
    std::vector<std::map<std::size_t, block_matrix<Matrix, SymmGroup> > > Unn;
    std::vector<MPO<Matrix, grp> > Umpo;
    
    NoopStorageMaster nossm;
};



#endif
