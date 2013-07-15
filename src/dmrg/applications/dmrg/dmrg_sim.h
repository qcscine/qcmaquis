/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef APP_DMRG_SIM_H
#define APP_DMRG_SIM_H

#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/stat.h>

#include <boost/shared_ptr.hpp>

#include "dmrg/sim/sim.h"

#include "dmrg/mp_tensors/optimize.h"

//template <class Matrix, class SymmGroup>
//struct do_fix_density {
//    
//    static void apply (MPS<Matrix, SymmGroup> & mps, std::vector<block_matrix<Matrix, SymmGroup> > const & dens_ops,
//                       typename SymmGroup::charge initc)
//    { }
//};
//template <class Matrix>
//struct do_fix_density<Matrix, TwoU1> {
//    
//    static void apply (MPS<Matrix, TwoU1> & mps, std::vector<block_matrix<Matrix, TwoU1> > const & dens_ops,
//                       TwoU1::charge initc)
//    {
//        std::vector<std::vector<double> > dens(2, std::vector<double>(mps.length(), 0.));
//        std::fill(dens[0].begin(), dens[0].end(), double(initc[0])/(mps.length()/2));
//        std::fill(dens[1].begin(), dens[1].end(), double(initc[1])/(mps.length()/2));
//        
//        fix_density(mps, dens_ops, dens);
//    }
//};

template <class Matrix, class SymmGroup>
class dmrg_sim : public sim<Matrix, SymmGroup> {
    
    typedef sim<Matrix, SymmGroup> base;
    typedef optimizer_base<Matrix, SymmGroup, storage::disk> opt_base_t;
    
    typedef std::vector<MPOTensor<Matrix, SymmGroup> > mpo_t;
    typedef Boundary<Matrix, SymmGroup> boundary_t;
    
public:
    
    dmrg_sim (DmrgParameters & parms_, ModelParameters & model_, bool measure_only = false)
    : base(parms_, model_)
    {
        if ( !measure_only &&  base::sweep < base::parms.template get<int>("nsweeps") )
        {
            if (parms_["optimization"] == "singlesite")
            {
                optimizer = 
                boost::shared_ptr<opt_base_t> ( new ss_optimize<Matrix, SymmGroup, storage::disk>
                                               (base::mps, base::mpoc,
                                                base::parms) );
            } 
            
            else if(parms_["optimization"] == "twosite")
            {
                optimizer =
                boost::shared_ptr<opt_base_t> ( new ts_optimize<Matrix, SymmGroup, storage::disk>
                                               (base::mps, base::mpoc, base::ts_cache_mpo,
                                                base::parms) );
            }

            else
                throw std::runtime_error("Don't know this optimizer");
        }
    }
    
    int do_sweep (double time_limit = -1)
    {
        int exit_site = optimizer->sweep(base::sweep, Both,
                                         base::site, time_limit);
        storage::disk::sync();
        
        base::mps = optimizer->get_current_mps();
        
        return exit_site;
    }

    ~dmrg_sim()
    {
        storage::disk::sync();
    }    
private:
    
    boost::shared_ptr<opt_base_t> optimizer;
    
};

#endif
