/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2012 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MAQUIS_DMRG_DMRG_INIT_SIM_H
#define MAQUIS_DMRG_DMRG_INIT_SIM_H

#include "dmrg_version.h"

#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/compression.h"
#include "dmrg/models/factory.h"

template <class Matrix, class SymmGroup>
class dmrg_init {    
public:
    typedef typename SymmGroup::charge charge;
    typedef std::pair<charge, size_t> local_state;
    typedef typename std::vector<local_state>::const_iterator states_iterator;
    
    dmrg_init(DmrgParameters const & parms_, ModelParameters const & model_)
    : parms(parms_)
    , model(model_)
    , chkpfile(parms.get<std::string>("chkpfile"))
    , num_states1(0)
    , num_states2(0)
    , totstates(0)
    {
        maquis::cout << DMRG_VERSION_STRING << std::endl;
        
        // TODO: insert boost::chrono timers
        
        model_parser<Matrix, SymmGroup>(parms.get<std::string>("lattice_library"), parms.get<std::string>("model_library"),
                                        model, lat, phys_model);
        initc = phys_model->initc(model);
        phys = phys_model->H().get_phys();
        L = lat->size();
    }
    
    void build()
    {
        build_fast();
        
        // final join of mps1 and mps2
        if (num_states2 > 0) {
            if (num_states1 == parms.get<size_t>("init_bond_dimension"))
                mps1.normalize_left();
            mps2.normalize_left();
            mps1 = join(mps1, mps2, std::sqrt(num_states1), std::sqrt(num_states2));
            mps1.normalize_left();
            mps1 = compression::l2r_compress(mps1, parms.get<size_t>("init_bond_dimension"), parms.get<double>("truncation_initial"));
            mps2 = MPS<Matrix, SymmGroup>();
            num_states1 += num_states2;
            num_states2 = 0;
        }
        
        mps1.normalize_left();
        
        // write parameters and mps
        alps::hdf5::archive h5ar(chkpfile, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
        h5ar << alps::make_pvp("/parameters", parms);
        h5ar << alps::make_pvp("/parameters", model);
        h5ar << alps::make_pvp("/version", DMRG_VERSION_STRING);
        h5ar << alps::make_pvp("/state", mps1);
        h5ar << alps::make_pvp("/status/sweep", 0);
    }
    
private:
    // slow version looping over all basis states, discarding those with N!=initc
    void build_slow()
    {
        std::vector<local_state> alllocal;
        for (size_t i=0; i<phys.size(); ++i)
            for (size_t j=0; j<phys[i].second; ++j)
                alllocal.push_back( local_state(phys[i].first, j) );
        
        std::vector<states_iterator> it(L);
        for (size_t i=0; i<L; ++i)
            it[i] = alllocal.begin();
        
        std::vector<local_state> basis;
        while (it[0] != alllocal.end()) {
            std::vector<local_state> state(L);
            for (size_t i=0; i<L; ++i)
                state[i] = *(it[i]);
            charge N = std::accumulate(state.begin(), state.end(), SymmGroup::IdentityCharge,
                                       boost::bind(static_cast<charge(*)(charge,charge)>(&SymmGroup::fuse), _1,  boost::bind(&local_state::first, _2)) );
            if (N == initc)
                add_state(state);
            
            ++it[L-1];
            for (int i = L-1; (i > 0) && (it[i] == alllocal.end()); --i) {
                it[i] = alllocal.begin();
                ++it[i-1];
            }
        }
    }
    
    // faster version looping only over basis states with N == initc
    void build_fast()
    {
        std::vector<local_state> alllocal;
        for (size_t i=0; i<phys.size(); ++i)
            for (size_t j=0; j<phys[i].second; ++j)
                alllocal.push_back( local_state(phys[i].first, j) );
        
        std::vector<states_iterator> it(L);
        for (size_t i=0; i<L; ++i)
            it[i] = alllocal.begin();
        
        std::vector<local_state> basis;
        while (it[0] != alllocal.end()) {
            charge N = SymmGroup::IdentityCharge;
            for (size_t i=0; i<L; ++i)
                 N = SymmGroup::fuse(N, it[i]->first);
            
            if (N == initc)
                permutate_states(it);
            
            ++it[L-1];
            for (int i = L-1; (i > 0) && (it[i] == alllocal.end()); --i) {
                if (++it[i-1] != alllocal.end())
                    for (int j=i; j<L; ++j)
                        it[j] = it[i-1];
            }
        }
    }
    
    void permutate_states(std::vector<states_iterator> its)
    {
        std::vector<local_state> state(L);
        do {
            for (size_t i=0; i<L; ++i)
                state[i] = *(its[i]);
            add_state(state);
        } while ( next_permutation(its.begin(), its.end()) );
    }
    
    MPS<Matrix, SymmGroup> state_mps(std::vector<local_state> const & state)
    {
        MPS<Matrix, SymmGroup> mps(state.size());
        
        Index<SymmGroup> curr_i;
        curr_i.insert(std::make_pair(SymmGroup::IdentityCharge, 1));
        size_t curr_b = 0;
        for (int i=0; i<state.size(); ++i)
        {
            charge newc = SymmGroup::fuse(curr_i[0].first, state[i].first);
            size_t news = (i == state.size()-1) ? 1 : phys[phys.position(state[i].first)].second;
            Index<SymmGroup> new_i;
            new_i.insert(std::make_pair(newc, news));
            ProductBasis<SymmGroup> left(phys, curr_i);
            mps[i] = MPSTensor<Matrix, SymmGroup>(phys, curr_i, new_i, false, 0);
            size_t b_in = left(state[i].first, curr_i[0].first) + state[i].second * curr_i[0].second + curr_b;
            size_t b_out = (i == state.size()-1) ? 0 : state[i].second;
            mps[i].make_left_paired();
            mps[i].data()(SymmGroup::fuse(curr_i[0].first, state[i].first), new_i[0].first)(b_in, b_out) = 1.;
            curr_i = new_i;
            curr_b = state[i].second;
        }
        mps.normalize_left();
        return mps;
    }
    
    void add_state(std::vector<local_state> const & state)
    {
        MPS<Matrix, SymmGroup> & curr = (num_states1 < parms.get<size_t>("init_bond_dimension")) ? mps1 : mps2;
        size_t & num_states = (num_states1 < parms.get<size_t>("init_bond_dimension")) ? num_states1 : num_states2;
        
        MPS<Matrix, SymmGroup> temp = state_mps(state);
        if (curr.length() > 1)
            curr = join(curr, temp);
        else
            swap(curr, temp);
        num_states += 1;
        totstates += 1;
        
        if (num_states2 > parms.get<size_t>("init_bond_dimension")) {
            if (num_states1 == parms.get<size_t>("init_bond_dimension"))
                mps1.normalize_left();
            mps2.normalize_left();
            mps1 = join(mps1, mps2, std::sqrt(num_states1), std::sqrt(num_states2));
            mps1.normalize_left();
            mps1 = compression::l2r_compress(mps1, parms.get<size_t>("init_bond_dimension"), parms.get<double>("truncation_initial"));
            mps2 = MPS<Matrix, SymmGroup>();
            num_states1 += num_states2;
            num_states2 = 0;
        }
    }
    
    
private:
    DmrgParameters parms;
    ModelParameters model;
    
    std::string chkpfile;
    
    Lattice_ptr lat;
    typename model_traits<Matrix, SymmGroup>::model_ptr phys_model;
    Index<SymmGroup> phys;
    charge initc;
    size_t L;
    MPS<Matrix, SymmGroup> mps1, mps2;
    size_t num_states1, num_states2, totstates;
};



#endif
