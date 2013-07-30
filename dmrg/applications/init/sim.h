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

#ifdef MAQUIS_OPENMP
#include <omp.h>
#endif

template <class Matrix, class SymmGroup>
class dmrg_init {    
public:
    typedef typename SymmGroup::charge charge;
    typedef std::pair<charge, size_t> local_state;
    typedef typename std::vector<local_state>::const_iterator states_iterator;
    
    dmrg_init(DmrgParameters const & parms_, ModelParameters const & model_)
    : parms(parms_)
    , model(model_)
    , chkpfile(parms["chkpfile"].str())
    , nthreads(1)
    {
        maquis::cout << DMRG_VERSION_STRING << std::endl;
        
        // TODO: insert boost::chrono timers
        
        model_parser<Matrix, SymmGroup>(parms["lattice_library"], parms["model_library"],
                                        model, lat, phys_model);
        initc = phys_model->initc(model);
        phys = phys_model->H().get_phys();
        L = lat->size();
        
#ifdef MAQUIS_OPENMP
        #pragma omp parallel
        {
            nthreads = omp_get_num_threads();
        }
#endif
        num_states1.resize(nthreads, 0);
        num_states2.resize(nthreads, 0);
        totstates.resize(nthreads, 0);
        vmps1.resize(nthreads);
        vmps2.resize(nthreads);
    }
    
    void build()
    {
#ifdef MAQUIS_OPENMP
        #pragma omp parallel
        {
            #pragma omp single nowait
            build_fast();
            #pragma omp barrier
            size_t rank = omp_get_thread_num();
#else
        {
            build_fast();
            size_t rank = 0;
#endif

            // final join of mps1 and mps2
            if (num_states2[rank] > 0) {
                for (int i=0; i<L; ++i)
                    assert(vmps2[rank][i].reasonable());
                maquis::cout << "final merge" << std::endl;
                if (num_states1[rank] == parms.get<size_t>("init_bond_dimension"))
                    vmps1[rank].normalize_left();
                vmps2[rank].normalize_left();
                vmps1[rank] = join(vmps1[rank], vmps2[rank], std::sqrt(num_states1[rank]), std::sqrt(num_states2[rank]));
                // compress & normalize
                vmps1[rank] = compression::l2r_compress(vmps1[rank], parms.get<size_t>("init_bond_dimension"), parms.get<double>("truncation_initial"));
                vmps2[rank] = MPS<Matrix, SymmGroup>();
                num_states1[rank] += num_states2[rank];
                num_states2[rank] = 0;
                maquis::cout << vmps1[rank].description();
            } else if (num_states1[rank] > 0) { // if the thread actually contributed
                vmps1[rank].normalize_left();
            }
        }
        
        MPS<Matrix, SymmGroup> & mps = vmps1[0];
        for (size_t n=0; n<nthreads; ++n) {
            if (num_states2[n] > 0) {
                mps = join(mps, vmps1[n], std::sqrt(num_states1[0]), std::sqrt(num_states1[n]));
                mps.normalize_left();
                mps = compression::l2r_compress(mps, parms.get<size_t>("init_bond_dimension"), parms.get<double>("truncation_initial"));
                num_states1[0] += num_states1[n];
            }
        }
        
        // write parameters and mps
        storage::archive ar(chkpfile, "w");
        ar["/parameters"] << parms;
        ar["/parameters"] << model;
        ar["/version"] << DMRG_VERSION_STRING;
        ar["/state"] << mps;
        ar["/status/sweep"] << 0;
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
            
            if (N == initc) {
                std::stringstream ss;
                ss << "permutations on [";
                for (size_t i=0; i<L; ++i)
                    ss << "  (" << it[i]->first << " - " << it[i]->second << ")";
                ss << "  ]" << std::endl;
                maquis::cout << ss.str();

                std::vector<states_iterator> * tmp = new std::vector<states_iterator>(it);
#ifdef MAQUIS_OPENMP
                #pragma omp task firstprivate(tmp) // ugly because gcc doesn't support copy of complex objects
#endif
                {
                    permutate_states(*tmp);
                    delete tmp;
                }
            }
            ++it[L-1];
            for (int i = L-1; (i > 0) && (it[i] == alllocal.end()); --i) {
                if (++it[i-1] != alllocal.end())
                    for (int j=i; j<L; ++j)
                        it[j] = it[i-1];
            }
        }
#ifdef MAQUIS_OPENMP
        #pragma omp taskwait
#endif
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
#ifdef MAQUIS_OPENMP
        size_t rank = omp_get_thread_num();
#else
        size_t rank = 0;
#endif
        
        MPS<Matrix, SymmGroup> & mps1=vmps1[rank], mps2=vmps2[rank];
        size_t nstates1=num_states1[rank], nstates2=num_states2[rank], tstates=totstates[rank];
        
        MPS<Matrix, SymmGroup> & curr = (nstates1 < parms.get<size_t>("init_bond_dimension")) ? vmps1[rank] : vmps2[rank];
        size_t & nstates = (nstates1 < parms.get<size_t>("init_bond_dimension")) ? nstates1 : nstates2;
        
        MPS<Matrix, SymmGroup> temp = state_mps(state);
        if (curr.length() > 1)
            curr = join(curr, temp);
        else
            swap(curr, temp);
        nstates += 1;
        tstates += 1;
        
        if (nstates2 > parms.get<size_t>("init_bond_dimension")) {
            if (nstates1 == parms.get<size_t>("init_bond_dimension"))
                mps1.normalize_left();
            mps2.normalize_left();
            mps1 = join(mps1, mps2, std::sqrt(nstates1), std::sqrt(nstates2));
            // compress & normalize
            mps1 = compression::l2r_compress(mps1, parms.get<size_t>("init_bond_dimension"), parms.get<double>("truncation_initial"));
            mps2 = MPS<Matrix, SymmGroup>();
            nstates1 += nstates2;
            nstates2 = 0;
        }

        num_states1[rank] = nstates1;
        num_states2[rank] = nstates2;
        totstates[rank] =  tstates;
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
    std::vector<MPS<Matrix, SymmGroup> > vmps1, vmps2;
    std::vector<size_t> num_states1, num_states2, totstates;
    size_t nthreads;
};



#endif
