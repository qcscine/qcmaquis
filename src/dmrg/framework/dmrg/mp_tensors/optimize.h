/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *
 *****************************************************************************/

#ifndef OPTIMIZE_H
#define OPTIMIZE_H

#include <boost/random.hpp>
#ifndef WIN32
#include <sys/time.h>
#define HAVE_GETTIMEOFDAY
#endif

#include <boost/algorithm/string.hpp>

#include "utils/sizeof.h"

#include "ietl_lanczos_solver.h"
#include "ietl_jacobi_davidson.h"

#include "dmrg/utils/BaseParameters.h"
#include "dmrg/utils/results_collector.h"
#include "dmrg/utils/storage.h"
#include "dmrg/utils/time_limit_exception.h"

template<class Matrix, class SymmGroup>
struct SiteProblem
{
    SiteProblem(Boundary<typename storage::constrained<Matrix>::type, SymmGroup> const & left_,
                Boundary<typename storage::constrained<Matrix>::type, SymmGroup> const & right_,
                MPOTensor<Matrix, SymmGroup> const & mpo_)
    : left(left_)
    , right(right_)
    , mpo(mpo_) 
    {
    }
    
    Boundary<typename storage::constrained<Matrix>::type, SymmGroup> const & left;
    Boundary<typename storage::constrained<Matrix>::type, SymmGroup> const & right;
    MPOTensor<Matrix, SymmGroup> const & mpo;
    double ortho_shift;
};

#ifdef HAVE_GETTIMEOFDAY
#define BEGIN_TIMING(name) \
gettimeofday(&now, NULL);
#define END_TIMING(name) \
gettimeofday(&then, NULL); \
maquis::cout << "Time elapsed in " << name << ": " << then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec) << std::endl;
#else
#define BEGIN_TIMING(name)
#define END_TIMING(name)
#endif

inline double log_interpolate(double y0, double y1, int N, int i)
{
    if (N < 2)
        return y1;
    if (y0 == 0)
        return 0;
    double x = log(y1/y0)/(N-1);
    return y0*exp(x*i);
}

enum OptimizeDirection { Both, LeftOnly, RightOnly };

template<class Matrix, class SymmGroup, class Storage>
class optimizer_base
{
public:
    optimizer_base(MPS<Matrix, SymmGroup> & mps_,
                   MPO<Matrix, SymmGroup> const & mpo_,
                   BaseParameters & parms_,
                   boost::function<bool ()> stop_callback_,
                   int initial_sweep_,
                   int initial_site_)
    : mps(mps_)
    , mpo(mpo_)
    , parms(parms_)
    , stop_callback(stop_callback_)
    , initial_site((initial_site_ < 0) ? 0 : initial_site_)
    {
        std::size_t L = mps.length();
        #ifdef AMBIENT_TRACKING
        ambient::overseer::log::region("optimizer_base::optimizer_base");
        for(int i = 0; i < L; ++i) ambient_track_array(mps, i);
        #endif
        int site = (initial_site < L) ? initial_site : 2*L-initial_site-1;
        
        mps.canonize(site);
        for(int i = 0; i < mps.length(); ++i)
            Storage::evict(mps[i]);

        northo = parms_["n_ortho_states"];
        maquis::cout << "Expecting " << northo << " states to orthogonalize to." << std::endl;
        ortho_mps.resize(northo);
        std::string files_ = parms_["ortho_states"];
        std::vector<std::string> files;
        boost::split(files, files_, boost::is_any_of(", "));
        for (int n = 0; n < northo; ++n) {
            maquis::cout << "Loading ortho state " << n << " from " << files[n] << std::endl;
            storage::archive ar(files[n]);
            ar["/state"] >> ortho_mps[n];
            maquis::cout << "Right end: " << ortho_mps[n][mps.length()-1].col_dim() << std::endl;
        }
        
        init_left_right(mpo, site);
        maquis::cout << "Done init_left_right" << std::endl;
    }
    
    virtual void sweep(int sweep, OptimizeDirection d = Both) = 0;
    
    results_collector const& iteration_results() const { return iteration_results_; }

protected:

    inline void boundary_left_step(MPO<Matrix, SymmGroup> const & mpo, int site)
    {
        left_[site+1] = contraction::overlap_mpo_left_step(mps[site], mps[site], left_[site], mpo[site]);
        Storage::pin(left_[site+1]);
        
        for (int n = 0; n < northo; ++n)
            ortho_left_[n][site+1] = contraction::overlap_left_step(mps[site], ortho_mps[n][site], ortho_left_[n][site]);
    }
    
    inline void boundary_right_step(MPO<Matrix, SymmGroup> const & mpo, int site)
    {
        right_[site] = contraction::overlap_mpo_right_step(mps[site], mps[site], right_[site+1], mpo[site]);
        Storage::pin(right_[site]);
        
        for (int n = 0; n < northo; ++n)
            ortho_right_[n][site] = contraction::overlap_right_step(mps[site], ortho_mps[n][site], ortho_right_[n][site+1]);
    }

    void init_left_right(MPO<Matrix, SymmGroup> const & mpo, int site)
    {
        std::size_t L = mps.length();
        
        left_.resize(mpo.length()+1);
        right_.resize(mpo.length()+1);
        
        ortho_left_.resize(northo);
        ortho_right_.resize(northo);
        for (int n = 0; n < northo; ++n) {
            ortho_left_[n].resize(L+1);
            ortho_right_[n].resize(L+1);
            
            ortho_left_[n][0] = mps.left_boundary()[0];
            ortho_right_[n][L] = mps.right_boundary()[0];
        }
        
        //Timer tlb("Init left boundaries"); tlb.begin();
        Storage::drop(left_[0]);
        left_[0] = mps.left_boundary();
        Storage::pin(left_[0]);
        
        for (int i = 0; i < site; ++i) {
            Storage::drop(left_[i+1]);
            boundary_left_step(mpo, i);
            Storage::evict(left_[i]);
        }
        Storage::evict(left_[site]);
        //tlb.end();
        
        //Timer trb("Init right boundaries"); trb.begin();
        Storage::drop(right_[L]);
        right_[L] = mps.right_boundary();
        Storage::pin(right_[L]);

        for (int i = L-1; i >= site; --i) {
            Storage::drop(right_[i]);
            boundary_right_step(mpo, i);
            Storage::evict(right_[i+1]);
        }
        Storage::evict(right_[site]);
        //trb.end();
    }
    
    double get_cutoff(int sweep) const
    {
        double cutoff;
        if (sweep >= parms.template get<int>("ngrowsweeps"))
            cutoff = parms.template get<double>("truncation_final");
        else
            cutoff = log_interpolate(parms.template get<double>("truncation_initial"), parms.template get<double>("truncation_final"), parms.template get<int>("ngrowsweeps"), sweep);
        return cutoff;
    }

    std::size_t get_Mmax(int sweep) const
    {
        std::size_t Mmax;
        if (parms.is_set("sweep_bond_dimensions")) {
            std::vector<std::size_t> ssizes = parms.template get<std::vector<std::size_t> >("sweep_bond_dimensions");
            if (sweep >= ssizes.size())
                Mmax = *ssizes.rbegin();
            else
                Mmax = ssizes[sweep];
        } else
            Mmax = parms.template get<std::size_t>("max_bond_dimension");
        return Mmax;
    }
    
    
    results_collector iteration_results_;
    
    MPS<Matrix, SymmGroup> & mps;
    MPO<Matrix, SymmGroup> const& mpo;
    
    BaseParameters & parms;
    boost::function<bool ()> stop_callback;
    int initial_site;

    std::vector<Boundary<typename storage::constrained<Matrix>::type, SymmGroup> > left_, right_;
    
    /* This is used for multi-state targeting */
    unsigned int northo;
    std::vector< std::vector<block_matrix<typename storage::constrained<Matrix>::type, SymmGroup> > > ortho_left_, ortho_right_;
    std::vector<MPS<Matrix, SymmGroup> > ortho_mps;
};

#include "ss_optimize.hpp"
#include "ts_optimize.hpp"

#endif
