/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
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

#ifndef OPTIMIZE_H
#define OPTIMIZE_H

#include <boost/random.hpp>
#if not defined(WIN32) && not defined(WIN64)
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
#include "dmrg/utils/placement.h"

template<class Matrix, class SymmGroup>
struct SiteProblem
{
    SiteProblem(Boundary<typename storage::constrained<Matrix>::type, SymmGroup> const & left_,
                Boundary<typename storage::constrained<Matrix>::type, SymmGroup> const & right_,
                MPOTensor<Matrix, SymmGroup> const & mpo_,
                boost::shared_ptr<contraction::Engine<Matrix, typename storage::constrained<Matrix>::type, SymmGroup> > eng_)
    : left(left_)
    , right(right_)
    , mpo(mpo_) 
    , eng(eng_)
    {
    }
    
    Boundary<typename storage::constrained<Matrix>::type, SymmGroup> const & left;
    Boundary<typename storage::constrained<Matrix>::type, SymmGroup> const & right;
    MPOTensor<Matrix, SymmGroup> const & mpo;
    boost::shared_ptr<contraction::Engine<Matrix, typename storage::constrained<Matrix>::type, SymmGroup> > eng;
    double ortho_shift;
};

#define BEGIN_TIMING(name) \
now = boost::chrono::high_resolution_clock::now();
#define END_TIMING(name) \
then = boost::chrono::high_resolution_clock::now(); \
maquis::cout << "Time elapsed in " << name << ": " << boost::chrono::duration<double>(then-now).count() << std::endl;

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
                   int site=0)
    : mps(mps_)
    , mpo(mpo_)
    , parms(parms_)
    , stop_callback(stop_callback_)
    {
        // Initialize the contraction engine
        contr = engine_factory<Matrix, typename storage::constrained<Matrix>::type, SymmGroup>(parms);

        std::size_t L = mps.length();
        #ifdef AMBIENT_TRACKING
        ambient::overseer::log::region("optimizer_base::optimizer_base");
        for(int i = 0; i < L; ++i) ambient_track_array(mps, i);
        #endif
        
        mps.canonize(site);
        for(int i = 0; i < mps.length(); ++i)
            Storage::evict(mps[i]);

        northo = parms_["n_ortho_states"];
        maquis::cout << "Expecting " << northo << " states to orthogonalize to." << std::endl;

        //if (northo > 0 && parms_["ortho_states"]=="")
        if (northo > 0 && !parms_.is_set("ortho_states"))
            throw std::runtime_error("Parameter \"ortho_states\" is not set\n");

        ortho_mps.resize(northo);
        std::string files_ = parms_["ortho_states"].str();
        std::vector<std::string> files;
        boost::split(files, files_, boost::is_any_of(", "));
        for (int n = 0; n < northo; ++n) {
            maquis::cout << "Loading ortho state " << n << " from " << files[n] << std::endl;
            load(files[n], ortho_mps[n]);
            maquis::cout << "Right end: " << ortho_mps[n][mps.length()-1].col_dim() << std::endl;
        }
        
        init_left_right(mpo, site);
        maquis::cout << "Done init_left_right" << std::endl;
    }
    
    virtual ~optimizer_base() {}
    
    virtual void sweep(int sweep, OptimizeDirection d = Both) = 0;
    
    results_collector const& iteration_results() const { return iteration_results_; }

protected:

    inline void boundary_left_step(MPO<Matrix, SymmGroup> const & mpo, int site)
    {
        left_[site+1] = contr->overlap_mpo_left_step(mps[site], mps[site], left_[site], mpo[site]);
        Storage::pin(left_[site+1]);
        
        for (int n = 0; n < northo; ++n)
            ortho_left_[n][site+1] = contr->overlap_left_step(mps[site], ortho_mps[n][site], ortho_left_[n][site]);
    }
    
    inline void boundary_right_step(MPO<Matrix, SymmGroup> const & mpo, int site)
    {
        right_[site] = contr->overlap_mpo_right_step(mps[site], mps[site], right_[site+1], mpo[site]);
        Storage::pin(right_[site]);
        
        for (int n = 0; n < northo; ++n)
            ortho_right_[n][site] = contr->overlap_right_step(mps[site], ortho_mps[n][site], ortho_right_[n][site+1]);
    }

    void init_left_right(MPO<Matrix, SymmGroup> const & mpo, int site)
    {
        construct_placements(mpo);
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
            #ifdef USE_AMBIENT
            ambient::sync(); // to scale down memory
            #endif
        }
        Storage::evict(left_[site]);
        //tlb.end();

        maquis::cout << "Boundaries are partially initialized...\n";
        
        //Timer trb("Init right boundaries"); trb.begin();
        Storage::drop(right_[L]);
        right_[L] = mps.right_boundary();
        Storage::pin(right_[L]);

        for (int i = L-1; i >= site; --i) {
            Storage::drop(right_[i]);
            boundary_right_step(mpo, i);
            Storage::evict(right_[i+1]);
            #ifdef USE_AMBIENT
            ambient::sync(); // to scale down memory
            #endif
        }
        Storage::evict(right_[site]);
        //trb.end();

        maquis::cout << "Boundaries are fully initialized...\n";
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

    boost::shared_ptr<contraction::Engine<Matrix, typename storage::constrained<Matrix>::type, SymmGroup> > contr;

    std::vector<Boundary<typename storage::constrained<Matrix>::type, SymmGroup> > left_, right_;
    
    /* This is used for multi-state targeting */
    unsigned int northo;
    std::vector< std::vector<block_matrix<typename storage::constrained<Matrix>::type, SymmGroup> > > ortho_left_, ortho_right_;
    std::vector<MPS<Matrix, SymmGroup> > ortho_mps;
};

#include "ss_optimize.hpp"
#include "ts_optimize.hpp"

#endif
