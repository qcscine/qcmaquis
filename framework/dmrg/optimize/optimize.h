/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *               2017 by Alberto Baiardi <alberto.baiardi@sns.it>
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
#include "ietl_davidson.h"

#include "dmrg/utils/BaseParameters.h"
#include "dmrg/utils/results_collector.h"
#include "dmrg/utils/storage.h"
#include "dmrg/utils/time_limit_exception.h"
#include "dmrg/utils/parallel/placement.hpp"
#include "dmrg/utils/checks.h"

#include "dmrg/optimize/partial_overlap.h"

//
// Site problem structure
// ----------------------
//
// This structure contains left and right boundary, together with the
// MPO tensor of the site of interest

template<class Matrix, class SymmGroup>
struct SiteProblem
{
    SiteProblem(Boundary<typename storage::constrained<Matrix>::type, SymmGroup> const & left_,
                Boundary<typename storage::constrained<Matrix>::type, SymmGroup> const & right_,
                MPOTensor<Matrix, SymmGroup> const & mpo_)
    : left(left_)
    , right(right_)
    , mpo(mpo_) 
    { }
    Boundary<typename storage::constrained<Matrix>::type, SymmGroup> const & left;
    Boundary<typename storage::constrained<Matrix>::type, SymmGroup> const & right;
    MPOTensor<Matrix, SymmGroup> const & mpo;
    double ortho_shift;
};

#define BEGIN_TIMING(name) \
now = boost::chrono::high_resolution_clock::now();
#define END_TIMING(name) \
then = boost::chrono::high_resolution_clock::now(); \
maquis::cout << " Time elapsed in " << name << ": " << boost::chrono::duration<double>(then-now).count() << std::endl;

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

//
// OPTIMIZER_BASE CLASS
// --------------------
//
// Virtual class to be used to inherit classess whose scope is optimizing an MPS.
// From here, the ss_optimizer and ts_optimizer are inherited

template<class Matrix, class SymmGroup, class Storage>
class optimizer_base
{
    typedef contraction::Engine<Matrix, typename storage::constrained<Matrix>::type, SymmGroup>   contr ;
    typedef std::vector<Boundary<typename storage::constrained<Matrix>::type, SymmGroup> >        boundaries_type ;
    typedef std::vector< boundaries_type>                                                         boundaries_vector ;
    typedef typename partial_overlap<Matrix,SymmGroup>::partial_overlap                           partial_overlap ;
public:
    optimizer_base(MPS<Matrix, SymmGroup> & mps_,
                   std::vector< MPS<Matrix, SymmGroup> > & mps_vector_ ,
                   MPO<Matrix, SymmGroup> const & mpo_,
                   BaseParameters & parms_,
                   boost::function<bool ()> stop_callback_,
                   int site=0)
    : mps(mps_)
    , mps_vector(mps_vector_)
    , mpo(mpo_)
    , parms(parms_)
    , stop_callback(stop_callback_)
    , do_root_homing_(false)
    {
        // Standard options
        mps2follow.resize(1) ;
        mps2follow[0].resize(0) ;
        std::size_t L = mps.length();
        // State-average calculation
        n_root_          = mps_vector.size() ;
        do_stateaverage_ = n_root_ > 0 ;
        if (do_stateaverage_) {
            // Reset the first guess for the mps
            for (int k = 0 ; k < L ; k++) {
                mps[k] = mps_vector[0][k] ;
                for (int i = 1; i < n_root_ ; i++)
                    mps[k] += mps_vector[i][k];
                mps[k] /= n_root_ ;
            }
        }
        //TODO ALB TO CHECK IF IT'S USEFUL OR NOT
        mps.canonize(site);
        for (int i = 0; i < mps.length(); ++i)
            Storage::evict(mps[i]);
        for (int i = 0 ; i < n_root_; ++i ) {
            mps_vector[i].canonize(site);
            for (int j = 0; j < mps_vector[i].length(); ++j)
                Storage::evict(mps_vector[i][j]);
        }
        // Root-homing criteria
        poverlap_vec_.resize(0) ;
        if (parms_["ietl_diag_homing_criterion"] == "") {
            do_root_homing_   = false ;
            root_homing_type_ = 0 ;
        } else if (parms_["ietl_diag_homing_criterion"] == "input") {
            do_root_homing_   = true ;
            root_homing_type_ = 1 ;
            if (do_stateaverage_) {
                // Extract the list of states to be followed
                std::vector<std::string> list_sa;
                std::string input_str = parms["follow_mps_stateaverage"].str();
                boost::split(list_sa, input_str, boost::is_any_of("|"));
                if (list_sa.size() != n_root_ )
                    throw std::runtime_error("Number of states to track != from number of SA MPS");
                for (int i = 0; i < list_sa.size(); i++) {
                    std::stringstream ss(list_sa[i]);
                    int ichar;
                    std::vector<int> tmp_vec;
                    while (ss >> ichar) {
                        tmp_vec.push_back(ichar);
                        ss.ignore(1);
                    }
                    if (tmp_vec.size() != L)
                        throw std::runtime_error("ONV to follow has a wrong length");
                    mps2follow.push_back(tmp_vec);
                }
            } else {
                mps2follow.push_back(parms_["follow_basis_state"].as<std::vector<int> >());
                if (mps2follow.size() != L)
                    throw std::runtime_error("ONV to follow has a wrong length");
            }
            // Builds the vector of partial overlap objects
            if (do_stateaverage_) {
                for (size_t i = 0 ; i < n_root_ ; i++){
                    partial_overlap poverlap(mps_vector[i],mps2follow[i]) ;
                    poverlap_vec_.push_back(poverlap) ;
                }
            } else {
                partial_overlap poverlap(mps,mps2follow[0]) ;
                poverlap_vec_.push_back(poverlap) ;
            }
        } else if (parms_["ietl_diag_homing_criterion"] == "last") {
            do_root_homing_   = true ;
            root_homing_type_ = 2 ;
        } else {
            throw std::runtime_error("Root homing criterion not recognized") ;
        }
        // Orthogonal states
        northo = parms_["n_ortho_states"];
        maquis::cout << "Expecting " << northo << " states to orthogonalize to." << std::endl;
        if (northo > 0 && !parms_.is_set("ortho_states"))
            throw std::runtime_error("Parameter \"ortho_states\" is not set\n");
        if (northo > 0 && do_root_homing_)
            throw std::runtime_error("Constrained optimization + SA NYI");
        ortho_mps.resize(northo);
        std::string files_ = parms_["ortho_states"].str();
        std::vector<std::string> files;
        boost::split(files, files_, boost::is_any_of(", "));
        for (int n = 0; n < northo; ++n) {
            maquis::cout << "Loading ortho state " << n << " from " << files[n] << std::endl;
            maquis::checks::symmetry_check(parms, files[n]);
            load(files[n], ortho_mps[n]);
            maquis::checks::right_end_check(files[n], ortho_mps[n], mps[mps.length()-1].col_dim()[0].first);
            maquis::cout << "Right end: " << ortho_mps[n][mps.length()-1].col_dim() << std::endl;
        }
        // Initialization of the MPO
        init_left_right(mpo, site);
        maquis::cout << "Done init_left_right" << std::endl;
    }
    
    virtual ~optimizer_base() {}
    virtual void sweep(int sweep, OptimizeDirection d = Both) = 0;
    results_collector const& iteration_results() const { return iteration_results_; }

protected:

    inline void boundary_left_step(MPO<Matrix, SymmGroup> const & mpo, int site)
    {
        left_[site+1] = contr::overlap_mpo_left_step(mps[site], mps[site], left_[site], mpo[site]);
        for (size_t i = 0 ; i < n_root_ ; i++)
            left_sa_[i][site+1] = contr::overlap_mpo_left_step(mps_vector[i][site], mps_vector[i][site], left_sa_[i][site], mpo[site]);
        for (int n = 0; n < northo; ++n)
            ortho_left_[n][site+1] = contr::overlap_left_step(mps[site], ortho_mps[n][site], ortho_left_[n][site]);
    }
    
    inline void boundary_right_step(MPO<Matrix, SymmGroup> const & mpo, int site)
    {
        right_[site] = contr::overlap_mpo_right_step(mps[site], mps[site], right_[site+1], mpo[site]);
        for (size_t i = 0 ; i < n_root_ ; i++)
            right_sa_[i][site] = contr::overlap_mpo_right_step(mps_vector[i][site], mps_vector[i][site], right_sa_[i][site+1], mpo[site]);
        for (int n = 0; n < northo; ++n)
            ortho_right_[n][site] = contr::overlap_right_step(mps[site], ortho_mps[n][site], ortho_right_[n][site+1]);
    }

    void init_left_right(MPO<Matrix, SymmGroup> const & mpo, int site)
    {
        parallel::construct_placements(mpo);
        std::size_t L = mps.length()  ;
        //
        left_.resize(mpo.length()+1)  ;
        right_.resize(mpo.length()+1) ;
        left_sa_.resize(n_root_) ;
        right_sa_.resize(n_root_) ;
        for (int i = 0 ; i < n_root_ ; i++){
            left_sa_[i].resize(mpo.length()+1) ;
            right_sa_[i].resize(mpo.length()+1) ;
        }
        ortho_left_.resize(northo);
        ortho_right_.resize(northo);
        for (int n = 0; n < northo; ++n) {
            ortho_left_[n].resize(L+1);
            ortho_right_[n].resize(L+1);
            ortho_left_[n][0] = mps.left_boundary()[0];
            ortho_right_[n][L] = mps.right_boundary()[0];
        }
        left_[0] = mps.left_boundary();
        for (int k = 0; k < n_root_ ; ++k)
            left_sa_[k][0] = mps_vector[k].left_boundary();
        for (int i = 0; i < site; ++i) {
            boundary_left_step(mpo, i);
            Storage::evict(left_[i]);
            parallel::sync(); // to scale down memory
        }
        Storage::evict(left_[site]);
        maquis::cout << "Boundaries are partially initialized...\n";
        //
        Storage::drop(right_[L]);
        right_[L] = mps.right_boundary();
        for (int k = 0; k < n_root_ ; ++k)
            right_sa_[k][L] = mps_vector[k].right_boundary();
        Storage::pin(right_[L]);
        for (int i = L-1; i >= site; --i) {
            Storage::drop(right_[i]);
            boundary_right_step(mpo, i);
            Storage::evict(right_[i+1]);
            parallel::sync(); // to scale down memory
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
    // -- GET_MMAX --
    // Routine to get the maximum value for the bond dimension
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
    // +----------+
    //  ATTRIBUTES
    // +----------+
    results_collector iteration_results_;
    MPS<Matrix, SymmGroup> & mps;
    MPO<Matrix, SymmGroup> const& mpo;
    BaseParameters & parms;
    boost::function<bool ()> stop_callback;
    boundaries_type   left_,    right_;
    boundaries_vector left_sa_, right_sa_;
    /* This is used for multi-state targeting */
    unsigned int northo;
    std::vector< std::vector<block_matrix<typename storage::constrained<Matrix>::type, SymmGroup> > > ortho_left_, ortho_right_;
    std::vector<MPS<Matrix, SymmGroup> > ortho_mps;
    // Root-homing procedure
    bool do_root_homing_ , do_stateaverage_ , do_shiftandinvert_ ;
    int root_homing_type_ ;
    std::vector<partial_overlap> poverlap_vec_;
    // State average
    std::vector< std::vector<int> > mps2follow;
    std::vector< MPS<Matrix, SymmGroup> > & mps_vector ;
    int n_root_ ;
};

#include "ss_optimize.hpp"
#include "ts_optimize.hpp"

#endif
