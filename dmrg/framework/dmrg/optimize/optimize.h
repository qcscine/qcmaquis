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

#include "siteproblem.h"
#include "ietl_lanczos_solver.h"
#include "ietl_jacobi_davidson.h"

#include "dmrg/utils/BaseParameters.h"
#include "dmrg/utils/results_collector.h"
#include "dmrg/utils/storage.h"
#include "dmrg/utils/time_limit_exception.h"
#include "dmrg/utils/parallel/placement.hpp"
#include "dmrg/utils/checks.h"
#include "dmrg/optimize/CorrectionEquation/correctionequation.h"
#include "dmrg/optimize/Finalizer/finalizer.h"
#include "dmrg/optimize/Orthogonalizer/orthogonalizer.h"


#include "dmrg/optimize/POverlap/partial_overlap.h"
#include "dmrg/optimize/siteproblem.h"
#include "dmrg/optimize/singlesitevs.h"
#include "dmrg/optimize/utils/bound_database.h"


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

// +--------------------+
//  OPTIMIZER_BASE CLASS
// +--------------------+
//
// Virtual class to be used to inherit classess whose scope is optimizing an MPS.
// From here, the ss_optimizer and ts_optimizer are inherited

template<class Matrix, class SymmGroup, class Storage>
class optimizer_base
{
    typedef contraction::Engine<Matrix, typename storage::constrained<Matrix>::type, SymmGroup>   contr ;
    typedef Boundary<typename storage::constrained<Matrix>::type, SymmGroup>                      boundary ;
    typedef std::vector< boundary >                                                               boundaries_type ;
    typedef std::vector< boundaries_type >                                                        boundaries_vector ;
    typedef typename partial_overlap<Matrix,SymmGroup>::partial_overlap                           partial_overlap ;
    typedef typename bound_database< MPS<Matrix, SymmGroup>, boundaries_type>::bound_database     bound_database ;
    typedef typename std::vector< std::pair<float, size_t> >                                      sorter_type ;
    typedef typename OptimizationAlgorithm<SiteProblem<Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup>,
             CorrectionEquation<SiteProblem<Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup> > >::OptimizationAlgorithm OptimizationAlgorithm;
    typedef typename Orthogonalizer<SingleSiteVS<Matrix, SymmGroup> >::Orthogonalizer             Orthogonalizer;
    typedef typename Finalizer<SiteProblem<Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup> >::Finalizer Finalizer;
public:
    optimizer_base(std::vector< MPS<Matrix, SymmGroup> > & mps_vector_ ,
                   MPO<Matrix, SymmGroup> const & mpo_,
                   BaseParameters & parms_,
                   boost::function<bool ()> stop_callback_,
                   int site = 0,
                   bool use_ortho = true, // ignore orthogonal states, is used e.g. in derived classes for the local Hamiltonian matrix element measurement
                   boost::optional<MPS<Matrix, SymmGroup> const &> aux_mps_ = boost::none, // Not needed for the optimiser, required only in the derived LocalHamiltonianInitialiser class for calculating sigma vectors
                   boost::optional<std::vector< MPSTensor<Matrix, SymmGroup> > &> mps_guess_ = boost::none,
                   boost::optional<std::vector< MPS<Matrix, SymmGroup> > &> mps_partial_overlap_ = boost::none,
                   boost::optional<MPO<Matrix, SymmGroup> const &> mpo_squared_ = boost::none
                   )
    : mps_vector(mps_vector_)
    , mps_guess_o(mps_guess_)
    , mpo(mpo_)
    , mpo_squared_o(mpo_squared_)
    , aux_mps_o(aux_mps_)
    , n_root_(mps_vector.size())
    , omega_shift_(0.)
    , parms(parms_)
    , stop_callback(stop_callback_)
    , northo(use_ortho ? parms["n_ortho_states"] : 0)
    , omega_vec(0)
    , do_root_homing_(false)
    , do_stateaverage_ (false)
    , do_shiftandinvert_(false)
    , do_H_squared_(false)
    , update_omega(false)
    , update_each_omega(false)
    , i_activate_last_overlap_(parms["activate_last_overlap"])
    , i_activate_constant_omega_(parms["activate_constant_omega"])
   	, is_folded_(false)
    , chebyshev_shift_(parms["chebyshev_shift"])
    , reshuffle_variance_(false)
    , track_variance_(false)
    {
        // Corrector
        init_algorithm(mps_vector_, parms_) ;
        // Standard options
        L_ = mps_vector[0].length() ;
        sorter_.resize(n_root_) ;
        // Chebyshev projection
        if (parms_["chebyshev_filter"] == "yes")
            do_chebyshev_ = true ;
        else if (parms_["chebyshev_filter"] == "no")
            do_chebyshev_ = false ;
        else
            throw std::runtime_error("Parameter chebyshev_filter not recognized") ;
        //
        // Shift-and-invert paramters
        // --------------------------
        double omega = parms["ietl_si_omega"] ;
        omega_shift_ = parms["si_omega_shift"] ;
	      i_activate_last_overlap_   = parms["activate_last_overlap"] ;
	      i_activate_constant_omega_ = parms["activate_constant_omega"] ;
        if (parms["ietl_si_operator"] == "yes" || parms["ietl_si_operator"] == "folded") {
            do_shiftandinvert_ = true ;
            omega_vec.resize(n_root_, omega) ;
            if (parms["si_omega_schedule" ]== "constant")
                update_omega = false ;
            else if (parms["si_omega_schedule"] == "update")
                update_omega = true ;
            else
                throw std::runtime_error("Scheduler for omega update not recognized") ;
        }
         //
        // Initialization of the MPS
        // -------------------------
        // The single SA are loaded inside the mps_vector, mps is initialized
        // with the average of the MPSs
        do_stateaverage_ = parms_["n_states_sa"].as<int>() > 0 ;
        for (size_t k = 0; k < n_root_; k++)
            order.push_back(k) ;

        for (int i = 0 ; i < n_root_; ++i ) {
            if (use_ortho)
                mps_vector[i].canonize(site);
            for (int j = 0; j < mps_vector[i].length(); ++j)
                Storage::evict(mps_vector[i][j]);
        }
        // Variance tracking
        if (parms["track_variance"] == "yes")
            track_variance_ = true ;
        if (parms["reshuffle_variance"] == "yes")
            reshuffle_variance_ = true ;
            do_H_squared_ = track_variance_ || reshuffle_variance_ ;
        //TODO ALB Hardcoded for the moment
        if (do_H_squared_)
            throw std::runtime_error("Variance tracking for electronic structure calculations not implemented") ;
        //
        // Root-homing criteria
        // --------------------
        poverlap_vec_.resize(0) ;
        if (parms_["ietl_diag_homing_criterion"] == "") {
            do_root_homing_   = false ;
            root_homing_type_ = 0 ;
        } else if (parms_["ietl_diag_homing_criterion"] == "input" || parms_["ietl_diag_homing_criterion"] == "both") {
            do_root_homing_   = true ;
            if (parms_["ietl_diag_homing_criterion"] == "input")
                root_homing_type_ = 1 ;
            else
                root_homing_type_ = 3 ;
            // Once the states to be followed have been initialized, builds the vector of
            // the partial overlap objects
            assert(mps_partial_overlap_);
            assert (mps_vector.size() == mps_partial_overlap_->size()) ;
            for (size_t i = 0 ; i < n_root_ ; i++) {
                partial_overlap poverlap(mps_vector[i], mps_partial_overlap_.get()[i]) ;
                poverlap_vec_.push_back(poverlap) ;
            }
        } else if (parms_["ietl_diag_homing_criterion"] == "last") {
            do_root_homing_   = true ;
            root_homing_type_ = 2 ;
        } else if (parms_["ietl_diag_homing_criterion"] == "variance") {
            throw std::runtime_error("Variance tracking NYI") ;
            //do_root_homing_   = true ;
            //root_homing_type_ = 4 ;
        } else {
            throw std::runtime_error("Root homing criterion not recognized") ;
        }
        //
        // Orthogonal states
        // -----------------
        // Loads the states already converged against which to perform the constrained optimization
        // for(int i = 0; i < mps.length(); ++i)
        //     Storage::evict(mps[i]);

        if (use_ortho)
        {
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
                maquis::checks::right_end_check(files[n], ortho_mps[n], mps_vector[0][mps_vector[0].length()-1].col_dim()[0].first);
                maquis::cout << "Right end: " << ortho_mps[n][mps_vector[0].length()-1].col_dim() << std::endl;
            }
        }

        //
        // Initialization of the boundaries
        // --------------------------------
        // Check also the consistency in the definition of sa_alg_
        sa_alg_ = parms["sa_algorithm"] ;
        init_bound() ;
        left_sa_.resize(n_bound_) ;
        right_sa_.resize(n_bound_) ;
        if (do_H_squared_) {
            assert(mpo_squared_o);
            left_squared_sa_.resize(n_bound_) ;
            right_squared_sa_.resize(n_bound_) ;
        }
        for (int i = 0 ; i < n_bound_ ; i++) {
            left_sa_[i].resize(mpo.length()+1) ;
            right_sa_[i].resize(mpo.length()+1) ;
            if (do_H_squared_) {
                left_squared_sa_[i].resize(mpo.length()+1) ;
                right_squared_sa_[i].resize(mpo.length()+1) ;
            }
        }
        boundaries_database_ = bound_database(mps_vector, left_sa_, right_sa_, left_squared_sa_, right_squared_sa_,
                                              sa_alg_, do_H_squared_) ;
        init_left_right(mpo, site);
        maquis::cout << "Done init_left_right" << std::endl;
    }
    //
    virtual ~optimizer_base() {}
    virtual void sweep(int sweep, OptimizeDirection d = Both) = 0;
    std::vector<results_collector> const& iteration_results() const { return iteration_results_; }
    //
protected:
    // +--------------+
    //  INIT_ALGORITHM
    // +--------------+
    // Routine used to initialize the parameters of the optimizer
    void init_algorithm(std::vector< MPS<Matrix, SymmGroup> > & mps_vector,
                        BaseParameters & parms)
    {
        // Initialize "empty" objects

        maquis::cout << "\n" ;
        maquis::cout << " Parameters of the simulation\n" ;
        maquis::cout << " ----------------------------\n\n" ;
        // Extract number of sites
        L_ = mps_vector[0].length() ;
        maquis::cout << " Number of sites - " << L_ << "\n" ;
        // S&I parameters
        if (parms["ietl_si_operator"] == "no") {
            // -- Standard simualtion --
            maquis::cout << " S&I formulation - deactivated\n" ;
            finalizer_ = std::unique_ptr<Finalizer>(new StandardFinalizer<SiteProblem<Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup> >());
            correction_equation.set_standard() ;
            maquis::cout << " Corrector op.   - standard\n" ;
            orthogonalizer_ = std::unique_ptr<Orthogonalizer>(new GS_ortho<SingleSiteVS<Matrix, SymmGroup> > (parms["ietl_ortho_refine"] == "yes"));
            maquis::cout << " Orthogonalizer  - Gram-Schmidt\n" ;
        } else if (parms["ietl_si_operator"] == "yes") {
            // -- S&I simulation --
            maquis::cout << " S&I formulation - activated\n" ;
            double omega = parms["ietl_si_omega"]  - mpo.getCoreEnergy();
            maquis::cout << " Omega parameter - " << omega << "\n" ;
            omega_shift_ = parms["si_omega_shift"] ;
            maquis::cout << " Shift parameter - " << omega_shift_ << "\n" ;
            omega_vec.resize(n_root_, omega) ;
            if (parms["si_omega_schedule" ]== "constant") {
                update_omega = false ;
                maquis::cout << " Omega update    - NO\n" ;
            } else if (parms["si_omega_schedule"] == "update") {
                update_omega = true ;
                update_each_omega = true ;
                maquis::cout << " Omega update    - YES\n" ;
            } else if (parms["si_omega_schedule"] == "updatesweep") {
                update_omega = true ;
                maquis::cout << " Omega update    - YES each sweep\n" ;
            } else {
                throw std::runtime_error("Scheduler for omega update not recognized") ;
            }
            // -- Set corrector --
            if (parms["ietl_corrector"] == "notSI") {
                finalizer_ = std::unique_ptr<Finalizer>(new StandardSIFinalizer<SiteProblem<Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup> >());
                orthogonalizer_ = std::unique_ptr<Orthogonalizer>(new GS_ortho_mod<SingleSiteVS<Matrix, SymmGroup> >(parms["ietl_ortho_refine"] == "yes"));
                maquis::cout << " Orthogonalizer  - Gram-Schmidt\n";
            } else if (parms["ietl_corrector"] == "skew") {
                finalizer_ = std::unique_ptr<Finalizer>(new SkewFinalizer<SiteProblem<Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup> >());
                correction_equation.set_skew() ;
                orthogonalizer_ = std::unique_ptr<Orthogonalizer>(new BI_ortho<SingleSiteVS<Matrix, SymmGroup> >());
                maquis::cout << " Orthogonalizer  - Biorthogonal\n" ;
            } else if (parms["ietl_corrector"] == "SI") {
                finalizer_ = std::unique_ptr<Finalizer>(new ModifiedFinalizer<SiteProblem<Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup> >());
                correction_equation.set_modified() ;
                orthogonalizer_ = std::unique_ptr<Orthogonalizer>(new GS_ortho<SingleSiteVS<Matrix, SymmGroup> >(parms["ietl_ortho_refine"] == "yes"));
                maquis::cout << " Orthogonalizer  - Gram-Schmidt\n" ;
            } else {
                throw std::runtime_error("Corrector operator not recognized") ;
            }
		    } else if (parms["ietl_si_operator"] == "folded") {
            maquis::cout << " S&I formulation - activated\n" ;
            double omega = parms["ietl_si_omega"]  - mpo.getCoreEnergy();
            maquis::cout << " Omega parameter - " << omega << "\n" ;
            omega_shift_ = parms["si_omega_shift"] ;
            maquis::cout << " Shift parameter - " << omega_shift_ << "\n" ;
            omega_vec.resize(n_root_, omega) ;
            if (parms["si_omega_schedule" ]== "constant") {
                update_omega = false ;
                maquis::cout << " Omega update    - NO\n" ;
            } else if (parms["si_omega_schedule"] == "update") {
                update_omega = true ;
                update_each_omega = true ;
                maquis::cout << " Omega update    - YES\n" ;
            } else if (parms["si_omega_schedule"] == "updatesweep") {
                update_omega = true ;
                maquis::cout << " Omega update    - YES each sweep\n" ;
            } else {
                throw std::runtime_error("Scheduler for omega update not recognized") ;
            }
			      is_folded_ = true ;
			      do_H_squared_ = true ;
                  finalizer_ = std::unique_ptr<Finalizer>(new FoldedFinalizer<SiteProblem<Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup> >());
			      correction_equation.set_folded() ;
			      orthogonalizer_ = std::unique_ptr<Orthogonalizer>(new GS_ortho<SingleSiteVS<Matrix, SymmGroup> >(parms["ietl_ortho_refine"] == "yes"));

        } else {
            throw std::runtime_error(" S&I modality not recognized\n") ;
        }
        // -- Set preconditioning --
        if (parms["ietl_precondition"] == "yes" || parms["eigensolver"] == std::string("IETL_DAVIDSON")) {
            correction_equation.activate_preconditioner();
            maquis::cout << " Preconditioner  - Activated\n";
        } else if (parms["ietl_precondition"] == "no" ) {
            maquis::cout << " Preconditioner  - Deactivated\n";
        } else {
            throw std::runtime_error("Precondition parameter not recognized");
        }

        // set the microoptimizer
        if (parms["eigensolver"] == std::string("IETL_JCD")) {
            if (parms["ietl_microopt_maxiter"] == 0)
            {
                micro_optimizer = std::unique_ptr<OptimizationAlgorithm>(new OP_optimizer<SiteProblem<Matrix, SymmGroup>,
                                          SingleSiteVS<Matrix, SymmGroup>,
                                          CorrectionEquation<SiteProblem<Matrix, SymmGroup>,
                                          SingleSiteVS<Matrix, SymmGroup> > >
                                          (correction_equation,
                                          parms["ietl_microopt_abstol"], parms["ietl_microopt_reltol"]
                                          ));
                maquis::cout << " Microoptimizer  - JD/Preconditioner only" << std::endl;
            }
            else
            {
                if (parms["ietl_microoptimizer"] == "GMRES") {
                    micro_optimizer = std::unique_ptr<OptimizationAlgorithm>(new GMRES_optimizer<SiteProblem<Matrix, SymmGroup>,
                                            SingleSiteVS<Matrix, SymmGroup>,
                                            CorrectionEquation<SiteProblem<Matrix, SymmGroup>,
                                            SingleSiteVS<Matrix, SymmGroup> > >
                                            (correction_equation,
                                            parms["ietl_microopt_abstol"], parms["ietl_microopt_reltol"],
                                            parms["ietl_microopt_maxiter"], parms["ietl_microopt_restart"]
                                            ));

                    maquis::cout << " Microoptimizer  - GMRES\n";
                } else if (parms["ietl_microoptimizer"] == "CG") {
                    micro_optimizer = std::unique_ptr<OptimizationAlgorithm>(new CG_optimizer<SiteProblem<Matrix, SymmGroup>,
                                            SingleSiteVS<Matrix, SymmGroup>,
                                            CorrectionEquation<SiteProblem<Matrix, SymmGroup>,
                                            SingleSiteVS<Matrix, SymmGroup> > >
                                            (correction_equation,
                                            parms["ietl_microopt_abstol"], parms["ietl_microopt_reltol"],
                                            parms["ietl_microopt_maxiter"], parms["ietl_microopt_restart"]
                                            ));
                    maquis::cout << " Microoptimizer  - CG\n";
                } else if (parms["ietl_microoptimizer"] == "BICGS") {
                    micro_optimizer = std::unique_ptr<OptimizationAlgorithm>(new BICGS_optimizer<SiteProblem<Matrix, SymmGroup>,
                                            SingleSiteVS<Matrix, SymmGroup>,
                                            CorrectionEquation<SiteProblem<Matrix, SymmGroup>,
                                            SingleSiteVS<Matrix, SymmGroup> > >
                                            (correction_equation,
                                            parms["ietl_microopt_abstol"], parms["ietl_microopt_reltol"],
                                            parms["ietl_microopt_maxiter"], parms["ietl_microopt_restart"]
                                            ));
                    maquis::cout << " Microoptimizer  - BICGS\n";
                } else
                    throw std::runtime_error("ietl_microoptimizer parameter not recognized");
            }
        } else if (parms["eigensolver"] == std::string("IETL_DAVIDSON")) {
            micro_optimizer = std::unique_ptr<OptimizationAlgorithm>(new OP_optimizer<SiteProblem<Matrix, SymmGroup>,
                                          SingleSiteVS<Matrix, SymmGroup>,
                                          CorrectionEquation<SiteProblem<Matrix, SymmGroup>,
                                          SingleSiteVS<Matrix, SymmGroup> > >
                                          (correction_equation,
                                          parms["ietl_microopt_abstol"], parms["ietl_microopt_reltol"]
                                          ));
            maquis::cout << " Microoptimizer  - Davidson/Preconditioner only\n";
        } else {
            throw std::runtime_error("I don't know this eigensolver.");
        }
        // Verbosity options
        if (parms["ietl_microopt_verbose"] == "yes")
            micro_optimizer->activate_verbosity() ;
        maquis::cout << "\n" ;
    }
    // +---------------+
    //  INIT_LEFT_RIGHT
    // +---------------+
    // Routine used to initialize the boundaries
    void init_left_right(MPO<Matrix, SymmGroup> const & mpo, int site)
    {
        const MPS<Matrix,SymmGroup> & mps = mps_vector[0];
        parallel::construct_placements(mpo);
        // -- Initialization of the various vectors --
        // Allocate and parially initialize the space for the left/right orthogonalization vectors
        ortho_left_.resize(northo);
        ortho_right_.resize(northo);
        for (int n = 0; n < northo; ++n) {
            ortho_left_[n].resize(L_+1) ;
            ortho_right_[n].resize(L_+1) ;
            ortho_left_[n][0] = mps.left_boundary()[0] ;
            ortho_right_[n][L_] = mps.right_boundary()[0] ;
        }
        // Allocate and parially initialize the space for the left/right orthogonalization vectors
        // for the state-average calculation
        vec_sa_left_.resize(n_root_) ;
        vec_sa_right_.resize(n_root_) ;
        for (int k = 0; k < n_root_; k++){
            vec_sa_left_[k].resize(n_root_) ;
            vec_sa_right_[k].resize(n_root_) ;
            for (int h = 0; h < n_root_; h++) {
                vec_sa_left_[k][h].resize(L_+1) ;
                vec_sa_right_[k][h].resize(L_+1) ;
                // Partial initialization
                vec_sa_left_[k][h][0]   = (*(boundaries_database_.get_mps(k))).left_boundary()[0];
                vec_sa_right_[k][h][L_] = (*(boundaries_database_.get_mps(h))).right_boundary()[0];
            }
        }
        // Complete initialization and builds all the boundaries objects
        if (!aux_mps_o) // if aux_mps is not present, construct boundaries normally
            for (size_t i = 0; i < n_bound_; i++)
            {
                (*(boundaries_database_.get_boundaries_left(i, false)))[0] =
                    (*(boundaries_database_.get_mps(i))).left_boundary();
                if (do_H_squared_)
                    (*(boundaries_database_.get_boundaries_left(i, true)))[0] =
                        (*(boundaries_database_.get_mps(i))).left_boundary();
            }
        else
        {
            for (size_t i = 0; i < n_bound_; i++)
            {
                (*(boundaries_database_.get_boundaries_left(i, false)))[0] =
                    make_left_boundary(aux_mps_o.get(), *(boundaries_database_.get_mps(i)));

                // for completeness only
                if (do_H_squared_)
                    (*(boundaries_database_.get_boundaries_left(i, true)))[0] =
                        make_left_boundary(aux_mps_o.get(), *(boundaries_database_.get_mps(i)));
            }
        }

        for (size_t i = 0; i < site; ++i) {
            boundary_left_step(mpo, i);
            parallel::sync();
        }
        maquis::cout << "Initialised left boundaries...\n";
        //
        if (!aux_mps_o) // if aux_mps is not present, construct boundaries normally
            for (size_t i = 0; i < n_bound_; i++) {
                (*(boundaries_database_.get_boundaries_right(i, false)))[L_] =
                    (*(boundaries_database_.get_mps(i))).right_boundary();
                if (do_H_squared_)
                    (*(boundaries_database_.get_boundaries_right(i, true)))[L_] =
                        (*(boundaries_database_.get_mps(i))).right_boundary();
            }
        else
        {
            for (size_t i = 0; i < n_bound_; i++)
            {
                (*(boundaries_database_.get_boundaries_right(i, false)))[L_] =
                    make_right_boundary(aux_mps_o.get(), *(boundaries_database_.get_mps(i)));

                // for completeness only
                if (do_H_squared_)
                    (*(boundaries_database_.get_boundaries_right(i, true)))[L_] =
                        make_right_boundary(aux_mps_o.get(), *(boundaries_database_.get_mps(i)));
            }
        }

        for (int i = L_-1; i >= site; --i) {
            boundary_right_step(mpo, i);
            parallel::sync(); // to scale down memory
        }
        //trb.end();
        maquis::cout << "Initialised right boundaries...\n";
    }
    // +----------+
    //  INIT_BOUND
    // +----------+
    void init_bound()
    {
        // Decides the number of boundaries to be stored
        if (n_root_ > 0) {
            if (sa_alg_ >= 0 || sa_alg_ == -3 || sa_alg_ == -1) {
                n_bound_ = 1 ;
            } else if (sa_alg_ == -2) {
                n_bound_ = n_root_ ;
            }
        } else {
            n_bound_ = 1 ;
        }
    }
    // +------------------+
    //  BOUNDARY_LEFT_STEP
    // +------------------+
    // Shifts the boundary one site to the left
    inline void boundary_left_step(MPO<Matrix, SymmGroup> const & mpo, int site)
    {
        // Variables definition
        const MPS<Matrix,SymmGroup> & mps = mps_vector[0];

        MPSTensor<Matrix, SymmGroup> tmp ;
        // Shifts the boundaries
        for (size_t i = 0 ; i < n_bound_ ; i++) {

            // use the bra MPS for constructing the boundaries if it is present
            const MPS<Matrix, SymmGroup> & mps_for_boundary = (aux_mps_o) ? aux_mps_o.get() : *(boundaries_database_.get_mps(i));
            (*(boundaries_database_.get_boundaries_left(i, false)))[site+1] =
                    contr::overlap_mpo_left_step(mps_for_boundary[site],
                                                *(boundaries_database_.get_mps(i,site)),
                                                (*(boundaries_database_.get_boundaries_left(i, false)))[site],
                                                 mpo[site]);
            if (do_H_squared_)
                (*(boundaries_database_.get_boundaries_left(i, true)))[site+1] =
                        contr::overlap_mpo_left_step(mps_for_boundary[site],
                                                     *(boundaries_database_.get_mps(i,site)),
                                                    (*(boundaries_database_.get_boundaries_left(i, true)))[site],
                                                     mpo_squared()[site]);
        }
        // Updates the orthogonal vectors
        for (int n = 0; n < northo; ++n)
            ortho_left_[n][site+1] = contr::overlap_left_step(mps[site], ortho_mps[n][site], ortho_left_[n][site]);
        for (int i = 0; i < n_root_ ; i++)
            for (int j = 0; j < n_root_; j++)
                if ( i != j )
                    vec_sa_left_[i][j][site+1] = contr::overlap_left_step(*(boundaries_database_.get_mps(i,site)),
                                                                          *(boundaries_database_.get_mps(j,site)),
                                                                          vec_sa_left_[i][j][site]);
    }
    // +-------------------+
    //  BOUNDARY_RIGHT_STEP
    // +-------------------+
    // Shifts the boundary one site to the right
    inline void boundary_right_step(MPO<Matrix, SymmGroup> const & mpo, int site)
    {
        const MPS<Matrix,SymmGroup> & mps = mps_vector[0];
        // Variables definition
        MPSTensor<Matrix, SymmGroup> tmp ;
        // Shifts the boundaries
        for (size_t i = 0 ; i < n_bound_ ; i++) {
            const MPS<Matrix, SymmGroup> & mps_for_boundary = (aux_mps_o) ? aux_mps_o.get() : *(boundaries_database_.get_mps(i));
            (*(boundaries_database_.get_boundaries_right(i, false)))[site] =
                    contr::overlap_mpo_right_step(mps_for_boundary[site],
                                                  *(boundaries_database_.get_mps(i,site)),
                                                 (*(boundaries_database_.get_boundaries_right(i, false)))[site+1],
                                                  mpo[site]);
            if (do_H_squared_)
                (*(boundaries_database_.get_boundaries_right(i, true)))[site] =
                        contr::overlap_mpo_right_step(mps_for_boundary[site],
                                                      *(boundaries_database_.get_mps(i,site)),
                                                     (*(boundaries_database_.get_boundaries_right(i, true)))[site+1],
                                                      mpo_squared()[site]);
        }
        // Updates the orthogonal vectors
        for (int n = 0; n < northo; ++n)
            ortho_right_[n][site] = contr::overlap_right_step(mps[site], ortho_mps[n][site], ortho_right_[n][site+1]);
        for (int i = 0; i < n_root_ ; i++)
            for (int j = 0; j < n_root_; j++)
                if ( i != j )
                    vec_sa_right_[i][j][site] = contr::overlap_right_step(*(boundaries_database_.get_mps(i,site)),
                                                                          *(boundaries_database_.get_mps(j,site)),
                                                                          vec_sa_right_[i][j][site+1]);
    }
    // +----------+
    //  GET_CUTOFF
    // +----------+
    // Simple routine to set the threshold for the truncation process
    double get_cutoff(int sweep) const
    {
        double cutoff;
        if (sweep >= parms.template get<int>("ngrowsweeps"))
            cutoff = parms.template get<double>("truncation_final");
        else
            cutoff = log_interpolate(parms.template get<double>("truncation_initial"), parms.template get<double>("truncation_final"), parms.template get<int>("ngrowsweeps"), sweep);
        std::cout << "Cutoff - " << cutoff << std::endl ;
        return cutoff;
    }

    // +----------+
    //   GET_MMAX
    // +----------+
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
    // +----------------------------------------+
    //   UPDATES THE ENERGY ORDER OF THE STATES
    // +----------------------------------------+
    void update_order(const std::vector< std::pair<float,std::size_t> >& sorter)
    {
        for (size_t i = 0; i < n_root_; i++)
            order[i]  = sorter[i].second ;
    }
    // +------------------------+
    //   UPDATES THE PARAMETERS
    // +------------------------+
    void update_parameters(const int& sweep)
    {
        if (sweep == i_activate_constant_omega_) {
            std::cout << "Omega update deactivated" << std::endl ;
            update_omega = false;
        }
        if (sweep == i_activate_last_overlap_) {
            std::cout << "Last homing deactivated" << std::endl;
            root_homing_type_ = 2;
        }
    }
    // +-----------------------+
    //   COMPUTES THE VARIANCE
    // +-----------------------+
    double compute_variance(MPSTensor<Matrix, SymmGroup> const &input,
                            std::size_t const& i_state,
                            std::size_t const& idx)
    {
        assert (track_variance_) ;
        MPSTensor<Matrix, SymmGroup> y;
        y = contraction::Engine<Matrix, Matrix, SymmGroup>::site_hamil2(input,
                                                                        left_squared_sa_[i_state][idx],
                                                                        right_squared_sa_[i_state][idx+1],
                                                                        mpo_squared()[idx]);
        return ietl::dot(input, y);
    }
    // +----------+
    //  ATTRIBUTES
    // +----------+
    std::vector< results_collector > iteration_results_;
    CorrectionEquation< SiteProblem<Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup> > correction_equation ;
    std::unique_ptr<Finalizer> finalizer_ ;
    std::unique_ptr<Orthogonalizer> orthogonalizer_ ;
    std::unique_ptr<OptimizationAlgorithm> micro_optimizer ;
    MPS<Matrix, SymmGroup> mps_average;
    MPO<Matrix, SymmGroup> const& mpo;
    BaseParameters & parms ;
    boost::function<bool ()> stop_callback;
    boundaries_vector left_sa_, right_sa_, left_squared_sa_, right_squared_sa_ ;
    sorter_type sorter_ ;
    double chebyshev_shift_ ;

    /* This is used for multi-state targeting */
    unsigned int northo;
    std::vector< std::vector<block_matrix<typename storage::constrained<Matrix>::type, SymmGroup> > > ortho_left_, ortho_right_;
    std::vector< std::vector< std::vector<block_matrix<typename storage::constrained<Matrix>::type, SymmGroup> > > > vec_sa_left_, vec_sa_right_;
    std::vector<MPS<Matrix, SymmGroup> > ortho_mps;
    // Root-homing procedure
    bool do_root_homing_ , do_stateaverage_ , do_shiftandinvert_, do_chebyshev_ ;
    int i_activate_last_overlap_, root_homing_type_ ;
    std::vector<partial_overlap> poverlap_vec_;
    // Energy-specific diagonalization
    bool update_omega, update_each_omega, is_folded_ ;
    double omega_shift_ ;
    int i_activate_constant_omega_ ;
    std::vector<double> omega_vec ;
    // Variance tracking
    bool reshuffle_variance_, track_variance_, do_H_squared_ ;
    // State average
    bound_database boundaries_database_ ;
    std::vector< MPS<Matrix, SymmGroup> > & mps_vector ;
    std::vector< std::size_t > order ;
    int sa_alg_ ;
    size_t n_bound_ , L_ , n_root_ ;

    boost::optional<std::vector<MPSTensor<Matrix, SymmGroup> > &> mps_guess_o ;
    boost::optional<MPO<Matrix, SymmGroup> const&> mpo_squared_o;

    // The bra MPS required in the derived LocalHamiltonianInitialiser class for calculating sigma vectors
    boost::optional<MPS<Matrix, SymmGroup> const&> aux_mps_o;

    std::vector<MPSTensor<Matrix, SymmGroup> > & mps_guess() { return mps_guess_o.get(); };
    MPO<Matrix, SymmGroup> const& mpo_squared() { return mpo_squared_o.get(); };
};

#include "ss_optimize.hpp"
#include "ts_optimize.hpp"

#endif
