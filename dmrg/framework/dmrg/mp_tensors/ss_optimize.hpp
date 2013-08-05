/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef SS_OPTIMIZE_H
#define SS_OPTIMIZE_H

#include "dmrg/mp_tensors/mpo_ops.h"
#include "dmrg/mp_tensors/optimize.h"


template<class Matrix, class SymmGroup, class Storage>
class ss_optimize : public optimizer_base<Matrix, SymmGroup, Storage>
{
public:

    typedef optimizer_base<Matrix, SymmGroup, Storage> base;
    using base::mpo;
    using base::mps;
    using base::left_;
    using base::right_;
    using base::parms;
    using base::iteration_results_;
    using base::initial_site;
    using base::stop_callback;

    ss_optimize(MPS<Matrix, SymmGroup> & mps_,
                MPO<Matrix, SymmGroup> const & mpo_,
                BaseParameters & parms_,
                boost::function<bool ()> stop_callback_,
                int initial_sweep_ = 0,
                int initial_site_ = 0)
    : base(mps_, mpo_, parms_, stop_callback_, initial_sweep_, initial_site_)
    { }
    
    void sweep(int sweep, OptimizeDirection d = Both)
    {
        timeval sweep_now, sweep_then;
        gettimeofday(&sweep_now, NULL);
        
        iteration_results_.clear();
        
        std::size_t L = mps.length();
        
        int _site = 0, site = 0;
        if (initial_site != -1) {
            _site = initial_site;
            site = (_site < L) ? _site : 2*L-_site-1;
        }
        
//        if (parms["beta_mode"] && sweep == 0 && resume_at < L) {
//            int site = (resume_at == -1) ? 0 : resume_at;
//            mpo = zero_after(mpo_orig, site+2);
//            mps.canonize(site);
//            this->init_left_right(mpo, site);
//        }
        
        Storage::prefetch(left_[site]);
        Storage::prefetch(right_[site+1]);
        
#ifndef NDEBUG
        maquis::cout << mps.description() << std::endl;
#endif
        for (; _site < 2*L; ++_site) {
            
            int lr;
            if (_site < L) {
                site = _site;
                lr = 1;
            } else {
                site = 2*L-_site-1;
                lr = -1;
            }
            
            if (lr == -1 && site == L-1) {
                maquis::cout << "Syncing storage" << std::endl;
                Storage::sync();
            }
            
            maquis::cout << "Sweep " << sweep << ", optimizing site " << site << std::endl;
            
//            mps[site].make_left_paired();
            
            // MD: some changes needed to re-enable it.
//            if (parms.["beta_mode"]) {
//                if (sweep == 0 && lr == 1) {
//                    mpo = zero_after(mpo_orig, 0);
//                    if (site == 0)
//                        this->init_left_right(mpo, 0);
//                } else if (sweep == 0 && lr == -1 && site == L-1) {
//                    mpo = mpo_orig;
//                    //this->init_left_right(mpo, site);
//                }
//            }
            
            
            Storage::fetch(left_[site]);
            Storage::fetch(right_[site+1]);
            
            if (lr == +1) Storage::prefetch(left_[site+1]);
            else          Storage::prefetch(right_[site]);
            
            assert( left_[site].reasonable() );    // in case something went wrong
            assert( right_[site+1].reasonable() ); // in case something went wrong
            
            
//            maquis::cout << "My size: " << std::endl;
//            maquis::cout << "  left_: " << utils::size_of(left_.begin(), left_.end())/1024.0/1024 << std::endl;
//            maquis::cout << "  right_: " << utils::size_of(right_.begin(), right_.end())/1024.0/1024 << std::endl;
//            maquis::cout << "  MPS: " << utils::size_of(mps.begin(), mps.end())/1024.0/1024 << std::endl;
//            maquis::cout << "  MPS[i]: " << utils::size_of(mps[site])/1024.0/1024 << std::endl;
            
            //SiteProblem<Matrix, SymmGroup> sp(mps[site], left_[site], right_[site+1], mpo[site]);
            
            timeval now, then;

            std::pair<double, MPSTensor<Matrix, SymmGroup> > res;
            SiteProblem<Matrix, SymmGroup> sp(left_[site], right_[site+1], mpo[site]);
            
            /// Compute orthogonal vectors
            std::vector<MPSTensor<Matrix, SymmGroup> > ortho_vecs(base::northo);
            for (int n = 0; n < base::northo; ++n) {
                ortho_vecs[n] = contraction::site_ortho_boundaries(mps[site], base::ortho_mps[n][site],
                                                                    base::ortho_left_[n][site], base::ortho_right_[n][site+1]);
            }

            if (d == Both ||
                (d == LeftOnly && lr == -1) ||
                (d == RightOnly && lr == +1))
            {
                if (parms["eigensolver"] == std::string("IETL")) {
                    BEGIN_TIMING("IETL")
                    res = solve_ietl_lanczos(sp, mps[site], parms);
                    END_TIMING("IETL")
                } else if (parms["eigensolver"] == std::string("IETL_JCD")) {
                    BEGIN_TIMING("JCD")
                    res = solve_ietl_jcd(sp, mps[site], parms, ortho_vecs);
                    END_TIMING("JCD")
                } else {
                    throw std::runtime_error("I don't know this eigensolver.");
                }
 
                mps[site] = res.second;
            }
            
#ifndef NDEBUG
            // Caution: this is an O(L) operation, so it really should be done only in debug mode
            for (int n = 0; n < base::northo; ++n)
                maquis::cout << "MPS overlap: " << overlap(mps, base::ortho_mps[n]) << std::endl;
#endif
            
            maquis::cout << "Energy " << lr << " " << res.first << std::endl;
//            maquis::cout << "Energy check " << maquis::real(expval(mps, mpo)) << std::endl;
            
            iteration_results_["Energy"] << res.first;
            
            double alpha;
            int ngs = parms.template get<int>("ngrowsweeps"), nms = parms.template get<int>("nmainsweeps");
            if (sweep < ngs)
                alpha = parms.template get<double>("alpha_initial");
            else if (sweep < ngs + nms)
                alpha = parms.template get<double>("alpha_main");
            else
                alpha = parms.template get<double>("alpha_final");
            
            double cutoff = this->get_cutoff(sweep);
            std::size_t Mmax = this->get_Mmax(sweep);
            truncation_results trunc;
            
            if (lr == +1) {
                if (site < L-1) {
                    maquis::cout << "Growing, alpha = " << alpha << std::endl;
                    trunc = mps.grow_l2r_sweep(mpo[site], left_[site], right_[site+1],
                                               site, alpha, cutoff, Mmax);
                } else {
                    block_matrix<Matrix, SymmGroup> t = mps[site].normalize_left(DefaultSolver());
                    if (site < L-1)
                        mps[site+1].multiply_from_left(t);
                }
                
                
                Storage::drop(left_[site+1]); // left_[site+1] is outdated
                this->boundary_left_step(mpo, site); // creating left_[site+1]
            } else if (lr == -1) {
                if (site > 0) {
                    maquis::cout << "Growing, alpha = " << alpha << std::endl;
                    // Invalid read occurs after this!\n
                    trunc = mps.grow_r2l_sweep(mpo[site], left_[site], right_[site+1],
                                               site, alpha, cutoff, Mmax);
                } else {
                    block_matrix<Matrix, SymmGroup> t = mps[site].normalize_right(DefaultSolver());
                    if (site > 0)
                        mps[site-1].multiply_from_right(t);
                }
                
                
                Storage::drop(right_[site]); // right_[site] is outdated
                this->boundary_right_step(mpo, site); // creating right_[site]
            }
            
        	Storage::evict(left_[site]); // move to out of core currently used boundary
        	Storage::evict(right_[site+1]); // move to out of core currently used boundary

            iteration_results_["BondDimension"]   << trunc.bond_dimension;
            iteration_results_["TruncatedWeight"] << trunc.truncated_weight;
            iteration_results_["SmallestEV"]      << trunc.smallest_ev;
            
            gettimeofday(&sweep_then, NULL);
            double elapsed = sweep_then.tv_sec-sweep_now.tv_sec + 1e-6 * (sweep_then.tv_usec-sweep_now.tv_usec);
            maquis::cout << "Sweep has been running for " << elapsed << " seconds." << std::endl;
            
            if (stop_callback())
                throw dmrg::time_limit(sweep, _site+1);
        }
    }
    
};

#endif

