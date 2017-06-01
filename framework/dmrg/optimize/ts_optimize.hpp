/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Bela Bauer <bauerb@phys.ethz.ch> 
 *	                          Sebastian Keller <sebkelle@phys.ethz.ch>
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

#ifndef TS_OPTIMIZE_H
#define TS_OPTIMIZE_H

#include "dmrg/optimize/optimize.h"

#include "dmrg/mp_tensors/twositetensor.h"
#include "dmrg/mp_tensors/mpo_ops.h"

#include <boost/tuple/tuple.hpp>

template<class Matrix, class SymmGroup, class Storage>
class ts_optimize : public optimizer_base<Matrix, SymmGroup, Storage>
{
public:

    typedef optimizer_base<Matrix, SymmGroup, Storage> base;
    using base::mpo;
    using base::mps;
    using base::left_;
    using base::right_;
    using base::parms;
    using base::iteration_results_;
    using base::stop_callback;
    // Constructor
    ts_optimize(MPS<Matrix, SymmGroup> & mps_,
                MPO<Matrix, SymmGroup> const & mpo_,
                BaseParameters & parms_,
                boost::function<bool ()> stop_callback_,
                int initial_site_ = 0)
    : base(mps_, mpo_, parms_, stop_callback_, to_site(mps_.length(), initial_site_))
    , initial_site((initial_site_ < 0) ? 0 : initial_site_)
    {
        parallel::guard::serial guard;
        make_ts_cache_mpo(mpo, ts_cache_mpo, mps);
    }
    // Inline function to convert from 2L to L the index of the site
    inline int to_site(const int L, const int i) const
    {
        if (i < 0) return 0;
        /// i, or (L-1) - (i - (L-1))
        return (i < L-1) ? i : 2*L - 2 - i;
    }
    // Function to perform a sweep
    void sweep(int sweep, OptimizeDirection d = Both)
    {
        boost::chrono::high_resolution_clock::time_point sweep_now = boost::chrono::high_resolution_clock::now();

        iteration_results_.clear();
        
        std::size_t L = mps.length();
        parallel::scheduler_balanced scheduler_mps(L);
        // Definition of the initial site
        int _site = 0, site = 0;
        if (initial_site != -1) {
            _site = initial_site;
            site = to_site(L, _site);
        }
        
        if (_site < L-1) {
            Storage::prefetch(left_[site]);
            Storage::prefetch(right_[site+2]);
        } else {
            Storage::prefetch(left_[site-1]);
            Storage::prefetch(right_[site+1]);
        }
        
        for (; _site < 2*L-2; ++_site) {
	/* (0,1), (1,2), ... , (L-1,L), (L-1,L), (L-2, L-1), ... , (0,1)
	    | |                        |
       site 1                      |
	      |         left to right  | right to left, lr = -1
	      site 2                   |                               */

            int lr, site1, site2;
            if (_site < L-1) {
                site = to_site(L, _site);
                lr = 1;
        		site1 = site;
        		site2 = site+1;
                ts_cache_mpo[site1].placement_l = mpo[site1].placement_l;
                ts_cache_mpo[site1].placement_r = parallel::get_right_placement(ts_cache_mpo[site1], mpo[site1].placement_l, mpo[site2].placement_r);
            } else {
                site = to_site(L, _site);
                lr = -1;
        		site1 = site-1;
        		site2 = site;
                ts_cache_mpo[site1].placement_l = parallel::get_left_placement(ts_cache_mpo[site1], mpo[site1].placement_l, mpo[site2].placement_r);
                ts_cache_mpo[site1].placement_r = mpo[site2].placement_r;
            }

            //if (lr == +1) mps.canonize(site1);
            //else          mps.canonize(site2);
            
    	    maquis::cout << std::endl;
            maquis::cout << "Sweep " << sweep << ", optimizing sites " << site1 << " and " << site2 << std::endl;

            // MD: some changes needed to re-enable it.
//            if (parms.template get<bool>("beta_mode")) {
//                if (sweep == 0 && lr == 1) {
//                    mpo = zero_after(mpo_orig, 0);
//                    if (site == 0)
//                        this->init_left_right(mpo, 0);
//                } else if (sweep == 0 && lr == -1 && site == L-1) {
//                    mpo = mpo_orig;
//                }
//            }
        
            if (_site != L-1)
            { 
                Storage::fetch(left_[site1]);
                Storage::fetch(right_[site2+1]);
            }

            if (lr == +1) {
                if (site2+2 < right_.size()){
                    Storage::prefetch(right_[site2+2]);
                }
            } else {
                if (site1 > 0){
                    Storage::prefetch(left_[site1-1]);
                }
            }


            boost::chrono::high_resolution_clock::time_point now, then;
            
    	    // Create TwoSite objects
    	    TwoSiteTensor<Matrix, SymmGroup> tst(mps[site1], mps[site2]);
    	    MPSTensor<Matrix, SymmGroup> twin_mps = tst.make_mps();
            tst.clear();
            SiteProblem<Matrix, SymmGroup> sp(left_[site1], right_[site2+1], ts_cache_mpo[site1]);
            
            /// Compute orthogonal vectors
            std::vector<MPSTensor<Matrix, SymmGroup> > ortho_vecs(base::northo);
            for (int n = 0; n < base::northo; ++n) {
                TwoSiteTensor<Matrix, SymmGroup> ts_ortho(base::ortho_mps[n][site1], base::ortho_mps[n][site2]);
                ortho_vecs[n] = contraction::site_ortho_boundaries(twin_mps, ts_ortho.make_mps(),
                                                                    base::ortho_left_[n][site1], base::ortho_right_[n][site2+1]);
            }

            std::pair<double, MPSTensor<Matrix, SymmGroup> > res;


            if (d == Both ||
                (d == LeftOnly && lr == -1) ||
                (d == RightOnly && lr == +1))
            {
                if (parms["eigensolver"] == std::string("IETL")) {
            	    BEGIN_TIMING("IETL")
                    res = solve_ietl_lanczos(sp, twin_mps, parms);
            	    END_TIMING("IETL")
                } else if (parms["eigensolver"] == std::string("IETL_JCD")) {
            	    BEGIN_TIMING("JCD")
                    // solve_ietl_jcd sta dentro ietl_jacobi_davidson
                    res = solve_ietl_jcd(sp, twin_mps, parms, ortho_vecs);
            	    END_TIMING("JCD")
                } else if (parms["eigensolver"] == std::string("IETL_DAVIDSON")) {
            	    BEGIN_TIMING("DAVIDSON")
                    // solve_ietl_davidson sta dentro ietd_davidson
                    res = solve_ietl_davidson(sp, twin_mps, parms, ortho_vecs);
            	    END_TIMING("DAVIDSON")
                } else if (parms["eigensolver"] == std::string("IETL_MODIFIED_DAVIDSON")) {
                    BEGIN_TIMING("MODIFIED_DAVIDSON")
                    res = solve_ietl_davidson_modified(sp, twin_mps, parms, parms["ietl_moddav_omega"], ortho_vecs);
                    END_TIMING("MODIFIED_DAVIDSON")
                } else if (parms["eigensolver"] == std::string("IETL_MODIFIED_JCD")) {
                    BEGIN_TIMING("MODIFIED_JCD")
                    res = solve_ietl_jcd_modified(sp, twin_mps, parms, parms["ietl_modjcd_omega"], ortho_vecs);
                    END_TIMING("MODIFIED_JCD")
                } else {
                    throw std::runtime_error("I don't know this eigensolver.");
                }

        		tst << res.second;
                res.second.clear();
            }
            twin_mps.clear();


#ifndef NDEBUG
            // Caution: this is an O(L) operation, so it really should be done only in debug mode
            for (int n = 0; n < base::northo; ++n)
                maquis::cout << "MPS overlap: " << overlap(mps, base::ortho_mps[n]) << std::endl;
#endif

            {
                int prec = maquis::cout.precision();
                maquis::cout.precision(15);
                maquis::cout << "Energy " << lr << " " << res.first + mpo.getCoreEnergy() << std::endl;
                maquis::cout.precision(prec);
            }
            iteration_results_["Energy"] << res.first + mpo.getCoreEnergy();
            
            
            double alpha;
            int ngs = parms["ngrowsweeps"], nms = parms["nmainsweeps"];
            if (sweep < ngs)
                alpha = parms["alpha_initial"];
            else if (sweep < ngs + nms)
                alpha = parms["alpha_main"];
            else
                alpha = parms["alpha_final"];

            double cutoff = this->get_cutoff(sweep);
            std::size_t Mmax = this->get_Mmax(sweep);
            truncation_results trunc;
            
    	    if (lr == +1)
    	    {
        		// Write back result from optimization
                BEGIN_TIMING("TRUNC")
                if (parms["twosite_truncation"] == "svd")
                    boost::tie(mps[site1], mps[site2], trunc) = tst.split_mps_l2r(Mmax, cutoff);
                else
                    boost::tie(mps[site1], mps[site2], trunc) = tst.predict_split_l2r(Mmax, cutoff, alpha, left_[site1], mpo[site1]);
                END_TIMING("TRUNC")
                tst.clear();


        		block_matrix<Matrix, SymmGroup> t;
		
        		//t = mps[site1].normalize_left(DefaultSolver());
        		//mps[site2].multiply_from_left(t);
        		//mps[site2].divide_by_scalar(mps[site2].scalar_norm());	

        		t = mps[site2].normalize_left(DefaultSolver());
                // MD: DEBUGGING OUTPUT
                maquis::cout << "Propagating t with norm " << t.norm() << std::endl;
        		if (site2 < L-1) mps[site2+1].multiply_from_left(t);

                if (site1 != L-2)
                    Storage::drop(right_[site2+1]);

                this->boundary_left_step(mpo, site1); // creating left_[site2]

                if (site1 != L-2){ 
                    if(site1 != 0){
                        #ifdef USE_AMBIENT
                        std::vector<int> placement_l = parallel::get_left_placement(ts_cache_mpo[site1], mpo[site1].placement_l, mpo[site2].placement_r);
                        parallel::scheduler_permute scheduler(placement_l, parallel::groups_granularity);
                        for(size_t b = 0; b < left_[site1].aux_dim(); ++b){
                            parallel::guard group(scheduler(b), parallel::groups_granularity);
                            storage::migrate(left_[site1][b], parallel::scheduler_size_indexed(left_[site1][b]));
                        }
                        parallel::sync();
                        #endif
                    }
                    Storage::evict(mps[site1]);
                    Storage::evict(left_[site1]);
                }
                { parallel::guard proc(scheduler_mps(site1)); storage::migrate(mps[site1]); }
                { parallel::guard proc(scheduler_mps(site2)); storage::migrate(mps[site2]); }
    	    }
    	    if (lr == -1){
        		// Write back result from optimization
                BEGIN_TIMING("TRUNC")
                if (parms["twosite_truncation"] == "svd")
                    boost::tie(mps[site1], mps[site2], trunc) = tst.split_mps_r2l(Mmax, cutoff);
                else
                    boost::tie(mps[site1], mps[site2], trunc) = tst.predict_split_r2l(Mmax, cutoff, alpha, right_[site2+1], mpo[site2]);
                END_TIMING("TRUNC")
                tst.clear();


        		block_matrix<Matrix, SymmGroup> t;

        		//t = mps[site2].normalize_right(DefaultSolver());
        		//mps[site1].multiply_from_right(t);
        		//mps[site1].divide_by_scalar(mps[site1].scalar_norm());	

        		t = mps[site1].normalize_right(DefaultSolver());
                // MD: DEBUGGING OUTPUT
                maquis::cout << "Propagating t with norm " << t.norm() << std::endl;
        		if (site1 > 0) mps[site1-1].multiply_from_right(t);

                if(site1 != 0)
                    Storage::drop(left_[site1]);

                this->boundary_right_step(mpo, site2); // creating right_[site2]

                if(site1 != 0){
                    if(site1 != L-2){
                        #ifdef USE_AMBIENT
                        std::vector<int> placement_r = parallel::get_right_placement(ts_cache_mpo[site1], mpo[site1].placement_l, mpo[site2].placement_r);
                        parallel::scheduler_permute scheduler(placement_r, parallel::groups_granularity);
                        for(size_t b = 0; b < right_[site2+1].aux_dim(); ++b){
                            parallel::guard group(scheduler(b), parallel::groups_granularity);
                            storage::migrate(right_[site2+1][b], parallel::scheduler_size_indexed(right_[site2+1][b]));
                        }
                        parallel::sync();
                        #endif
                    }
                    Storage::evict(mps[site2]);
                    Storage::evict(right_[site2+1]); 
                }
                { parallel::guard proc(scheduler_mps(site1)); storage::migrate(mps[site1]); }
                { parallel::guard proc(scheduler_mps(site2)); storage::migrate(mps[site2]); }
    	    }
            
            iteration_results_["BondDimension"]     << trunc.bond_dimension;
            iteration_results_["TruncatedWeight"]   << trunc.truncated_weight;
            iteration_results_["TruncatedFraction"] << trunc.truncated_fraction;
            iteration_results_["SmallestEV"]        << trunc.smallest_ev;
            
            parallel::meminfo();
            
            boost::chrono::high_resolution_clock::time_point sweep_then = boost::chrono::high_resolution_clock::now();
            double elapsed = boost::chrono::duration<double>(sweep_then - sweep_now).count();
            maquis::cout << "Sweep has been running for " << elapsed << " seconds." << std::endl;
            
            if (stop_callback())
                throw dmrg::time_limit(sweep, _site+1);

    	} // for sites
        initial_site = -1;
    } // sweep

private:
    int initial_site;
    MPO<Matrix, SymmGroup> ts_cache_mpo;
};

#endif
