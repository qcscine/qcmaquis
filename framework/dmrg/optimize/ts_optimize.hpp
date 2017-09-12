/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Bela Bauer <bauerb@phys.ethz.ch> 
 *	                          Sebastian Keller <sebkelle@phys.ethz.ch>
 *	             2017 by Alberto Baiardi <alberto.baiardi@sns.it>
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
#include "dmrg/optimize/partial_overlap.h"
#include "dmrg/optimize/vectorset.h"
#include <boost/tuple/tuple.hpp>

//
// TS_OPTIMIZE CLASS
// -----------------
// Inherited from the virtual optimizer_base class
template<class Matrix, class SymmGroup, class Storage>
class ts_optimize : public optimizer_base<Matrix, SymmGroup, Storage>
{
public:
    typedef optimizer_base<Matrix, SymmGroup, Storage> base;
    typedef typename partial_overlap<Matrix, SymmGroup>::partial_overlap partial_overlap;
    using base::bound_vec_pnt_ ;
    using base::do_root_homing_ ;
    using base::iteration_results_;
    using base::left_;
    using base::left_sa_;
    using base::mpo;
    using base::mps;
    using base::mps2follow ;
    using base::mps_vector;
    using base::n_root_ ;
    using base::order ;
    using base::parms;
    using base::poverlap_vec_ ;
    using base::right_;
    using base::right_sa_ ;
    using base::root_homing_type_ ;
    using base::stop_callback ;
    using base::vec_sa_left_ ;
    using base::vec_sa_right_ ;
    // Constructor
    ts_optimize(MPS<Matrix, SymmGroup> & mps_,
                std::vector< MPS<Matrix, SymmGroup> > & mps_sa_,
                MPO<Matrix, SymmGroup> const & mpo_,
                BaseParameters & parms_,
                boost::function<bool ()> stop_callback_,
                int initial_site_ = 0)
    : base(mps_, mps_sa_, mpo_, parms_, stop_callback_, to_site(mps_.length(), initial_site_))
    , initial_site((initial_site_ < 0) ? 0 : initial_site_)
    {
        parallel::guard::serial guard;
        make_ts_cache_mpo(mpo, ts_cache_mpo, mps);
    }
    // Inline function to convert from 2L to L the index of the site
    inline int to_site(const int L, const int i) const
    {
        if (i < 0) return 0;
        // i, or (L-1) - (i - (L-1))
        return (i < L-1) ? i : 2*L - 2 - i;
    }
    // Function to perform a sweep
    void sweep(int sweep, OptimizeDirection d = Both)
    {
        // Initialization
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
        std::vector< std::pair<float,int> > sorter ;
        sorter.resize(n_root_) ;
        if (_site < L-1) {
            Storage::prefetch(left_[site]);
            Storage::prefetch(right_[site+2]);
        } else {
            Storage::prefetch(left_[site-1]);
            Storage::prefetch(right_[site+1]);
        }
        // +------------------------------+
        //  MAIN LOOP - SWEEP OPTIMIZATION
        // +------------------------------+
        for (; _site < 2*L-2; ++_site) {
            //
            double i ;
            // -- GENERATES THE TWO SITE MPO --
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
    	    maquis::cout << std::endl;
            maquis::cout << "Sweep " << sweep << ", optimizing sites " << site1 << " and " << site2 << std::endl;
            // -- STORAGE MANAGEMENT --
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
    	    // Create TwoSite objects. For SA calculations, creates a vector
    	    TwoSiteTensor<Matrix, SymmGroup> tst(mps[site1], mps[site2]);
            MPSTensor<Matrix, SymmGroup> twin_mps = tst.make_mps();
            tst.clear();
            //
            std::vector< MPSTensor<Matrix, SymmGroup> > tst_vec ;
            std::vector< TwoSiteTensor<Matrix, SymmGroup> > two_vec ;
            for (int i = 0 ; i < n_root_ ; i++){
                TwoSiteTensor<Matrix, SymmGroup> tst_tmp(mps_vector[i][site1],mps_vector[i][site2]) ;
                two_vec.push_back(tst_tmp) ;
                twin_mps = tst_tmp.make_mps() ;
                tst_vec.push_back(twin_mps) ;
                tst_tmp.clear();
            }
            SiteProblem<Matrix, SymmGroup> sp(left_sa_, right_sa_, ts_cache_mpo[site1], site1, site2+1, bound_vec_pnt_);
            VectorSet<Matrix,SymmGroup> vector_set(tst_vec) ;
            // Compute orthogonal vectors
            std::vector<MPSTensor<Matrix, SymmGroup> > ortho_vecs(base::northo);
            for (int n = 0; n < base::northo; ++n) {
                TwoSiteTensor<Matrix, SymmGroup> ts_ortho(base::ortho_mps[n][site1], base::ortho_mps[n][site2]);
                ortho_vecs[n] = contraction::site_ortho_boundaries(twin_mps,
                                                                   ts_ortho.make_mps(),
                                                                   base::ortho_left_[n][site1],
                                                                   base::ortho_right_[n][site2+1] );
            }
            // +-----------------------------+
            //  MAIN PART: performs the sweep
            // +-----------------------------+
            std::vector<std::pair<double, MPSTensor<Matrix, SymmGroup> > > res;
            if (d == Both ||
                (d == LeftOnly && lr == -1) ||
                (d == RightOnly && lr == +1))
            {
                if (parms["eigensolver"] == std::string("IETL")) {
            	    BEGIN_TIMING("IETL")
                    //res = solve_ietl_lanczos(sp, twin_mps, parms);
            	    END_TIMING("IETL")
                } else if (parms["eigensolver"] == std::string("IETL_JCD")) {
            	    BEGIN_TIMING("JCD")
                    res = solve_ietl_jcd(sp, vector_set, parms, poverlap_vec_, 2, site1, site2, root_homing_type_,
                                         vec_sa_left_, vec_sa_right_, order, ortho_vecs);
            	    END_TIMING("JCD")
                } else if (parms["eigensolver"] == std::string("IETL_DAVIDSON")) {
                    BEGIN_TIMING("DAVIDSON")
                    //res = solve_ietl_davidson(sp, vector_set, parms, poverlap_vec_, 2, site1, root_homing_type_, ortho_vecs, site2);
                    END_TIMING("DAVIDSON")
                } else {
                    throw std::runtime_error("I don't know this eigensolver.");
                }
        		tst << res[0].second;
                if (n_root_ > 0) {
                    two_vec[0] << res[0].second ;
                    res[0].second.clear();
                    for (size_t k = 1; k < n_root_; k++) {
                        two_vec[k] << res[k].second;
                        res[k].second.clear();
                    }
                }
            }
            // +---------------------+
            //  Collection of results
            // +---------------------+
            twin_mps.clear();
            {
                int prec = maquis::cout.precision();
                maquis::cout.precision(15);
                for (size_t k = 0 ; k < n_root_ ; k++)
                    maquis::cout << " Converged energy - state " << k << " = " << res[k].first + mpo.getCoreEnergy() << std::endl;
                maquis::cout.precision(prec);
            }
            iteration_results_["Energy"] << res[0].first + mpo.getCoreEnergy();
            // +---------------------+
            //  Truncation of the MPS
            // +---------------------+
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
                for (std::size_t i = 0; i < n_root_; i++) {
                    BEGIN_TIMING("TRUNC")
                    if (parms["twosite_truncation"] == "svd")
                        boost::tie(mps_vector[i][site1], mps_vector[i][site2], trunc) = two_vec[i].split_mps_l2r(Mmax, cutoff);
                    else
                        boost::tie(mps_vector[i][site1], mps_vector[i][site2], trunc) = two_vec[i].predict_split_l2r(Mmax, cutoff, alpha, left_[site1], mpo[site1]);
                    END_TIMING("TRUNC")
                    two_vec[i].clear();
        		    block_matrix<Matrix, SymmGroup> t;
        		    t = mps_vector[i][site2].normalize_left(DefaultSolver());
                    // MD: DEBUGGING OUTPUT
                    maquis::cout << "Propagating t with norm " << t.norm() << std::endl;
        		    if (site2 < L-1) mps_vector[i][site2+1].multiply_from_left(t);
                }
                if (site1 != L-2)
                    Storage::drop(right_[site2+1]);
                this->boundary_left_step(mpo, site1); // creating left_[site2]
                if (site1 != L-2){
                    Storage::evict(mps[site1]);
                    Storage::evict(left_[site1]);
                }
                { parallel::guard proc(scheduler_mps(site1)); storage::migrate(mps[site1]); }
                { parallel::guard proc(scheduler_mps(site2)); storage::migrate(mps[site2]); }
    	    }
    	    if (lr == -1){
        		// Write back result from optimization
                for (std::size_t i = 0; i < n_root_; i++) {
                    BEGIN_TIMING("TRUNC")
                    if (parms["twosite_truncation"] == "svd")
                        boost::tie(mps_vector[i][site1], mps_vector[i][site2], trunc) = two_vec[i].split_mps_r2l(Mmax, cutoff);
                    else
                        boost::tie(mps_vector[i][site1], mps_vector[i][site2], trunc) = two_vec[i].predict_split_r2l(Mmax, cutoff, alpha, right_[site2+1], mpo[site2]);
                    END_TIMING("TRUNC")
                    two_vec[i].clear();
        		    block_matrix<Matrix, SymmGroup> t;
        		    t = mps_vector[i][site1].normalize_right(DefaultSolver());
                    // MD: DEBUGGING OUTPUT
                    maquis::cout << "Propagating t with norm " << t.norm() << std::endl;
        		    if (site1 > 0) mps_vector[i][site1-1].multiply_from_right(t);
                }
                if(site1 != 0)
                    Storage::drop(left_[site1]);
                this->boundary_right_step(mpo, site2); // creating right_[site2]
                if(site1 != 0){
                    Storage::evict(mps[site2]);
                    Storage::evict(right_[site2+1]); 
                }
                { parallel::guard proc(scheduler_mps(site1)); storage::migrate(mps[site1]); }
                { parallel::guard proc(scheduler_mps(site2)); storage::migrate(mps[site2]); }
    	    }
            // +------------+
            //  FINALIZATION
            // +------------+
            for (size_t i = 0; i < n_root_; i++) {
                sorter[i].first  = res[i].first ;
                sorter[i].second = i ;
            }
            std::sort(sorter.begin(), sorter.end()) ;
            this->update_order(sorter) ;
            if (root_homing_type_ == 1)
                for (size_t k = 0 ; k < n_root_ ; k++)
                    poverlap_vec_[k].update(mps_vector[k], site, lr) ;
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
