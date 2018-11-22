/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Bela Bauer <bauerb@phys.ethz.ch>
 *	                          Sebastian Keller <sebkelle@phys.ethz.ch>
 *	             2017-2018 by Alberto Baiardi <alberto.baiardi@sns.it>
 *               2018 by Leon Freitag <lefreita@ethz.ch>
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
#include "dmrg/optimize/POverlap/partial_overlap.h"
#include "dmrg/optimize/vectorset.h"
#include <boost/tuple/tuple.hpp>

// ===================
//  TS_OPTIMIZE CLASS
// ===================
//
template<class Matrix, class SymmGroup, class Storage>
class ts_optimize : public optimizer_base<Matrix, SymmGroup, Storage>
{
public:
    typedef typename Matrix::value_type value_type;
    typedef optimizer_base<Matrix, SymmGroup, Storage> base;
    typedef typename partial_overlap<Matrix, SymmGroup>::partial_overlap partial_overlap;
    using base::boundaries_database_ ;
    using base::chebyshev_shift_ ;
    using base::correction_equation ;
    using base::do_chebyshev_ ;
    using base::do_H_squared_ ;
    using base::do_root_homing_ ;
    using base::do_shiftandinvert_ ;
    using base::finalizer_ ;
   	using base::is_folded_ ;
    using base::iteration_results_;
    using base::L_ ;
    using base::left_sa_;
    using base::left_squared_sa_ ;
    using base::micro_optimizer ;
    using base::mpo;
    using base::mpo_squared;
    using base::mps_guess ;
    using base::mps_vector ;
    using base::n_bound_ ;
    using base::n_root_ ;
    using base::omega_shift_ ;
    using base::omega_vec ;
    using base::order ;
    using base::orthogonalizer_ ;
    using base::parms ;
    using base::poverlap_vec_ ;
    using base::right_sa_ ;
    using base::right_squared_sa_ ;
    using base::root_homing_type_ ;
    using base::sa_alg_ ;
    using base::sorter_ ;
    using base::stop_callback ;
    using base::reshuffle_variance_ ;
    using base::track_variance_ ;
    using base::update_each_omega ;
    using base::update_omega ;
    using base::update_order ;
    using base::vec_sa_left_ ;
    using base::vec_sa_right_ ;
    // Constructor
    ts_optimize(std::vector< MPS<Matrix, SymmGroup> > & mps_vector_,
                std::vector< MPSTensor<Matrix, SymmGroup> > & mps_guess_,
                std::vector< MPS<Matrix, SymmGroup> > & mps_partial_overlap_,
                MPO<Matrix, SymmGroup> const & mpo_,
                MPO<Matrix, SymmGroup> const & mpo_squared_,
                BaseParameters & parms_,
                boost::function<bool ()> stop_callback_,
                int initial_site_ = 0)
    : base(mps_vector_, mps_guess_, mps_partial_overlap_, mpo_, mpo_squared_, parms_, stop_callback_,
           to_site(mps_vector_[0].length(), initial_site_))
    , initial_site((initial_site_ < 0) ? 0 : initial_site_)
    {
        parallel::guard::serial guard;
        make_ts_cache_mpo(mpo, ts_cache_mpo, mps_vector[0]) ;
        make_ts_cache_mpo(mpo_squared, ts_cache_mpo_squared, mps_vector[0]) ;
        iteration_results_.resize(n_root_);
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

        for (auto&& it : iteration_results_)
            it.clear();
        parallel::scheduler_balanced scheduler_mps(L_);
        // Definition of the initial site
        int _site = 0;
        int site = 0 ;
        if (initial_site != -1) {
            _site = initial_site;
            site = to_site(L_, _site);
        }
        // +------------------------------+
        //  MAIN LOOP - SWEEP OPTIMIZATION
        // +------------------------------+
        this->update_parameters(sweep) ;
        for (; _site < 2*L_-2; ++_site) {
            //
            double i ;
            // Debug printing
            if (do_shiftandinvert_) {
                std::cout << " -- VALUES OF OMEGA -- " << std::endl ;
                for (size_t idx = 0; idx < n_root_; idx++)
                    std::cout << omega_vec[idx] << std::endl ;
            }
            // -- GENERATES THE TWO SITE MPO --
            int lr, site1, site2;
            if (_site < L_-1) {
                site = to_site(L_, _site);
                lr = 1;
                site1 = site ;
                site2 = site+1 ;
                ts_cache_mpo[site1].placement_l = mpo[site1].placement_l;
                ts_cache_mpo[site1].placement_r = parallel::get_right_placement(ts_cache_mpo[site1], mpo[site1].placement_l, mpo[site2].placement_r);
            } else {
                site = to_site(L_, _site);
                lr = -1;
                site1 = site-1 ;
                site2 = site ;
                ts_cache_mpo[site1].placement_l = parallel::get_left_placement(ts_cache_mpo[site1], mpo[site1].placement_l, mpo[site2].placement_r);
                ts_cache_mpo[site1].placement_r = mpo[site2].placement_r;
            }
            // Printing
            print_header(sweep, site1, site2, lr) ;
            boost::chrono::high_resolution_clock::time_point now, then;
       	    // Create TwoSite objects. For SA calculations, creates a vector
            MPSTensor<Matrix, SymmGroup> twin_mps ;
            std::vector<MPSTensor<Matrix, SymmGroup> > tst_vec ;
            std::vector<TwoSiteTensor<Matrix, SymmGroup> > two_vec ;
            BEGIN_TIMING("TWO SITE TENSOR DEFINITION")
            for (int i = 0 ; i < n_root_ ; i++){
                TwoSiteTensor<Matrix, SymmGroup> tst_tmp(mps_vector[i][site1],mps_vector[i][site2]) ;
                two_vec.push_back(tst_tmp) ;
                twin_mps = tst_tmp.make_mps() ;
                tst_vec.push_back(twin_mps) ;
                tst_tmp.clear();
            }
            END_TIMING("TWO SITE TENSOR DEFINITION")
            SiteProblem<Matrix, SymmGroup> sp(tst_vec, ts_cache_mpo[site1], ts_cache_mpo_squared[site1], site1, site2+1,
                                              boundaries_database_, do_H_squared_);
            VectorSet<Matrix,SymmGroup> vector_set(tst_vec) ;
            // Compute orthogonal vectors
            std::vector<MPSTensor<Matrix, SymmGroup> > ortho_vecs(base::northo);
            for (int n = 0; n < base::northo; ++n) {
                TwoSiteTensor<Matrix, SymmGroup> ts_ortho(base::ortho_mps[n][site1], base::ortho_mps[n][site2]);
                ortho_vecs[n] = contraction::site_ortho_boundaries(tst_vec[0],
                                                                   ts_ortho.make_mps(),
                                                                   base::ortho_left_[n][site1],
                                                                   base::ortho_right_[n][site2+1] );
            }
            // +-----------------------------+
            //  MAIN PART: performs the sweep
            // +-----------------------------+
            std::vector<std::pair<typename maquis::traits::real_type<value_type>::type, MPSTensor<Matrix, SymmGroup> > > res;
            for (std::size_t idx = 0; idx < poverlap_vec_.size(); idx++)
                poverlap_vec_[idx].prepare(tst_vec[idx], site1, site2) ;
            if (d == Both ||
                (d == LeftOnly && lr == -1) ||
                (d == RightOnly && lr == +1))
            {
                BEGIN_TIMING("LOCAL DIAGONALIZATION")
                res = solve_ietl_jcd(sp, vector_set, correction_equation, micro_optimizer, finalizer_,
                                     orthogonalizer_, parms, poverlap_vec_, site1, site2, root_homing_type_, do_root_homing_,
                                     do_shiftandinvert_, do_chebyshev_, chebyshev_shift_, do_H_squared_, reshuffle_variance_,
                                     track_variance_, is_folded_, vec_sa_left_, vec_sa_right_, order, boundaries_database_,
                                     sa_alg_, omega_vec, ortho_vecs);
                END_TIMING("LOCAL DIAGONALIZATION")
                // Correct the energies
                if (sa_alg_ == -1)
                    for (size_t k = 0; k < n_root_; k++)
                       res[k].first = ietl::get_energy(sp, res[k].second, k, false) ;
                // Collects the results
                two_vec[0] << res[0].second ;
                res[0].second.clear();
                for (size_t k = 1; k < n_root_; k++) {
                    two_vec[k] << res[k].second;
                    res[k].second.clear();
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
                {
                    maquis::cout << " Converged energy - state " << k << " = " << res[k].first + mpo.getCoreEnergy() << std::endl;
                    iteration_results_[k]["Energy"] << res[k].first + mpo.getCoreEnergy();
                }
                maquis::cout.precision(prec);
            }
            // +------------------------------------+
            //  Setting up parameters for truncation
            // +------------------------------------+
            double alpha;
            int ngs = parms["ngrowsweeps"], nms = parms["nmainsweeps"];
            if (sweep < ngs)
                alpha = parms["alpha_initial"];
            else if (sweep < ngs + nms)
                alpha = parms["alpha_main"];
            else
                alpha = parms["alpha_final"];
            double cutoff = this->get_cutoff(sweep) ;
            std::size_t Mmax = this->get_Mmax(sweep) ;
            // WARNING & TODO
            // for sa_alg == -3 only trunc[0] is used!
            // sa_alg=-1 WIP
            std::vector<truncation_results> trunc(n_root_) ;
            // +--------------------------------------+
            //  Truncation of the TwoSiteTensor object
            // +--------------------------------------+
            BEGIN_TIMING("MPS TRUNCATION")
            // -------------
            // Forward sweep
            // -------------
            if (lr == +1) {
                // State-averaged optimisation
                if (sa_alg_ == -3 || sa_alg_ == -1) {
                    //
                    std::vector< MPSTensor<Matrix,SymmGroup> > vec_mps ;
                    MPSTensor<Matrix, SymmGroup> tst_mps ;
                    if (sa_alg_ == -3) // New state-average model (truncation based on a super-DM)
                    {
                        std::tie(tst_mps, vec_mps, trunc[0]) = split_mps_l2r_vector(two_vec, Mmax, cutoff);
                        int Mval = trunc[0].bond_dimension ;
                        std::cout << "Bond dimension " << Mval << std::endl ;
                    }
                    else //sa_alg_ == -1: Standard, state-average model (truncation based on an average DM)
                        std::tie(tst_mps, vec_mps, trunc[0]) = predict_split_l2r_avg(two_vec, Mmax, cutoff,
                             alpha, (*(boundaries_database_.get_boundaries_left(0, false)))[site1], mpo[site1]);

                    // Final update of the MPS
                    for (std::size_t idx = 0; idx < n_root_; idx++) {
                        two_vec[idx].clear() ;
                        block_matrix<Matrix, SymmGroup> t ;
                        (*(boundaries_database_.get_mps(idx)))[site1] = tst_mps ;
                        (*(boundaries_database_.get_mps(idx)))[site2] = vec_mps[idx] ;
                        t = (*(boundaries_database_.get_mps(idx)))[site2].normalize_left(DefaultSolver()) ;
                        if (site2 < L_-1)
                            (*(boundaries_database_.get_mps(idx)))[site2+1].multiply_from_left(t) ;
                    }
                }
//                else if (sa_alg_ == -1) {
//                     // old sa_alg_ == -1 code that never worked correctly. if the above works, it should be deleted

//                     if (parms["twosite_truncation"] != "svd") throw std::runtime_error("twosite_truncation != svd with SA optimization not supported yet!");

//                     // Construct the two-site tensor averaged over all states
//                     TwoSiteTensor<Matrix, SymmGroup> avg_tst = std::accumulate(std::next(two_vec.begin()), two_vec.end(), *two_vec.begin());
//                     avg_tst /= n_root_;

//                     // Truncation results for the average two-site tensor
//                     truncation_results avg_truncation;

//                     // sa_alg=-1 can be carried out in several variants:
//                     // - the average TwoSiteTensor is split into two sites with the SVD [ P = U.S.V ]
//                     // - U becomes the new MPSTensor at site 1 and S.V at site 2 for _all states_ [ for the left sweep, for the right sweep it is U.S and V, accordingly ]
//                     // We call this variant "USV" and it can be enabled by uncommenting the following define:
// // #define _SAALG_VAR_USV_
//                     // Additionally, we may perform a state-specific TwoSiteTensor splitting and take either U from the state-specific TST [ or V for the right sweep ]
//                     // leading to the "SV" variant, enabled by uncommenting the followind define:
// #define _SAALG_VAR_SV_
//                     // or both U and V from the state-specific TST, where only S remains from the state-average TST:
// // #define _SAALG_VAR_S_
//                     // We enforce the truncation with the same number of retained eigenvalues in each state-specific truncated SVD
//                     // by taking the number of eigenvalues retained in the SVD of the average TST
//                     // The following table illustrates from which splitting the matrices are used in the optimization (A -- state-averaged, S -- state-specific)
//                     //
//                     //         Forward sweep           Backward sweep
//                     //
//                     //      Variant | U | S | V      Variant | U | S | V
//                     //     ----------------------   ---------------------
//                     //        USV   | A | A | A        USV   | A | A | A
//                     //        SV    | S | A | A        SV    | A | A | S
//                     //        S     | S | A | S        S     | S | A | S
// #if ( defined(_SAALG_VAR_S_) && defined(_SAALG_VAR_SV_) ) || ( defined(_SAALG_VAR_SV_) && defined(_SAALG_VAR_USV_) ) || ( defined(_SAALG_VAR_S_) && defined(_SAALG_VAR_USV_) )
// #error Only one of _SAALG_VAR_S_, _SAALG_VAR_SV_ or _SAALG_VAR_USV_ may be defined!
// #endif

// #ifdef _SAALG_VAR_S_
//                     // sa_alg=-1 variant S ////////////
//                     // Perform truncation of the average TwoSiteTensor and obtain the S (diagonal) matrix

//                     typename TwoSiteTensor<Matrix, SymmGroup>::block_diag_matrix s_avg;
//                     std::tie(s_avg, avg_truncation) = avg_tst.get_S(Mmax, cutoff);
// #endif

// #if defined(_SAALG_VAR_SV_) || defined (_SAALG_VAR_USV_)
//                     // sa_alg=-1 variants SV and USV
//                     // Simple splitting of the average MPS tensor, to be applied to all states
//                     // Get U,S,V from the average
//                     MPSTensor<Matrix, SymmGroup> split_site1, split_site2;
//                     std::tie(split_site1, split_site2, avg_truncation) = avg_tst.split_mps_l2r(Mmax, cutoff);

// #endif
//                     for (size_t idx = 0; idx < n_root_; idx++) {
// #ifdef _SAALG_VAR_S_
//                         // sa_alg=-1 variant S
//                         // for an external S split_mps_l2r is overloaded and doesn't return truncation_results
//                         std::tie(mps_vector[idx][site1], mps_vector[idx][site2]) = two_vec[idx].split_mps_l2r(s_avg, avg_truncation.keeps);
// #endif

// #ifdef _SAALG_VAR_USV_
//                         // sa_alg=-1 variant USV
//                         // Simple splitting of the average MPS tensor, to be applied to all states
//                         mps_vector[idx][site1] = split_site1;
//                         mps_vector[idx][site2] = split_site2;
// #endif

// #ifdef _SAALG_VAR_SV_
//                         // sa_alg=-1 variant SV
//                         // U is obtained from the state-specific TST but S*V from the state-average TST
//                         MPSTensor<Matrix, SymmGroup> jnk;
//                         std::tie(mps_vector[idx][site1], jnk, trunc[idx]) = two_vec[idx].split_mps_l2r(Mmax, cutoff, avg_truncation.keeps);
//                         mps_vector[idx][site2] = split_site2;
// #endif

//                         two_vec[idx].clear();
//                         block_matrix<Matrix, SymmGroup> t = mps_vector[idx][site2].normalize_left(DefaultSolver());
//                         if (site2 < L_-1)
//                             mps_vector[idx][site2+1].multiply_from_left(t);
//                    }
                else if (sa_alg_ == -2) {
                    for (size_t idx = 0; idx < n_root_; idx++) {
                        if (parms["twosite_truncation"] == "svd")
                            std::tie(mps_vector[idx][site1], mps_vector[idx][site2], trunc[idx])
                                    = two_vec[idx].split_mps_l2r(Mmax, cutoff);
                        else
                            std::tie(mps_vector[idx][site1], mps_vector[idx][site2], trunc[idx])
                                    = two_vec[idx].predict_split_l2r(Mmax, cutoff, alpha, (*(boundaries_database_.get_boundaries_left(idx, false)))[site1], mpo[site1]);
                        two_vec[idx].clear();
                        block_matrix<Matrix, SymmGroup> t = mps_vector[idx][site2].normalize_left(DefaultSolver());
                        if (site2 < L_-1)
                            mps_vector[idx][site2+1].multiply_from_left(t);
                    }
                // A single root taken as reference for the truncation
                } else if (sa_alg_ > -1) {
                    // TODO: Leon to ALB: Shouldn't this use state # sa_alg_ as the reference state for the truncation and NOT the ground state?
                    // TODO: yes it should, see ss_optimize.hpp
                    // Truncation of the reference state
                    if (parms["twosite_truncation"] == "svd")
                        std::tie(mps_vector[0][site1], mps_vector[0][site2], trunc[0]) = two_vec[0].split_mps_l2r(Mmax, cutoff);
                    else
                        std::tie(mps_vector[0][site1], mps_vector[0][site2], trunc[0])
                                = two_vec[0].predict_split_l2r(Mmax, cutoff, alpha, (*(boundaries_database_.get_boundaries_left(0, false)))[site1], mpo[site1]);
                    two_vec[0].clear();
                    //Get # of eigenvalues to keep per block from state 0 and use the same number to truncate the other states
                    std::vector<size_t>& keeps = trunc[0].keeps;
                    block_matrix<Matrix, SymmGroup> t = mps_vector[0][site2].normalize_left(DefaultSolver());
                    if (site2 < L_-1)
                        mps_vector[0][site2+1].multiply_from_left(t);
                    // Truncation of the other states
                    for (size_t idx = 1; idx < n_root_; idx++) {
                        if (parms["twosite_truncation"] == "svd")
                            std::tie(mps_vector[idx][site1], mps_vector[idx][site2], trunc[idx])
                                    = two_vec[idx].split_mps_l2r(Mmax, cutoff, keeps);
                        else
                            std::tie(mps_vector[idx][site1], mps_vector[idx][site2], trunc[idx])
                                    = two_vec[idx].predict_split_l2r(Mmax, cutoff, alpha, (*(boundaries_database_.get_boundaries_left(0, false)))[site1], mpo[site1], keeps);
                        two_vec[idx].clear();
                        block_matrix<Matrix, SymmGroup> t = mps_vector[idx][site2].normalize_left(DefaultSolver());
                        if (site2 < L_-1)
                            mps_vector[idx][site2+1].multiply_from_left(t);
                    }
                }
                this->boundary_left_step(mpo, site1); // creating left_[site2]
    	      }
            // +--------------+
            //  Backward sweep
            // +--------------+
       	    if (lr == -1){
                if (sa_alg_ == -3 || sa_alg_ == -1) {
                    std::vector< MPSTensor<Matrix,SymmGroup> > vec_mps ;
                    MPSTensor<Matrix, SymmGroup> tst_mps ;
                    if (sa_alg_ == -3)
                    {
                        std::tie(vec_mps, tst_mps, trunc[0]) = split_mps_r2l_vector(two_vec, Mmax, cutoff);
                        int Mval = trunc[0].bond_dimension ;
                        std::cout << "Bond dimension " << Mval << std::endl ;
                    }
                    else //sa_alg_ == -1
                        std::tie(vec_mps, tst_mps, trunc[0]) = predict_split_r2l_avg(two_vec, Mmax, cutoff,
                            alpha, (*(boundaries_database_.get_boundaries_right(0, false)))[site2+1], mpo[site2]);

                    // Final update of the MPS
                    for (std::size_t idx = 0; idx < n_root_; idx++) {
                        two_vec[idx].clear() ;
                        block_matrix<Matrix, SymmGroup> t ;
                        (*(boundaries_database_.get_mps(idx)))[site2] = tst_mps ;
                        (*(boundaries_database_.get_mps(idx)))[site1] = vec_mps[idx] ;
                        t = (*(boundaries_database_.get_mps(idx)))[site1].normalize_right(DefaultSolver()) ;
                        if (site1 > 0)
                            (*(boundaries_database_.get_mps(idx)))[site1-1].multiply_from_right(t) ;
                    }

                }
//                 else if (sa_alg_ == -1) {
//                     // old sa_alg_ == -1 code that never worked correctly. if the above works, it should be deleted

//                     // Construct the two-site tensor averaged over all states
//                     TwoSiteTensor<Matrix, SymmGroup> avg_tst = std::accumulate(std::next(two_vec.begin()), two_vec.end(), *two_vec.begin());
//                     avg_tst /= n_root_;

//                     truncation_results avg_truncation;

//                     // See the comment in the forward sweep for the explanation of the defines
// #ifdef _SAALG_VAR_S_
//                     // sa_alg=-1 variant S
//                     // Perform truncation of the average TwoSiteTensor and obtain the S (diagonal) matrix

//                     typename TwoSiteTensor<Matrix, SymmGroup>::block_diag_matrix s_avg;
//                     std::tie(s_avg, avg_truncation) = avg_tst.get_S(Mmax, cutoff);
// #endif

// #if defined(_SAALG_VAR_SV_) || defined (_SAALG_VAR_USV_)
//                     // sa_alg=-1 variants SV and USV
//                     // Simple splitting of the average MPS tensor, to be applied to all states
//                     // Get U,S,V from the average
//                     MPSTensor<Matrix, SymmGroup> split_site1, split_site2;
//                     std::tie(split_site1, split_site2, avg_truncation) = avg_tst.split_mps_r2l(Mmax, cutoff);

// #endif

//                     for (size_t idx = 0; idx < n_root_; idx++) {
// #ifdef _SAALG_VAR_S_
//                         // sa_alg=-1 variant S
//                         // for an external S split_mps_l2r is overloaded and doesn't return truncation_results
//                         std::tie(mps_vector[idx][site1], mps_vector[idx][site2]) = two_vec[idx].split_mps_r2l(s_avg, avg_truncation.keeps);
// #endif

// #ifdef _SAALG_VAR_USV_
//                         // sa_alg=-1 variant USV
//                         // Simple splitting of the average MPS tensor, to be applied to all states
//                         mps_vector[idx][site1] = split_site1;
//                         mps_vector[idx][site2] = split_site2;
// #endif

// #ifdef _SAALG_VAR_SV_
//                         //sa_alg=-1 variant SV
//                         // V is obtained from the state-specific TST but U*S from the state-average TST
//                         MPSTensor<Matrix, SymmGroup> jnk;
//                         std::tie(jnk, mps_vector[idx][site2], trunc[idx]) = two_vec[idx].split_mps_r2l(Mmax, cutoff, avg_truncation.keeps);
//                         mps_vector[idx][site1] = split_site1;
// #endif

//                         two_vec[idx].clear();
//                         block_matrix<Matrix, SymmGroup> t = mps_vector[idx][site1].normalize_right(DefaultSolver());
//                         if (site1 > 0)
//                             mps_vector[idx][site1-1].multiply_from_right(t);
//                    }
                else if (sa_alg_ == -2) {
                    for (size_t idx = 0; idx < n_root_; idx++) {
                        if (parms["twosite_truncation"] == "svd")
                            std::tie(mps_vector[idx][site1], mps_vector[idx][site2], trunc[idx])
                                    = two_vec[idx].split_mps_r2l(Mmax, cutoff);
                        else
                            std::tie(mps_vector[idx][site1], mps_vector[idx][site2], trunc[idx])
                                    = two_vec[idx].predict_split_r2l(Mmax, cutoff, alpha, (*(boundaries_database_.get_boundaries_right(idx, false)))[site2+1], mpo[site2]);
                        two_vec[idx].clear();
                        block_matrix<Matrix, SymmGroup> t;
                        t = mps_vector[idx][site1].normalize_right(DefaultSolver());
                        if (site1 > 0)
                            mps_vector[idx][site1-1].multiply_from_right(t);
                    }
                } else if (sa_alg_ > -1) {
                    // Truncation of the reference state
                    if (parms["twosite_truncation"] == "svd")
                        std::tie(mps_vector[0][site1], mps_vector[0][site2], trunc[0])
                                = two_vec[0].split_mps_r2l(Mmax, cutoff);
                    else
                        std::tie(mps_vector[0][site1], mps_vector[0][site2], trunc[0])
                                = two_vec[0].predict_split_r2l(Mmax, cutoff, alpha, (*(boundaries_database_.get_boundaries_right(0, false)))[site2+1], mpo[site2]);
                    two_vec[0].clear();
                    //Get # of eigenvalues to keep per block from state 0 and use the same number to truncate the other states
                    std::vector<size_t>& keeps = trunc[0].keeps;
                    block_matrix<Matrix, SymmGroup> t = mps_vector[0][site1].normalize_right(DefaultSolver());
                    if (site2 < L_-1)
                        mps_vector[0][site1-1].multiply_from_right(t);
                    // Truncation of the other states
                    for (size_t idx = 1; idx < n_root_; idx++) {
                        if (parms["twosite_truncation"] == "svd")
                            std::tie(mps_vector[idx][site1], mps_vector[idx][site2], trunc[idx])
                                    = two_vec[idx].split_mps_r2l(Mmax, cutoff, keeps);
                        else
                            std::tie(mps_vector[idx][site1], mps_vector[idx][site2], trunc[idx])
                                    = two_vec[idx].predict_split_r2l(Mmax, cutoff, alpha, (*(boundaries_database_.get_boundaries_right(idx, false)))[site2+1], mpo[site1], keeps);
                        two_vec[idx].clear();
                        block_matrix<Matrix, SymmGroup> t = mps_vector[idx][site1].normalize_right(DefaultSolver());
                        if (site1 > 0)
                            mps_vector[idx][site1-1].multiply_from_right(t);
                    }
                }
                //
                this->boundary_right_step(mpo, site2);
            }
            END_TIMING("MPS TRUNCATION")
            // +------------+
            //  FINALIZATION
            // +------------+
            for (size_t i = 0; i < n_root_; i++) {
                sorter_[i].first  = res[i].first ;
                sorter_[i].second = i ;
                if (update_each_omega || update_omega && _site == 0 )
                    omega_vec[i] = res[i].first - omega_shift_ ;
            }
            std::sort(sorter_.begin(), sorter_.end()) ;
            this->update_order(sorter_) ;
            BEGIN_TIMING("PARTIAL OVERLAP UPDATE")
            if (root_homing_type_ == 1 || root_homing_type_ == 3) {
                for (size_t k = 0; k < n_root_; k++) {
                    for (size_t k1 = 0; k1 < L_; k1++) {
                        poverlap_vec_[k].update(mps_vector[k], k1, 1);
                        poverlap_vec_[k].update(mps_vector[k], L_-k1-1, -1);
                    }
                }
            }
// Leon: root_homing_type != 0 not supported yet, therefore below is commented out
/*            if (root_homing_type_ == 1 || root_homing_type_ == 3)
                for (size_t k = 0 ; k < n_root_ ; k++)
                    poverlap_vec_[k].update(mps_vector[k], site, lr) ;*/
            for (size_t k = 0 ; k < n_root_ ; k++)
            {
                iteration_results_[k]["BondDimension"]     << trunc[k].bond_dimension;
                iteration_results_[k]["TruncatedWeight"]   << trunc[k].truncated_weight;
                iteration_results_[k]["TruncatedFraction"] << trunc[k].truncated_fraction;
                iteration_results_[k]["SmallestEV"]        << trunc[k].smallest_ev;
            }
            parallel::meminfo();
            boost::chrono::high_resolution_clock::time_point sweep_then = boost::chrono::high_resolution_clock::now();
            double elapsed = boost::chrono::duration<double>(sweep_then - sweep_now).count();
            maquis::cout << "Sweep has been running for " << elapsed << " seconds." << std::endl;
            if (stop_callback())
                throw dmrg::time_limit(sweep, _site+1);
        }
        initial_site = -1;
    } // sweep

private:
    // --------------------
    //  Private attributes
    // --------------------
    int initial_site;
    MPO<Matrix, SymmGroup> ts_cache_mpo, ts_cache_mpo_squared;
    // -------------------------------------
    //  Simple function to print the header
    // -------------------------------------
    void print_header(int& sweep, int& site1, int& site2, int& lr){
        char buffer[50] ;
        int n , a ;
        if (lr == 1) {
            a = 2*sweep+1 ;
            n = sprintf(buffer, "  Sweep number %3d - site numbers %3d and %3d", a, site1, site2);
        } else {
            a = 2*sweep+2 ;
            n = sprintf(buffer, "  Sweep number %3d - site numbers %3d and %3d", a, site1, site2);
        }
        std::cout << " +--------------------------------------------+" << std::endl ;
        std::cout << buffer << std::endl ;
        std::cout << " +--------------------------------------------+" << std::endl ;
    }

};

#endif
