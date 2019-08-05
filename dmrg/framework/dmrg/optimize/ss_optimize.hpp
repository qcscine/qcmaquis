/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *               2017-2018 by Alberto Baiardi <alberto.baiardi@sns.it>
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

#ifndef SS_OPTIMIZE_H
#define SS_OPTIMIZE_H

#include "dmrg/mp_tensors/mpo_ops.h"
#include "dmrg/optimize/optimize.h"
#include "dmrg/optimize/vectorset.h"
#include "dmrg/utils/DmrgParameters.h"

// +-------------------+
//   SS_OPTIMIZE CLASS
// +-------------------+
// Inherited from the virtual optimizer_base class
template<class Matrix, class SymmGroup, class Storage>
class ss_optimize : public optimizer_base<Matrix, SymmGroup, Storage>
{
public:
    typedef typename Matrix::value_type value_type;
    // Inherits several data from the "mother" class
    typedef optimizer_base<Matrix, SymmGroup, Storage> base;
    //
    using base::boundaries_database_ ;
    using base::chebyshev_shift_ ;
    using base::correction_equation ;
    using base::do_chebyshev_ ;
    using base::do_H_squared_ ;
    using base::do_root_homing_ ;
    using base::do_shiftandinvert_ ;
    using base::do_stateaverage_ ;
    using base::finalizer_ ;
    using base::is_folded_ ;
    using base::iteration_results_ ;
    using base::left_sa_ ;
    using base::left_squared_sa_ ;
    using base::L_ ;
    using base::micro_optimizer ;
    using base::mpo ;
    using base::mpo_squared ;
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

    // Constructor declaration
    ss_optimize(std::vector< MPS<Matrix, SymmGroup> > & mps_vector_ ,
                MPO<Matrix, SymmGroup> const & mpo_,
                BaseParameters & parms_,
                boost::function<bool ()> stop_callback_,
                int initial_site_,
                std::vector< MPSTensor<Matrix, SymmGroup> > & mps_guess_ ,
                std::vector< MPS<Matrix, SymmGroup> > & mps_partial_overlap_ ,
                MPO<Matrix, SymmGroup> const & mpo_squared_)
    : base(mps_vector_, mpo_, parms_, stop_callback_,
          to_site(mps_vector_[0].length(), initial_site_), true, boost::none,
           mps_guess_, mps_partial_overlap_, mpo_squared_)
    , initial_site((initial_site_ < 0) ? 0 : initial_site_)
    {
        iteration_results_.resize(n_root_);
    };
    // Inline function to get the site index modulo 2
    inline int to_site(const int L, const int i) const
    {
        if (i < 0) return 0;
        return (i < L) ? i : 2*L - 1 - i;
    }
    //
    // SWEEP ROUTINE
    // -------------
    void sweep(int sweep, OptimizeDirection d = Both) {
        // Initialization
        boost::chrono::high_resolution_clock::time_point sweep_now = boost::chrono::high_resolution_clock::now();
        for (auto&& it : iteration_results_)
            it.clear();
        std::size_t L = mps_vector[0].length();
        for (size_t i = 1; i < n_root_; i++)
            assert(mps_vector[i].length() == L);
        //
        int _site = 0, site = 0;
        if (initial_site != -1) {
            _site = initial_site;
            site = to_site(L, _site);
        }
        // Main loop
        this->update_parameters(sweep) ;
        for (; _site < 2 * L; ++_site) {
            boost::chrono::high_resolution_clock::time_point now, then;
            BEGIN_TIMING("PRELIMINARY OPERATIONS")
            // lr indicates the direction of the sweep
            int lr = (_site < L) ? +1 : -1;
            site = to_site(L, _site);
            print_header(sweep, site, lr);
            if (lr == -1 && site == L - 1) {
                maquis::cout << "Syncing storage" << std::endl;
                Storage::sync();
            }
            // Storage
            Storage::fetch(left_sa_[0][site]) ;
            Storage::fetch(right_sa_[0][site+1]) ;
            if (lr == +1 && site+2 <= L)
                Storage::prefetch(right_sa_[0][site+2]);
            if (lr == -1 && site > 0)
                Storage::prefetch(left_sa_[0][site-1]);
            assert( left_sa_[0][site].reasonable() );    // in case something went wrong
            assert( right_sa_[0][site+1].reasonable() ); // in case something went wrong
            // Debug printing
            if (do_shiftandinvert_) {
                std::cout << " -- VALUES OF OMEGA -- " << std::endl ;
                for (size_t idx = 0; idx < n_root_; idx++)
                    std::cout << omega_vec[idx] << std::endl ;
            }
            std::vector< std::pair<typename maquis::traits::real_type<value_type>::type, MPSTensor<Matrix, SymmGroup> > > res;
            std::vector< MPSTensor<Matrix, SymmGroup> > site_vec(n_root_) ;
            for (std::size_t idx = 0; idx < n_root_; idx++)
                site_vec[idx] = mps_vector[idx][site] ;
            SiteProblem<Matrix, SymmGroup> sp(site_vec, mpo[site], mpo_squared()[site], site, site+1, boundaries_database_,
                                              do_H_squared_);
            if (sa_alg_ != -3 || sa_alg_ != -1) {
                mps_guess().clear() ;
                for (std::size_t idx = 0; idx < mps_vector.size(); idx++)
                    mps_guess().push_back(mps_vector[idx][site]) ;
            }
            // Generates the vectorset object
            VectorSet<Matrix, SymmGroup> vector_set(mps_guess());
            // Compute orthogonal vectors
            std::vector<MPSTensor<Matrix, SymmGroup> > ortho_vecs(base::northo);
            for (int n = 0; n < base::northo; ++n) {
                ortho_vecs[n] = contraction::site_ortho_boundaries(mps_vector[0][site],
                                                                   base::ortho_mps[n][site],
                                                                   base::ortho_left_[n][site],
                                                                   base::ortho_right_[n][site + 1]);
            }
            END_TIMING("PRELIMINARY OPERATIONS")
            // +-----------------------------+
            //  MAIN PART: performs the sweep
            // +-----------------------------+
            if (d == Both || (d == LeftOnly && lr == -1) || (d == RightOnly && lr == +1)) {
                BEGIN_TIMING("LOCAL DIAGONALIZATION")
                res = solve_ietl_jcd(sp, vector_set, correction_equation, micro_optimizer, *finalizer_,
                                     orthogonalizer_, parms, poverlap_vec_, site, site, root_homing_type_, do_root_homing_,
                                     do_shiftandinvert_, do_chebyshev_, chebyshev_shift_, do_H_squared_, reshuffle_variance_,
                                     track_variance_, is_folded_, vec_sa_left_, vec_sa_right_, order, boundaries_database_,
									 sa_alg_, omega_vec, ortho_vecs) ;
                END_TIMING("LOCAL DIAGONALIZATION")
                BEGIN_TIMING("MPS UPDATE")
                // Correct the energies
                for (size_t k = 0; k < n_root_; k++) {
                    res[k].first = ietl::get_energy(sp, res[k].second, k, false);
                    mps_vector[k][site] = res[k].second;
                    // TODO: usage of mps_guess() is very confusing and most likely also redundant -- get rid of it
                    mps_guess()[k] = res[k].second ;
                }
                END_TIMING("MPS UPDATE")
            }
            // +---------------------+
            //  Collection of results
            // +---------------------+
            int prec = maquis::cout.precision();
            for (size_t k = 0; k < n_root_; k++)
            {
                maquis::cout.precision(15);
                maquis::cout << " Converged energy - state " << k << " = " << res[k].first + mpo.getCoreEnergy()
                        << std::endl;
                maquis::cout.precision(prec);
                iteration_results_[k]["Energy"] << res[k].first + mpo.getCoreEnergy();

            }
            // Loads the alpha parameter
            double alpha;
            int ngs = parms.template get<int>("ngrowsweeps"), nms = parms.template get<int>("nmainsweeps");
            if (sweep < ngs)
                alpha = parms.template get<double>("alpha_initial");
            else if (sweep < ngs + nms)
                alpha = parms.template get<double>("alpha_main");
            else
                alpha = parms.template get<double>("alpha_final");
            // Update of the MPS
            double cutoff = this->get_cutoff(sweep);
            std::size_t Mmax = this->get_Mmax(sweep), Mval;
            std::vector<truncation_results> trunc(n_bound_);
            // +---------------------+
            //  TRUNCATION OF THE MPS
            // +---------------------+
            // -- FORWARD SWEEP --
            BEGIN_TIMING("MPS TRUNCATION")
            if (lr == +1) {
                if (site < L - 1) {
                    maquis::cout << " Alpha = " << alpha << std::endl;
                    if (sa_alg_ == -3) {
                        trunc[0] = (*(boundaries_database_.get_mps(0))).grow_l2r_sweep_vec(mps_guess(), mpo[site],
                                                  (*(boundaries_database_.get_boundaries_left(0, false)))[site],
                                                  (*(boundaries_database_.get_boundaries_right(0, false)))[site+1],
                                                  site, alpha, cutoff, Mmax);
                        for (std::size_t idx = 1; idx < n_root_; idx++) {
                            // after calling grow_l2r_sweep_, mps_guess[idx] contains the correct mpstensor at [site+1]
                            (*(boundaries_database_.get_mps(idx)))[site + 1] = mps_guess()[idx];
                            (*(boundaries_database_.get_mps(idx)))[site] = (*(boundaries_database_.get_mps(0)))[site] ;
                        }

                    } else if (sa_alg_ == -1) {
                        // similar to sa_alg=-3 but use grow_l2r_sweep_avg instead of grow_l2r_sweep_vec
                        // TODO: reduce the level of copy-paste
                        trunc[0] = (*(boundaries_database_.get_mps(0))).grow_l2r_sweep_avg(mps_guess(), mpo[site],
                                                  (*(boundaries_database_.get_boundaries_left(0, false)))[site],
                                                  (*(boundaries_database_.get_boundaries_right(0, false)))[site+1],
                                                  site, alpha, cutoff, Mmax);
                        for (std::size_t idx = 1; idx < n_root_; idx++) {
                            (*(boundaries_database_.get_mps(idx)))[site + 1] = mps_guess()[idx];
                            (*(boundaries_database_.get_mps(idx)))[site] = (*(boundaries_database_.get_mps(0)))[site] ;
                        }
                    } else if (sa_alg_ == -2) {
                        for (size_t idx = 0; idx < n_bound_; idx++)
                            trunc[idx] = (*(boundaries_database_.get_mps(idx))).grow_l2r_sweep(mpo[site],
                                                                                               (*(boundaries_database_.get_boundaries_left(idx, false)))[site],
                                                                                               (*(boundaries_database_.get_boundaries_right(idx, false)))[site+1],
                                                                                               site, alpha, cutoff, Mmax);
                    // -- SINGLE STATE TAKEN AS REFERENCE FOR BOUNDARY CALCULATION --
                    } else if (sa_alg_ > -1) {
                        for (size_t idx = 0; idx < n_bound_; idx++)
                            trunc[idx] = (*(boundaries_database_.get_mps(idx))).grow_l2r_sweep(mpo[site],
                                                                                               (*(boundaries_database_.get_boundaries_left(idx, false)))[site],
                                                                                               (*(boundaries_database_.get_boundaries_right(idx, false)))[site+1],
                                                                                               site, alpha, cutoff, Mmax);
                        std::vector<size_t>& keeps = trunc[0].keeps;
                        for (size_t k = 0; k < n_root_; k++) {
                            if (k != sa_alg_)
                                trunc.push_back(mps_vector[k].grow_l2r_sweep(mpo[site],
                                                                             (*(boundaries_database_.get_boundaries_left(k, false)))[site],
                                                                             (*(boundaries_database_.get_boundaries_right(k,false)))[site+1],
                                                                             site, alpha, cutoff, Mmax, keeps));
                        }
                    }
                } else {
                    if (sa_alg_ == -3 || sa_alg_ == -1) {
                        for (std::size_t idx = 0; idx < n_root_; idx++) {
                            (*(boundaries_database_.get_mps(idx)))[site] = mps_guess()[idx];
                                mps_vector[idx][site].normalize_left(DefaultSolver());
                        }
                    } else {
                        for (std::size_t idx = 0; idx < n_root_; idx++)
                            mps_vector[idx][site].normalize_left(DefaultSolver());
                    }
                }
                // Update the left boundary
                this->boundary_left_step(mpo, site);
            // -- BACKWARD SWEEP --
            } else if (lr == -1) {
                if (site > 0) {
                    maquis::cout << " Alpha = " << alpha << std::endl;
                    if (sa_alg_ == -3) {
                        trunc[0] = (*(boundaries_database_.get_mps(0))).grow_r2l_sweep_vec(mps_guess(), mpo[site],
                                                                                           (*(boundaries_database_.get_boundaries_left(0, false)))[site],
                                                                                           (*(boundaries_database_.get_boundaries_right(0, false)))[site+1],
                                                                                           site, alpha, cutoff, Mmax);
                        for (std::size_t idx = 1; idx < n_root_; idx++) {
                            // after calling grow_l2r_sweep_, mps_guess[idx] contains the correct mpstensor at [site+1], not at [site]
                            (*(boundaries_database_.get_mps(idx)))[site-1] = mps_guess()[idx];
                            (*(boundaries_database_.get_mps(idx)))[site] = (*(boundaries_database_.get_mps(0)))[site] ;
                        }
                    } else if (sa_alg_ == -1) {
                        trunc[0] = (*(boundaries_database_.get_mps(0))).grow_r2l_sweep_avg(mps_guess(), mpo[site],
                                                                                           (*(boundaries_database_.get_boundaries_left(0, false)))[site],
                                                                                           (*(boundaries_database_.get_boundaries_right(0, false)))[site+1],
                                                                                           site, alpha, cutoff, Mmax);
                        for (std::size_t idx = 1; idx < n_root_; idx++) {
                            (*(boundaries_database_.get_mps(idx)))[site-1] = mps_guess()[idx];
                            (*(boundaries_database_.get_mps(idx)))[site] = (*(boundaries_database_.get_mps(0)))[site] ;
                        }
                    } else if (sa_alg_ == -2) {
                        for (size_t idx = 0; idx < n_bound_; idx++)
                            trunc[idx] = (*(boundaries_database_.get_mps(idx))).grow_r2l_sweep(mpo[site],
                                                                                               (*(boundaries_database_.get_boundaries_left(idx, false)))[site],
                                                                                               (*(boundaries_database_.get_boundaries_right(idx, false)))[site+1],
                                                                                               site, alpha, cutoff, Mmax);
                    } else if (sa_alg_ > -1) {
                        for (size_t idx = 0; idx < n_bound_; idx++)
                            trunc[idx] = (*(boundaries_database_.get_mps(idx))).grow_r2l_sweep(mpo[site],
                                                                                               (*(boundaries_database_.get_boundaries_left(idx, false)))[site],
                                                                                               (*(boundaries_database_.get_boundaries_right(idx, false)))[site+1],
                                                                                               site, alpha, cutoff, Mmax);
                        std::vector<size_t>& keeps = trunc[0].keeps;
                        for (size_t k = 0; k < n_root_; k++)
                            if (k != sa_alg_) {
                                trunc.push_back(mps_vector[k].grow_r2l_sweep(mpo[site],
                                                                             (*(boundaries_database_.get_boundaries_left(k, false)))[site],
                                                                             (*(boundaries_database_.get_boundaries_right(k, false)))[site+1],
                                                                             site, alpha, cutoff, Mmax, keeps));
                            }
                    }
                } else {
                    // -- TERMINAL SITE OF THE DMRG CHAIN --
                    if (sa_alg_ == -3 || sa_alg_ == -1) {
                        // here mps_guess contains MPSTensor at site [site]
                        for (std::size_t idx = 0; idx < n_root_; idx++) {
                            (*(boundaries_database_.get_mps(idx)))[site] = mps_guess()[idx];
                                mps_vector[idx][site].normalize_right(DefaultSolver());
                        }
                    } else {
                        for (std::size_t idx = 0; idx < n_root_; idx++)
                            mps_vector[idx][site].normalize_right(DefaultSolver());
                    }
                }
                // Update the right boundary
                this->boundary_right_step(mpo, site);
            }
            END_TIMING("MPS TRUNCATION")
            // +----------------+
            //  Final operations
            // +----------------+
            BEGIN_TIMING("FINAL OPERATIONS")
            for (size_t i = 0; i < n_root_; i++) {
                sorter_[i].first = res[i].first;
                sorter_[i].second = i;
                if ( (update_each_omega || update_omega && _site == 0) )
                    omega_vec[i] = res[i].first - omega_shift_ ;
            }
            std::sort(sorter_.begin(), sorter_.end());
            this->update_order(sorter_);
            //
            if (root_homing_type_ == 1 || root_homing_type_ == 3) {
                for (size_t k = 0; k < n_root_; k++) {
                    for (size_t k1 = 0; k1 < L_; k1++) {
                        poverlap_vec_[k].update(mps_vector[k], k1, 1);
                        poverlap_vec_[k].update(mps_vector[k], L_-k1-1, -1);
                    }
                }
            }
            END_TIMING("FINAL OPERATIONS")
            //
            for (size_t k = 0; k < n_root_; k++)
            {
                iteration_results_[k]["BondDimension"] << trunc[k].bond_dimension;
                iteration_results_[k]["TruncatedWeight"] << trunc[k].truncated_weight;
                iteration_results_[k]["SmallestEV"] << trunc[k].smallest_ev;
            }
            boost::chrono::high_resolution_clock::time_point sweep_then = boost::chrono::high_resolution_clock::now();
            double elapsed = boost::chrono::duration<double>(sweep_then - sweep_now).count();
            maquis::cout << "Sweep has been running for " << elapsed << " seconds." << std::endl;
            if (stop_callback())
                throw dmrg::time_limit(sweep, _site + 1);
        }
        initial_site = -1;
    }
    //
private:
    // --------------------
    //  Private attributes
    // --------------------
    int initial_site ;
    //
    // Simple function to print the header
    // -----------------------------------
    void print_header(int& sweep, int& site, int& lr){
        char buffer[40] ;
        int n , a ;
        if (lr == 1) {
            a = 2*sweep+1 ;
            n = sprintf(buffer, "  Sweep number %3d - site number %3d", a, site);
        } else {
            a = 2*sweep+2 ;
            n = sprintf(buffer, "  Sweep number %3d - site number %3d", a, site);
        }
        std::cout << " +-----------------------------------+" << std::endl ;
        std::cout << buffer << std::endl ;
        std::cout << " +-----------------------------------+" << std::endl ;
    }
};

#endif

