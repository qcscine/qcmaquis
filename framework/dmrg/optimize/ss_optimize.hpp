/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *               2017 by Alberto Baiardi <alberto.baiardi@sns.it>
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

//
// SS_OPTIMIZE CLASS
// -----------------
// Inherited from the virtual optimizer_base class
template<class Matrix, class SymmGroup, class Storage>
class ss_optimize : public optimizer_base<Matrix, SymmGroup, Storage>
{
public:
    // Inherits several data from the "mother" class
    typedef optimizer_base<Matrix, SymmGroup, Storage> base;
    //
    using base::do_root_homing_ ;
    using base::do_stateaverage_ ;
    using base::mpo ;
    using base::mps ;
    using base::mps2follow ;
    using base::mps_vector ;
    using base::n_root_ ;
    using base::order ;
    using base::left_ ;
    using base::left_sa_ ;
    using base::right_ ;
    using base::right_sa_ ;
    using base::parms ;
    using base::poverlap_vec_ ;
    using base::iteration_results_ ;
    using base::stop_callback ;
    using base::root_homing_type_ ;
    using base::update_order ;
    using base::vec_sa_left_ ;
    using base::vec_sa_right_ ;
    // Constructor declaration
    ss_optimize(MPS<Matrix, SymmGroup> & mps_, 
                std::vector< MPS<Matrix, SymmGroup> > & mps_vector_ ,
                MPO<Matrix, SymmGroup> const & mpo_,
                BaseParameters & parms_,
                boost::function<bool ()> stop_callback_,
                int initial_site_ = 0)
    : base(mps_, mps_vector_, mpo_, parms_, stop_callback_, to_site(mps_vector_[0].length(), initial_site_))
    , initial_site((initial_site_ < 0) ? 0 : initial_site_)
    { };
    // Inline function to get the site index modulo 2
    inline int to_site(const int L, const int i) const
    {
        if (i < 0) return 0;
        return (i < L) ? i : 2*L - 1 - i;
    }
    //
    // SWEEP ROUTINE
    // -------------
    void sweep(int sweep, OptimizeDirection d = Both)
    {
        // Initialization
        boost::chrono::high_resolution_clock::time_point sweep_now = boost::chrono::high_resolution_clock::now();
        iteration_results_.clear();
        std::size_t L = mps_vector[0].length();
        for (size_t i = 1 ; i < n_root_ ; i++)
            assert(mps_vector[i].length() == L);
        //
        int _site = 0, site = 0;
        if (initial_site != -1) {
            _site = initial_site;
            site = to_site(L, _site);
        }
        //TODO ALB This should be done with pointers, very stupid
        std::vector< std::pair<float,int> > sorter ;
        sorter.resize(n_root_) ;
        // Initialization of the overlap object
        Storage::prefetch(left_[site]) ;
        Storage::prefetch(right_[site+1]) ;
        // Main loop
        for (; _site < 2*L; ++_site) {
            boost::chrono::high_resolution_clock::time_point now, then;
            BEGIN_TIMING("PRELIMINARY OPERATIONS")
            // lr indicates the direction of the sweep
            int lr = (_site < L) ? +1 : -1;
            site = to_site(L, _site);
            print_header(sweep, site, lr) ;
            if (lr == -1 && site == L-1) {
                maquis::cout << "Syncing storage" << std::endl;
                Storage::sync();
            }
            Storage::fetch(left_[site]);
            Storage::fetch(right_[site+1]);
            if (lr == +1 && site+2 <= L) Storage::prefetch(right_[site+2]);
            if (lr == -1 && site > 0)    Storage::prefetch(left_[site-1]);
            assert( left_[site].reasonable() );    // in case something went wrong
            assert( right_[site+1].reasonable() ); // in case something went wrong
            std::vector<std::pair<double, MPSTensor<Matrix, SymmGroup> > > res;
            SiteProblem<Matrix, SymmGroup> sp(left_sa_, right_sa_, mpo[site], site, site+1);
            // Generates the vectorset object
            VectorSet<Matrix,SymmGroup> vector_set(mps_vector, site) ;
            // Compute orthogonal vectors
            std::vector<MPSTensor<Matrix, SymmGroup> > ortho_vecs(base::northo);
            for (int n = 0; n < base::northo; ++n) {
                ortho_vecs[n] = contraction::site_ortho_boundaries(mps_vector[0][site],
                                                                   base::ortho_mps[n][site],
                                                                   base::ortho_left_[n][site],
                                                                   base::ortho_right_[n][site+1]);
            }
            END_TIMING("PRELIMINARY OPERATIONS")
            // +-----------------------------+
            //  MAIN PART: performs the sweep
            // +-----------------------------+
            if (d == Both ||
                (d == LeftOnly && lr == -1) ||
                (d == RightOnly && lr == +1))
            {
                if (parms["eigensolver"] == std::string("IETL")) {
                    BEGIN_TIMING("IETL")
                    //res = solve_ietl_lanczos(sp, mps[site], parms);
                    END_TIMING("IETL")
                } else if (parms["eigensolver"] == std::string("IETL_JCD")) {
                    BEGIN_TIMING("JCD")
                    res = solve_ietl_jcd(sp, vector_set, parms, poverlap_vec_, 1, site, site, root_homing_type_,
                                         vec_sa_left_, vec_sa_right_, order, ortho_vecs);
                    END_TIMING("JCD")
                } else if (parms["eigensolver"] == std::string("IETL_DAVIDSON")) {
                    BEGIN_TIMING("DAVIDSON")
                    res = solve_ietl_davidson(sp, vector_set, parms, poverlap_vec_, 1, site, site, root_homing_type_,
                                              vec_sa_left_, vec_sa_right_, order, ortho_vecs);
                    END_TIMING("DAVIDSON")
                } else {
                    throw std::runtime_error("I don't know this eigensolver.");
                }
                BEGIN_TIMING("MPS UPDATE")
                // Collects the results
                if (n_root_ > 0) {
                    mps[site]  = res[0].second ;
                    mps_vector[0][site]  = res[0].second ;
                    for (size_t k = 1; k < n_root_; k++) {
                        mps_vector[k][site] = res[k].second;
                        //mps[site] += res[k].second;
                    }
                    //mps[site] /= n_root_ ;
                }
                END_TIMING("MPS UPDATE")
            }
            // +---------------------+
            //  Collection of results
            // +---------------------+
            int prec = maquis::cout.precision();
            maquis::cout.precision(15);
            for (size_t k = 0 ; k < n_root_ ; k++)
                maquis::cout << " Converged energy - state " << k << " = " << res[k].first + mpo.getCoreEnergy() << std::endl;
            maquis::cout.precision(prec);
            iteration_results_["Energy"] << res[0].first + mpo.getCoreEnergy();
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
            truncation_results trunc;
            std::vector<truncation_results> trunc_sa ;
            // +---------------------+
            //  Truncation of the MPS
            // +---------------------+
            BEGIN_TIMING("MPS TRUNCATION")
            if (lr == +1) {
                if (site < L-1) {
                    maquis::cout << " Alpha = " << alpha << std::endl;
                    trunc = mps.grow_l2r_sweep(mpo[site], left_[site], right_[site+1], site, alpha, cutoff, Mmax) ;
                    Mval = trunc.bond_dimension ;
                    for (size_t k = 0 ; k < n_root_ ; k++) {
                        //trunc_sa.push_back(mps_vector[k].grow_l2r_sweep(mpo[site], left_sa_[k][site], right_sa_[k][site+1], site, alpha, cutoff, Mmax, Mval));
                        trunc_sa.push_back(mps_vector[k].grow_l2r_sweep(mpo[site], left_sa_[k][site], right_sa_[k][site+1], site, alpha, cutoff, Mmax));
                    }
                } else {
                    block_matrix<Matrix, SymmGroup> t = mps[site].normalize_left(DefaultSolver());
                    if (site < L-1)
                        mps[site+1].multiply_from_left(t);
                    for (size_t k = 0 ; k < n_root_ ; k++) {
                        t = mps_vector[k][site].normalize_left(DefaultSolver());
                        if (site < L-1)
                            mps_vector[k][site+1].multiply_from_left(t);
                    }
                }
                // Update the left boundary
                this->boundary_left_step(mpo, site);
                if (site != L-1) {
                    Storage::drop(right_[site+1]);
                    Storage::evict(left_[site]);
                }
            } else if (lr == -1) {
                if (site > 0) {
                    maquis::cout << " Alpha = " << alpha << std::endl;
                    trunc = mps.grow_r2l_sweep(mpo[site], left_[site], right_[site+1], site, alpha, cutoff, Mmax);
                    Mval = trunc.bond_dimension ;
                    for (size_t k = 0 ; k < n_root_ ; k++)
                        //trunc_sa.push_back(mps_vector[k].grow_r2l_sweep(mpo[site], left_sa_[k][site], right_sa_[k][site+1], site, alpha, cutoff, Mmax, Mval)) ;
                        trunc_sa.push_back(mps_vector[k].grow_r2l_sweep(mpo[site], left_sa_[k][site], right_sa_[k][site+1], site, alpha, cutoff, Mmax));
                } else {
                    block_matrix<Matrix, SymmGroup> t = mps[site].normalize_right(DefaultSolver());
                    if (site > 0)
                        mps[site-1].multiply_from_right(t);
                    for (size_t k = 0 ; k < n_root_ ; k++){
                        t = mps_vector[k][site].normalize_right(DefaultSolver());
                        if (site > 0)
                            mps_vector[k][site-1].multiply_from_right(t);
                    }
                }
                // Update the right boundary
                this->boundary_right_step(mpo, site);
                if (site > 0) {
                    Storage::drop(left_[site]);
                    Storage::evict(right_[site+1]);
                }
            }
            END_TIMING("MPS TRUNCATION")
            // +----------------+
            //  Final operations
            // +----------------+
            BEGIN_TIMING("FINAL OPERATIONS")
            for (size_t i = 0; i < n_root_; i++) {
                sorter[i].first  = res[i].first ;
                sorter[i].second = i ;
            }
            std::sort(sorter.begin(), sorter.end()) ;
            this->update_order(sorter) ;
            //
            if (root_homing_type_ == 1)
                for (size_t k = 0 ; k < n_root_ ; k++)
                    poverlap_vec_[k].update(mps_vector[k], site, lr) ;
            END_TIMING("FINAL OPERATIONS")
            //
            iteration_results_["BondDimension"]   << trunc.bond_dimension;
            iteration_results_["TruncatedWeight"] << trunc.truncated_weight;
            iteration_results_["SmallestEV"]      << trunc.smallest_ev;
            boost::chrono::high_resolution_clock::time_point sweep_then = boost::chrono::high_resolution_clock::now();
            double elapsed = boost::chrono::duration<double>(sweep_then - sweep_now).count();
            maquis::cout << " Sweep has been running for " << elapsed << " seconds. \n" << std::endl;
            if (stop_callback())
                throw dmrg::time_limit(sweep, _site+1);
        }
        initial_site = -1;
    }
    //
    // Simple function to print the header
    // -----------------------------------
    void print_header(int& sweep, int& site, int& lr){
        char buffer[40] ;
	    int n , a ;
        if (lr == 1) {
	        a = 2*sweep+1 ;
            n = sprintf(buffer, "  Sweeep number %3d - site number %3d", a, site);
        } else {
	        a = 2*sweep+2 ;
            n = sprintf(buffer, "  Sweeep number %3d - site number %3d", a, site);
	    }
        std::cout << " +-----------------------------------+" << std::endl ;
        std::cout << buffer << std::endl ;
        std::cout << " +-----------------------------------+" << std::endl ;
    }
// +----------+
//  Attributes
// +----------+
private:
    int initial_site;
};

#endif

