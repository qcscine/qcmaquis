/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
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
#include "dmrg/optimize/partial_overlap.h"
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
    using base::mpo ;
    using base::mps ;
    using base::left_ ;
    using base::right_ ;
    using base::parms ;
    using base::iteration_results_ ;
    using base::stop_callback ;
    using base::mps2follow ;
    // Constructor declaration
    ss_optimize(MPS<Matrix, SymmGroup> & mps_,
                MPO<Matrix, SymmGroup> const & mpo_,
                BaseParameters & parms_,
                boost::function<bool ()> stop_callback_,
                int initial_site_ = 0)
    : base(mps_, mpo_, parms_, stop_callback_, to_site(mps_.length(), initial_site_))
    , initial_site((initial_site_ < 0) ? 0 : initial_site_)
    { };
    // Inline function to get the site index modulo 2
    inline int to_site(const int L, const int i) const
    {
        if (i < 0) return 0;
        return (i < L) ? i : 2*L - 1 - i;
    }
    typedef typename partial_overlap<Matrix,SymmGroup>::partial_overlap partial_overlap ;
    //
    // SWEEP ROUTINE
    // -------------
    void sweep(int sweep, OptimizeDirection d = Both)
    {
        // Some initialization
        boost::chrono::high_resolution_clock::time_point sweep_now = boost::chrono::high_resolution_clock::now();
        iteration_results_.clear();
        std::size_t L = mps.length();
        int _site = 0, site = 0;
        if (initial_site != -1) {
            _site = initial_site;
            site = to_site(L, _site);
        }
        // Initialization of the overlap object
        partial_overlap poverlap(mps,mps2follow) ;
        Storage::prefetch(left_[site]) ;
        Storage::prefetch(right_[site+1]) ;
        // Main loop
        for (; _site < 2*L; ++_site) {
            //
            double i ;
            if (poverlap.is_defined()) {
                i = poverlap.overlap(site);
                std::cout << "Overlap " << i << std::endl;
            }
            // lr indicates the direction of the sweep
            int lr = (_site < L) ? +1 : -1;
            site = to_site(L, _site);
            if (lr == -1 && site == L-1) {
                maquis::cout << "Syncing storage" << std::endl;
                Storage::sync();
            }
            maquis::cout << "Sweep " << sweep << ", optimizing site " << site << std::endl;
            Storage::fetch(left_[site]);
            Storage::fetch(right_[site+1]);
            if (lr == +1 && site+2 <= L) Storage::prefetch(right_[site+2]);
            if (lr == -1 && site > 0)    Storage::prefetch(left_[site-1]);
            assert( left_[site].reasonable() );    // in case something went wrong
            assert( right_[site+1].reasonable() ); // in case something went wrong
            boost::chrono::high_resolution_clock::time_point now, then;
            std::pair<double, MPSTensor<Matrix, SymmGroup> > res;
            SiteProblem<Matrix, SymmGroup> sp(left_[site], right_[site+1], mpo[site]);
            // Compute orthogonal vectors
            std::vector<MPSTensor<Matrix, SymmGroup> > ortho_vecs(base::northo);
            for (int n = 0; n < base::northo; ++n) {
                ortho_vecs[n] = contraction::site_ortho_boundaries(mps[site],
                                                                   base::ortho_mps[n][site],
                                                                   base::ortho_left_[n][site],
                                                                   base::ortho_right_[n][site+1]);
            }
            //
            // MAIN PART: performs the sweep
            // -----------------------------
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
                    res = solve_ietl_jcd(sp, mps[site], parms, site, poverlap, ortho_vecs);
                    END_TIMING("JCD")
                } else if (parms["eigensolver"] == std::string("IETL_DAVIDSON")) {
                    BEGIN_TIMING("DAVIDSON")
                    res = solve_ietl_davidson(sp, mps[site], parms, site, poverlap, ortho_vecs);
                    END_TIMING("DAVIDSON")
                } else {
                    throw std::runtime_error("I don't know this eigensolver.");
                }
                mps[site] = res.second;
            }
            //
            // Collection of results
            // ---------------------
            int prec = maquis::cout.precision();
            maquis::cout.precision(15);
            maquis::cout << "Energy " << lr << " " << res.first + mpo.getCoreEnergy()<< std::endl;
            maquis::cout.precision(prec);
            iteration_results_["Energy"] << res.first + mpo.getCoreEnergy();
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
            std::size_t Mmax = this->get_Mmax(sweep);
            truncation_results trunc;
            //
            if (lr == +1) {
                if (site < L-1) {
                    maquis::cout << "Growing, alpha = " << alpha << std::endl;
                    trunc = mps.grow_l2r_sweep(mpo[site], left_[site], right_[site+1], site, alpha, cutoff, Mmax);
                } else {
                    block_matrix<Matrix, SymmGroup> t = mps[site].normalize_left(DefaultSolver());
                    if (site < L-1)
                        mps[site+1].multiply_from_left(t);
                }
                
                this->boundary_left_step(mpo, site); // creating left_[site+1]
                if (site != L-1) {
                    Storage::drop(right_[site+1]);
                    Storage::evict(left_[site]);
                }
            } else if (lr == -1) {
                if (site > 0) {
                    maquis::cout << "Growing, alpha = " << alpha << std::endl;
                    // Invalid read occurs after this!\n
                    trunc = mps.grow_r2l_sweep(mpo[site], left_[site], right_[site+1], site, alpha, cutoff, Mmax);
                } else {
                    block_matrix<Matrix, SymmGroup> t = mps[site].normalize_right(DefaultSolver());
                    if (site > 0)
                        mps[site-1].multiply_from_right(t);
                }
                this->boundary_right_step(mpo, site); // creating right_[site]
                if (site > 0) {
                    Storage::drop(left_[site]);
                    Storage::evict(right_[site+1]);
                }
            }
            if (poverlap.is_defined())
                poverlap.update(mps, site, lr);
            iteration_results_["BondDimension"]   << trunc.bond_dimension;
            iteration_results_["TruncatedWeight"] << trunc.truncated_weight;
            iteration_results_["SmallestEV"]      << trunc.smallest_ev;
            boost::chrono::high_resolution_clock::time_point sweep_then = boost::chrono::high_resolution_clock::now();
            double elapsed = boost::chrono::duration<double>(sweep_then - sweep_now).count();
            maquis::cout << "Sweep has been running for " << elapsed << " seconds." << std::endl;
            if (stop_callback())
                throw dmrg::time_limit(sweep, _site+1);
        }
        initial_site = -1;
    }
    
private:
    int initial_site;
};

#endif

