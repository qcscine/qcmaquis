/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef SS_OPTIMIZE_H
#define SS_OPTIMIZE_H

#include "dmrg/mp_tensors/mpo_ops.h"
#include "dmrg/optimize/optimize.h"


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
    using base::stop_callback;

    ss_optimize(MPS<Matrix, SymmGroup> & mps_,
                MPO<Matrix, SymmGroup> const & mpo_,
                BaseParameters & parms_,
                boost::function<bool ()> stop_callback_,
                const Lattice& lat,
                int initial_site_ = 0)
    : base(mps_, mpo_, parms_, stop_callback_, to_site(mps_.length(), initial_site_)), lat_(lat)
    , initial_site((initial_site_ < 0) ? 0 : initial_site_)
    { }

    inline int to_site(const int L, const int i) const
    {
        if (i < 0) return 0;
        /// i, or (L-1) - (i - L)
        return (i < L) ? i : 2*L - 1 - i;
    }

    void sweep(int sweep, OptimizeDirection d = Both)
    {
        std::chrono::high_resolution_clock::time_point sweep_now = std::chrono::high_resolution_clock::now();

        iteration_results_.clear();

        std::size_t L = mps.length();

        int _site = 0, site = 0;
        if (initial_site != -1) {
            _site = initial_site;
            site = to_site(L, _site);
        }

        Storage::prefetch(left_[site]);
        Storage::prefetch(right_[site+1]);

        for (; _site < 2*L; ++_site) {

            int lr = (_site < L) ? +1 : -1;
            site = to_site(L, _site);

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

            if (lr == +1 && site+2 <= L)
                Storage::prefetch(right_[site+2]);
            if (lr == -1 && site > 0)
                Storage::prefetch(left_[site-1]);

            assert( left_[site].reasonable() );    // in case something went wrong
            assert( right_[site+1].reasonable() ); // in case something went wrong

            //SiteProblem<Matrix, SymmGroup> sp(mps[site], left_[site], right_[site+1], mpo[site]);

            std::chrono::high_resolution_clock::time_point now, then;

            std::pair<typename MPSTensor<Matrix, SymmGroup>::value_type, MPSTensor<Matrix, SymmGroup> > res;
            SiteProblem<Matrix, SymmGroup> sp(left_[site], right_[site+1], mpo[site]);

            /// Compute orthogonal vectors
            std::vector<MPSTensor<Matrix, SymmGroup> > ortho_vecs(base::northo);
            for (int n = 0; n < base::northo; ++n)
                ortho_vecs[n] = contraction::site_ortho_boundaries(mps[site], base::ortho_mps[n][site],
                                                                    base::ortho_left_[n][site], base::ortho_right_[n][site+1]);

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

            {
                int prec = maquis::cout.precision();
                maquis::cout.precision(15);
                maquis::cout << "Energy " << lr << " " << res.first + mpo.getCoreEnergy()<< std::endl;
                maquis::cout.precision(prec);
            }

            iteration_results_["Energy"] << res.first + mpo.getCoreEnergy();

            double alpha;
            int ngs = parms.template get<int>("ngrowsweeps"), nms = parms.template get<int>("nmainsweeps");
            if (sweep < ngs)
                alpha = parms.template get<double>("alpha_initial");
            else if (sweep < ngs + nms)
                alpha = parms.template get<double>("alpha_main");
            else
                alpha = parms.template get<double>("alpha_final");

            double cutoff = this->get_cutoff(sweep);
            std::size_t Mmax;
            if (!parms.is_set("PreBO_MaxBondDimVector"))
                Mmax = this->get_Mmax(sweep);
            else {
                Mmax = lat_.template get_prop<size_t>("Mmax", {lat_.template get_prop<int>("type", {site}) });
                std::cout << "Mmax is set to " << Mmax << std::endl;
            }
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

                this->boundary_left_step(mpo, site); // creating left_[site+1]
                if (site != L-1) {
                    Storage::drop(right_[site+1]);
                    Storage::StoreToFile(left_[site]);
                }
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

                this->boundary_right_step(mpo, site); // creating right_[site]
                if (site > 0) {
                    Storage::drop(left_[site]);
                    Storage::StoreToFile(right_[site+1]);
                }
            }

            iteration_results_["BondDimension"]   << trunc.bond_dimension;
            iteration_results_["TruncatedWeight"] << trunc.truncated_weight;
            iteration_results_["SmallestEV"]      << trunc.smallest_ev;

            std::chrono::high_resolution_clock::time_point sweep_then = std::chrono::high_resolution_clock::now();
            double elapsed = std::chrono::duration<double>(sweep_then - sweep_now).count();
            maquis::cout << "Sweep has been running for " << elapsed << " seconds." << std::endl;

            if (stop_callback())
                throw dmrg::time_limit(sweep, _site+1);
        }
        initial_site = -1;
    }

private:
    int initial_site;
    const Lattice &lat_;
};

#endif

