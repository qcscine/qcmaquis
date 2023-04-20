/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

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
    typedef typename Matrix::value_type value_type;

    typedef optimizer_base<Matrix, SymmGroup, Storage> base;
    using base::mpo;
    using base::mps;
    using base::left_;
    using base::right_;
    using base::parms;
    using base::iteration_results_;
    using base::stop_callback;

    ts_optimize(MPS<Matrix, SymmGroup> & mps_,
                MPO<Matrix, SymmGroup> const & mpo_,
                BaseParameters & parms_,
                boost::function<bool ()> stop_callback_,
                const Lattice& lat,
                int initial_site_ = 0)
    : base(mps_, mpo_, parms_, stop_callback_, to_site(mps_.length(), initial_site_)), lat_(lat)
    , initial_site((initial_site_ < 0) ? 0 : initial_site_)
    {
        parallel::guard::serial guard;
        make_ts_cache_mpo(mpo, ts_cache_mpo, mps);
    }

    inline int to_site(const int L, const int i) const
    {
        if (i < 0) return 0;
        /// i, or (L-1) - (i - (L-1))
        return (i < L-1) ? i : 2*L - 2 - i;
    }
    void sweep(int sweep, OptimizeDirection d = Both)
    {
        std::chrono::high_resolution_clock::time_point sweep_now = std::chrono::high_resolution_clock::now();

        iteration_results_.clear();

        std::size_t L = mps.length();
        parallel::scheduler_balanced scheduler_mps(L);

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
            int lr, site1, site2;
            if (_site < L-1) {
                site = to_site(L, _site);
                lr = 1;
                site1 = site;
                site2 = site+1;
            } else {
                site = to_site(L, _site);
                lr = -1;
                site1 = site-1;
                site2 = site;
            }

            //if (lr == +1) mps.canonize(site1);
            //else          mps.canonize(site2);

            maquis::cout << std::endl;
            maquis::cout << "Sweep " << sweep << ", optimizing sites " << site1 << " and " << site2 << std::endl;

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


            std::chrono::high_resolution_clock::time_point now, then;

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

            std::pair<value_type, MPSTensor<Matrix, SymmGroup> > res;

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
                    res = solve_ietl_jcd(sp, twin_mps, parms, ortho_vecs);
                    END_TIMING("JCD")
                } else if (parms["eigensolver"] == std::string("IETL_DAVIDSON")) {
                    BEGIN_TIMING("DAVIDSON")
                    res = solve_ietl_davidson(sp, twin_mps, parms, ortho_vecs);
                    END_TIMING("DAVIDSON")
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
            std::size_t Mmax;
            if (!parms.is_set("PreBO_MaxBondDimVector"))
                Mmax = this->get_Mmax(sweep);
            else {
                auto m1 = lat_.template get_prop<size_t>("Mmax", {lat_.template get_prop<int>("type", {site1}) });
                auto m2 = lat_.template get_prop<size_t>("Mmax", {lat_.template get_prop<int>("type", {site2}) });
                Mmax = (m1>m2) ? m1 : m2;
                std::cout << "Mmax is set to " << Mmax << std::endl;
            }
            truncation_results trunc;

            if (lr == +1)
            {
                // Write back result from optimization
                BEGIN_TIMING("TRUNC")
                if (parms["twosite_truncation"] == "svd")
                    boost::tie(mps[site1], mps[site2], trunc) = tst.split_mps_l2r(Mmax, cutoff);
                else
                    boost::tie(mps[site1], mps[site2], trunc) = tst.predict_split_l2r(Mmax, cutoff, alpha, left_[site1], mpo[site1], true);
                END_TIMING("TRUNC")
                tst.clear();


                block_matrix<Matrix, SymmGroup> t;
                //t = mps[site1].normalize_left(DefaultSolver());
                //mps[site2].multiply_from_left(t);
                //mps[site2].divide_by_scalar(mps[site2].scalar_norm());
                t = mps[site2].leftNormalizeAndReturn(DefaultSolver());
                // MD: DEBUGGING OUTPUT
                maquis::cout << "Propagating t with norm " << t.norm() << std::endl;
                if (site2 < L-1)
                    mps[site2+1].multiply_from_left(t);

                if (site1 != L-2)
                    Storage::drop(right_[site2+1]);

                this->boundary_left_step(mpo, site1); // creating left_[site2]

                if (site1 != L-2){
                    Storage::StoreToFile(mps[site1]);
                    Storage::StoreToFile(left_[site1]);
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
                    boost::tie(mps[site1], mps[site2], trunc) = tst.predict_split_r2l(Mmax, cutoff, alpha, right_[site2+1], mpo[site2], true);
                END_TIMING("TRUNC")
                tst.clear();


                block_matrix<Matrix, SymmGroup> t;
                //t = mps[site2].normalize_right(DefaultSolver());
                //mps[site1].multiply_from_right(t);
                //mps[site1].divide_by_scalar(mps[site1].scalar_norm());
                t = mps[site1].rightNormalizeAndReturn(DefaultSolver());
                // MD: DEBUGGING OUTPUT
                maquis::cout << "Propagating t with norm " << t.norm() << std::endl;
                if (site1 > 0) mps[site1-1].multiply_from_right(t);

                if(site1 != 0)
                    Storage::drop(left_[site1]);

                this->boundary_right_step(mpo, site2); // creating right_[site2]

                if(site1 != 0){
                    Storage::StoreToFile(mps[site2]);
                    Storage::StoreToFile(right_[site2+1]);
                }
                { parallel::guard proc(scheduler_mps(site1)); storage::migrate(mps[site1]); }
                { parallel::guard proc(scheduler_mps(site2)); storage::migrate(mps[site2]); }
            }

            iteration_results_["BondDimension"]     << trunc.bond_dimension;
            iteration_results_["TruncatedWeight"]   << trunc.truncated_weight;
            iteration_results_["TruncatedFraction"] << trunc.truncated_fraction;
            iteration_results_["SmallestEV"]        << trunc.smallest_ev;

            parallel::meminfo();

            std::chrono::high_resolution_clock::time_point sweep_then = std::chrono::high_resolution_clock::now();
            double elapsed = std::chrono::duration<double>(sweep_then - sweep_now).count();
            maquis::cout << "Sweep has been running for " << elapsed << " seconds." << std::endl;

            if (stop_callback())
                throw dmrg::time_limit(sweep, _site+1);

        } // for sites
        initial_site = -1;
    } // sweep

private:
    int initial_site;
    const Lattice& lat_;
    MPO<Matrix, SymmGroup> ts_cache_mpo;
};

#endif
