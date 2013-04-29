/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch> 
 *			      Sebastian Keller <sebkelle@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef TS_OPTIMIZE_H
#define TS_OPTIMIZE_H

#include <boost/tuple/tuple.hpp>

#include "dmrg/mp_tensors/twositetensor.h"
#include "dmrg/mp_tensors/mpo_ops.h"

template<class Matrix, class SymmGroup, class StorageMaster>
class ts_optimize : public optimizer_base<Matrix, SymmGroup, StorageMaster>
{
public:

    typedef optimizer_base<Matrix, SymmGroup, StorageMaster> base;
    using base::mpo;
    using base::mpo_orig;
    using base::mps;
    using base::left_;
    using base::left_stores_;
    using base::right_;
    using base::right_stores_;
    using base::parms;

    ts_optimize(MPS<Matrix, SymmGroup> const & mps_,
                MPO<Matrix, SymmGroup> const & mpo_,
                MPO<Matrix, SymmGroup> const & ts_mpo_,
                BaseParameters & parms_,
                StorageMaster & sm)
    : ts_cache_mpo(ts_mpo_), base(mps_, mpo_, parms_, sm) { }

    int sweep(int sweep, Logger & iteration_log,
               OptimizeDirection d = Both,
               int resume_at = -1,
               int max_secs = -1)
    {
        // Sebastian: Is this still necessary?
        mpo = mpo_orig;

    	timeval sweep_now, sweep_then;
    	gettimeofday(&sweep_now, NULL);

        
        std::size_t L = mps.length();

        if (resume_at != -1)
        {
            int site;
            if (resume_at < L)
                site = resume_at;
            else
                site = 2*L-resume_at-2;
            mps.canonize(site);
            this->init_left_right(mpo, site);
        }

        storage::prefetch(left_[0], left_stores_[0]);
        storage::prefetch(right_[2], right_stores_[2]);

#ifndef NDEBUG
    	maquis::cout << mps.description() << std::endl;
#endif
        for (int _site = (resume_at == -1 ? 0 : resume_at); _site < 2*L-2; ++_site) {
	/* (0,1), (1,2), ... , (L-1,L), (L-1,L), (L-2, L-1), ... , (0,1)
	    | |                        |
       site 1 |                        |
	      |         left to right  | right to left, lr = -1
	      site 2                   |                               */

            int site, lr, site1, site2;
            if (_site < L-1) {
                site = _site;
                lr = 1;
        		site1 = site;
        		site2 = site+1;
            } else {
                site = 2*L-_site-2;
                lr = -1;
        		site1 = site-1;
        		site2 = site;
            }
            
    	    maquis::cout << std::endl;
            maquis::cout << "Sweep " << sweep << ", optimizing sites " << site1 << " and " << site2 << std::endl;

            if (parms.template get<bool>("beta_mode")) {
                if (sweep == 0 && lr == 1) {
                    mpo = zero_after(mpo_orig, 0);
                    if (site == 0)
                        this->init_left_right(mpo, 0);
                } else if (sweep == 0 && lr == -1 && site == L-1) {
                    mpo = mpo_orig;
                    // this is not needed, but when enabled is showing a bug in StreamStorage
                    // todo: need more investigation
                    //this->init_left_right(mpo, site);
                }
            }
        
        
            storage::load(left_[site1], left_stores_[site1]);
            storage::load(right_[site2+1], right_stores_[site2+1]);

            if (lr == +1) {
                if (site2+2 < right_.size())
                    storage::prefetch(right_[site2+2], right_stores_[site2+2]);
            } else {
                if (site1 > 1)
                    storage::prefetch(left_[site1-1], left_stores_[site1-1]);
            }


    	    timeval now, then;

    	    // Create TwoSite objects
    	    TwoSiteTensor<Matrix, SymmGroup> tst(mps[site1], mps[site2]);
    	    MPSTensor<Matrix, SymmGroup> twin_mps = tst.make_mps();
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
                if (parms.template get<std::string>("eigensolver") == std::string("IETL")) {
            	    BEGIN_TIMING("IETL")
                    res = solve_ietl_lanczos(sp, twin_mps, parms);
            	    END_TIMING("IETL")
                } else if (parms.template get<std::string>("eigensolver") == std::string("IETL_JCD")) {
            	    BEGIN_TIMING("JCD")
                    res = solve_ietl_jcd(sp, twin_mps, parms, ortho_vecs);
            	    END_TIMING("JCD")
                } else {
                    throw std::runtime_error("I don't know this eigensolver.");
                }

        		tst << res.second;
            }

#ifndef NDEBUG
            // Caution: this is an O(L) operation, so it really should be done only in debug mode
            for (int n = 0; n < base::northo; ++n)
                maquis::cout << "MPS overlap: " << overlap(mps, base::ortho_mps[n]) << std::endl;
#endif

            maquis::cout << "Energy " << lr << " " << res.first << std::endl;
            iteration_log << make_log("Energy", res.first);
            
            double cutoff = this->get_cutoff(sweep);
            std::size_t Mmax = this->get_Mmax(sweep);
            
    	    if (lr == +1)
    	    {
        		// Write back result from optimization
        		boost::tie(mps[site1], mps[site2]) = tst.split_mps_l2r(Mmax, cutoff, &iteration_log);

        		block_matrix<Matrix, SymmGroup> t;
		
        		//t = mps[site1].normalize_left(DefaultSolver());
        		//mps[site2].multiply_from_left(t);
        		//mps[site2].divide_by_scalar(mps[site2].scalar_norm());	

        		t = mps[site2].normalize_left(DefaultSolver());
                // MD: DEBUGGING OUTPUT
                maquis::cout << "Propagating t with norm " << t.norm() << std::endl;
        		if (site2 < L-1) mps[site2+1].multiply_from_left(t);

                storage::reset(left_stores_[site2]); // left_stores_[site2] is outdated
                this->boundary_left_step(mpo, site1); // creating left_[site2]
    	    }
    	    if (lr == -1){
        		// Write back result from optimization
        		boost::tie(mps[site1], mps[site2]) = tst.split_mps_r2l(Mmax, cutoff, &iteration_log);

        		block_matrix<Matrix, SymmGroup> t;

        		//t = mps[site2].normalize_right(DefaultSolver());
        		//mps[site1].multiply_from_right(t);
        		//mps[site1].divide_by_scalar(mps[site1].scalar_norm());	

        		t = mps[site1].normalize_right(DefaultSolver());
                // MD: DEBUGGING OUTPUT
                maquis::cout << "Propagating t with norm " << t.norm() << std::endl;
        		if (site1 > 0) mps[site1-1].multiply_from_right(t);

                storage::reset(right_stores_[site2]); // right_stores_[site2] is outdated
                this->boundary_right_step(mpo, site2); // creating right_[site2]
    	    }
            
            if (_site != L-1)
            { 
                storage::store(left_[site1], left_stores_[site1]); // store currently used boundary
                storage::store(right_[site2+1], right_stores_[site2+1]); // store currently used boundary
            }
            
            gettimeofday(&sweep_then, NULL);
            double elapsed = sweep_then.tv_sec-sweep_now.tv_sec + 1e-6 * (sweep_then.tv_usec-sweep_now.tv_usec);
            maquis::cout << "Sweep has been running for " << elapsed << " seconds." << std::endl;
            if (max_secs != -1 && elapsed > max_secs && _site+1<2*L) {
                return _site+1;
            } else
                maquis::cout << max_secs - elapsed << " seconds left." << std::endl;
    	} // for sites

        return -1;
    } // sweep

private:
    MPO<Matrix, SymmGroup> const & ts_cache_mpo;
};

#endif
