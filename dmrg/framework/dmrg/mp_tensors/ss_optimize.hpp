/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef SS_OPTIMIZE_H
#define SS_OPTIMIZE_H


template<class Matrix, class SymmGroup, class StorageMaster>
class ss_optimize : public optimizer_base<Matrix, SymmGroup, StorageMaster>
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

    ss_optimize(MPS<Matrix, SymmGroup> const & mps_,
                MPO<Matrix, SymmGroup> const & mpo_,
                BaseParameters & parms_,
                StorageMaster & sm)
    : base(mps_, mpo_, parms_, sm) { }
    
    int sweep(int sweep, Logger & iteration_log,
               OptimizeDirection d = Both,
               int resume_at = -1,
               int max_secs = -1)
    {
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
                site = 2*L-resume_at-1;
            mps.canonize(site);
            this->init_left_right(mpo, site);
        }

//        if (parms.template <bool>("beta_mode") && sweep == 0 && resume_at < L) {
//            int site = (resume_at == -1) ? 0 : resume_at;
//            mpo = zero_after(mpo_orig, site+2);
//            mps.canonize(site);
//            this->init_left_right(mpo, site);
//        }
        
        storage::prefetch(left_[0], left_stores_[0]);
        storage::prefetch(right_[1], right_stores_[1]);
        
#ifndef NDEBUG
        maquis::cout << mps.description() << std::endl;
#endif
        for (int _site = (resume_at == -1 ? 0 : resume_at);
             _site < 2*L; ++_site) {
            
            int site, lr;
            if (_site < L) {
                site = _site;
                lr = 1;
            } else {
                site = 2*L-_site-1;
                lr = -1;
            }
            
            maquis::cout << "Sweep " << sweep << ", optimizing site " << site << std::endl;
//            storage_master.print_size();
            
//            mps[site].make_left_paired();
            
            if (parms.template get<bool>("beta_mode")) {
                if (sweep == 0 && lr == 1) {
                    mpo = zero_after(mpo_orig, 0);
                    if (site == 0)
                        this->init_left_right(mpo, 0);
                } else if (sweep == 0 && lr == -1 && site == L-1) {
                    mpo = mpo_orig;
                    //this->init_left_right(mpo, site);
                }
            }
            
            
            storage::load(left_[site], left_stores_[site]);
            storage::load(right_[site+1], right_stores_[site+1]);
            
            if (lr == +1) {
                //storage::prefetch(left_[site+1], left_stores_[site+1]);
                if (site+2 < right_.size())
                    storage::prefetch(right_[site+2], right_stores_[site+2]);
            } else {
//                storage::prefetch(right_[site], right_stores_[site]);
                if (site > 1)
                    storage::prefetch(left_[site-1], left_stores_[site-1]);
            }
            
//            maquis::cout << "My size: " << std::endl;
//            maquis::cout << "  left_: " << utils::size_of(left_.begin(), left_.end())/1024.0/1024 << std::endl;
//            maquis::cout << "  right_: " << utils::size_of(right_.begin(), right_.end())/1024.0/1024 << std::endl;
//            maquis::cout << "  MPS: " << utils::size_of(mps.begin(), mps.end())/1024.0/1024 << std::endl;
//            maquis::cout << "  MPS[i]: " << utils::size_of(mps[site])/1024.0/1024 << std::endl;
            
            SiteProblem<Matrix, SymmGroup> sp(mps[site], left_[site], right_[site+1], mpo[site]);
            
            timeval now, then;

//            { // LAUSANNE
//                MPSTensor<Matrix, SymmGroup> vec1 = mps[site], vec2;
//                vec1.make_left_paired(); vec2.make_left_paired();
//                maquis::cout << vec1 << " " << vec2 << std::endl;
//                ietl::mult(sp, vec1, vec2);
//                vec1.make_left_paired(); vec2.make_left_paired();
//                maquis::cout << vec1 << " " << vec2 << std::endl;
//                maquis::cout << "Initial energy guess " << ietl::dot(vec1, vec2) << std::endl;
//            } // LAUSANNE
            
            std::pair<double, MPSTensor<Matrix, SymmGroup> > res;
           
            if (d == Both ||
                (d == LeftOnly && lr == -1) ||
                (d == RightOnly && lr == +1))
            {
                if (parms.template get<std::string>("eigensolver") == std::string("IETL")) {
                    BEGIN_TIMING("IETL")
                    res = solve_ietl_lanczos(sp, mps[site], parms);
                    END_TIMING("IETL")
                } else if (parms.template get<std::string>("eigensolver") == std::string("IETL_JCD")) {
                    BEGIN_TIMING("JCD")
                    res = solve_ietl_jcd(sp, mps[site], parms);
                    END_TIMING("JCD")
                } /* else if (parms.template <std::string>("eigensolver") == std::string("IETL_NEW_JCD")) {
                    BEGIN_TIMING("JD")
                    res = solve_ietl_new_jd(sp, mps[site], parms);
                    END_TIMING("JD")
                } */ else {
                    throw std::runtime_error("I don't know this eigensolver.");
                }
 
//                {
//                    ietl::mult(sp, mps[site], res.second);
//                    res.first = ietl::dot(res.second, mps[site]);
//                    res.second = mps[site];
//                }
                
                mps[site] = res.second;
            }
            
            
            maquis::cout << "Energy " << lr << " " << res.first << std::endl;
//            maquis::cout << "Energy check " << expval(mps, mpo) << std::endl;
            
            iteration_log << make_log("Energy", res.first);
            
            double alpha;
//            if (sweep < parms.template <int>("ngrowsweeps"))
//                alpha = parms.template <double>("alpha_initial");
//            else
//                alpha = log_interpolate(parms.template <double>("alpha_initial"), parms.template <double>("alpha_final"),
//                                        parms.template <int>("nsweeps")-parms.template <int>("ngrowsweeps"),
//                                        sweep-parms.template <int>("ngrowsweeps"));
            int ngs = parms.template get<int>("ngrowsweeps"), nms = parms.template get<int>("nmainsweeps");
            if (sweep < ngs)
                alpha = parms.template get<double>("alpha_initial");
            else if (sweep < ngs + nms)
                alpha = parms.template get<double>("alpha_main");
            else
                alpha = parms.template get<double>("alpha_final");
            
            
            double cutoff;
            if (sweep >= parms.template get<int>("ngrowsweeps"))
                cutoff = parms.template get<double>("truncation_final");
            else
                cutoff = log_interpolate(parms.template get<double>("truncation_initial"), parms.template get<double>("truncation_final"), parms.template get<int>("ngrowsweeps"), sweep);
            
            std::size_t Mmax;
            if (parms.is_set("sweep_bond_dimensions")) {
                std::vector<std::size_t> ssizes = parms.template get<std::vector<std::size_t> >("sweep_bond_dimensions");
                if (sweep >= ssizes.size())
                    Mmax = *ssizes.rbegin();
                else
                    Mmax = ssizes[sweep];
            } else
                Mmax = parms.template get<std::size_t>("max_bond_dimension");
            
            std::pair<std::size_t, double> trunc;
            
                
            if (lr == +1) {
                if (site < L-1) {
                    maquis::cout << "Growing, alpha = " << alpha << std::endl;
                    mps.grow_l2r_sweep(mpo[site], left_[site], right_[site+1],
                                       site, alpha, cutoff, Mmax, iteration_log);
                } else {
                    block_matrix<Matrix, SymmGroup> t = mps[site].normalize_left(SVD);
                    if (site < L-1)
                        mps[site+1].multiply_from_left(t);
                }
                
                
                storage::reset(left_stores_[site+1]); // left_stores_[site+1] is outdated
                this->boundary_left_step(mpo, site); // creating left_[site+1]
            } else if (lr == -1) {
                if (site > 0) {
                    maquis::cout << "Growing, alpha = " << alpha << std::endl;
                    // Invalid read occurs after this!\n
                    mps.grow_r2l_sweep(mpo[site], left_[site], right_[site+1],
                                       site, alpha, cutoff, Mmax, iteration_log);
                } else {
                    block_matrix<Matrix, SymmGroup> t = mps[site].normalize_right(SVD);
                    if (site > 0)
                        mps[site-1].multiply_from_right(t);
                }
                
                
                storage::reset(right_stores_[site]); // right_stores_[site] is outdated
                this->boundary_right_step(mpo, site); // creating right_[site]
            }
            
            if (! (lr == +1 && site == L-1))
            {
            	storage::store(left_[site], left_stores_[site]); // store currently used boundary
            	storage::store(right_[site+1], right_stores_[site+1]); // store currently used boundary
            }

            
            
            gettimeofday(&sweep_then, NULL);
            double elapsed = sweep_then.tv_sec-sweep_now.tv_sec + 1e-6 * (sweep_then.tv_usec-sweep_now.tv_usec);
            maquis::cout << "Sweep has been running for " << elapsed << " seconds." << std::endl;
            if (max_secs != -1 && elapsed > max_secs && _site+1<2*L) {
                return _site+1;
            }
            else
               maquis::cout << max_secs - elapsed << " seconds left." << std::endl;
        }
        
        return -1;
    }
    
};

#endif

