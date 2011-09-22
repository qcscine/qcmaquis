/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef SS_OPTIMIZE_H
#define SS_OPTIMIZE_H

#include <boost/random.hpp>
#include <sys/time.h>

#include "utils/zout.hpp"
#include "utils/sizeof.h"

#include "ietl_lanczos_solver.h"
#include "ietl_jacobi_davidson.h"

#include "utils/BaseParameters.h"
#include "utils/logger.h"
#include "utils/stream_storage.h"

template<class Matrix, class SymmGroup>
struct SiteProblem
{
    SiteProblem(MPSTensor<Matrix, SymmGroup> const & ket_tensor_,
                Boundary<Matrix, SymmGroup> const & left_,
                Boundary<Matrix, SymmGroup> const & right_,
                MPOTensor<Matrix, SymmGroup> const & mpo_)
    : ket_tensor(ket_tensor_)
    , left(left_)
    , right(right_)
    , mpo(mpo_) { }
    
    MPSTensor<Matrix, SymmGroup> const & ket_tensor;
    Boundary<Matrix, SymmGroup> const & left;
    Boundary<Matrix, SymmGroup> const & right;
    MPOTensor<Matrix, SymmGroup> const & mpo;
};

#define BEGIN_TIMING(name) \
gettimeofday(&now, NULL);
#define END_TIMING(name) \
gettimeofday(&then, NULL); \
zout << "Time elapsed in " << name << ": " << then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec) << endl;

double log_interpolate(double y0, double y1, int N, int i)
{
    if (N < 2)
        return y1;
    if (y0 == 0)
        return 0;
    double x = log(y1/y0)/(N-1);
    return y0*exp(x*i);
}

enum OptimizeDirection { Both, LeftOnly, RightOnly };

template<class Matrix, class SymmGroup, class StorageMaster>
class ss_optimize
{
public:
    ss_optimize(MPS<Matrix, SymmGroup> const & mps_,
                MPO<Matrix, SymmGroup> const & mpo_,
                BaseParameters & parms_,
                StorageMaster & sm)
    : mps(mps_)
    , mpo(mpo_)
    , parms(parms_)
    , storage_master(sm)
    {
//        mps.normalize_right();
        mps.canonize(0);
        init_left_right(mpo);
        zout << "Done init_left_right" << endl;
    }
    
    int sweep(int sweep, Logger & iteration_log,
               OptimizeDirection d = Both,
               int resume_at = -1,
               int max_secs = -1)
    {
        timeval sweep_now, sweep_then;
        gettimeofday(&sweep_now, NULL);
        
        static Timer
        t_io("sweep_io"),
        t_solver("sweep_solver"),
        t_grow("sweep_grow");
        
        std::size_t L = mps.length();
        
        if (resume_at != -1)
        {
            int site;
            if (resume_at < L)
                site = resume_at;
            else
                site = 2*L-resume_at-1;
            mps.canonize(site);
            init_left_right(mpo);
        }

        storage::prefetch(left_[0], left_stores_[0]);
        storage::prefetch(right_[1], right_stores_[1]);
        
        zout << mps.description() << endl;
        for (int _site = (resume_at == -1 ? 0 : resume_at);
             _site < 2*L; ++_site) {
            Timer iteration_t("Iteration took");
            iteration_t.begin();
            
            int site, lr;
            if (_site < L) {
                site = _site;
                lr = 1;
            } else {
                site = 2*L-_site-1;
                lr = -1;
            }
            
            zout << "Sweep " << sweep << ", optimizing site " << site << endl;
//            storage_master.print_size();
            
//            mps[site].make_left_paired();
            
            t_io.begin();
            
            storage::load(left_[site], left_stores_[site]);
            storage::load(right_[site+1], right_stores_[site+1]);
            
            if (lr == +1) {
                storage::prefetch(left_[site+1], left_stores_[site+1]);
                if (site+2 < right_.size())
                    storage::prefetch(right_[site+2], right_stores_[site+2]);
            } else {
                storage::prefetch(right_[site], right_stores_[site]);
                if (site > 1)
                    storage::prefetch(left_[site-1], left_stores_[site-1]);
            }
            
//            cout << "My size: " << endl;
//            cout << "  left_: " << utils::size_of(left_.begin(), left_.end())/1024.0/1024 << endl;
//            cout << "  right_: " << utils::size_of(right_.begin(), right_.end())/1024.0/1024 << endl;
//            cout << "  MPS: " << utils::size_of(mps.begin(), mps.end())/1024.0/1024 << endl;
//            cout << "  MPS[i]: " << utils::size_of(mps[site])/1024.0/1024 << endl;
            
            t_io.end();
            t_solver.begin();
            SiteProblem<Matrix, SymmGroup> sp(mps[site], left_[site], right_[site+1], mpo[site]);
            
            timeval now, then;

//            { // LAUSANNE
//                MPSTensor<Matrix, SymmGroup> vec1 = mps[site], vec2;
//                vec1.make_left_paired(); vec2.make_left_paired();
//                cout << vec1 << " " << vec2 << endl;
//                ietl::mult(sp, vec1, vec2);
//                vec1.make_left_paired(); vec2.make_left_paired();
//                cout << vec1 << " " << vec2 << endl;
//                cout << "Initial energy guess " << ietl::dot(vec1, vec2) << endl;;
//            } // LAUSANNE
            
            std::pair<double, MPSTensor<Matrix, SymmGroup> > res;
           
            if (d == Both ||
                (d == LeftOnly && lr == -1) ||
                (d == RightOnly && lr == +1))
            {
                if (parms.get<std::string>("eigensolver") == std::string("IETL")) {
                    BEGIN_TIMING("IETL")
                    res = solve_ietl_lanczos(sp, mps[site], parms);
                    END_TIMING("IETL")
                } else if (parms.get<std::string>("eigensolver") == std::string("IETL_JCD")) {
                    BEGIN_TIMING("JCD")
                    res = solve_ietl_jcd(sp, mps[site], parms);
                    END_TIMING("JCD")
                } else if (parms.get<std::string>("eigensolver") == std::string("IETL_NEW_JCD")) {
                    BEGIN_TIMING("JD")
                    res = solve_ietl_new_jd(sp, mps[site], parms);
                    END_TIMING("JD")
                } else {
                    throw std::runtime_error("I don't know this eigensolver.");
                }
 
//                {
//                    ietl::mult(sp, mps[site], res.second);
//                    res.first = ietl::dot(res.second, mps[site]);
//                    res.second = mps[site];
//                }
                
                mps[site] = res.second;
            }
            
            t_solver.end();
            
            zout << "Energy " << lr << " " << res.first << endl;
            
            iteration_log << make_log("Energy", res.first);
            
            double alpha;
//            if (sweep < parms.get<int>("ngrowsweeps"))
//                alpha = parms.get<double>("alpha_initial");
//            else
//                alpha = log_interpolate(parms.get<double>("alpha_initial"), parms.get<double>("alpha_final"),
//                                        parms.get<int>("nsweeps")-parms.get<int>("ngrowsweeps"),
//                                        sweep-parms.get<int>("ngrowsweeps"));
            int ngs = parms.get<int>("ngrowsweeps"), nms = parms.get<int>("nmainsweeps");
            if (sweep < ngs)
                alpha = parms.get<double>("alpha_initial");
            else if (sweep < ngs + nms)
                alpha = parms.get<double>("alpha_main");
            else
                alpha = parms.get<double>("alpha_final");
            
            
            double cutoff;
            if (sweep >= parms.get<int>("ngrowsweeps"))
                cutoff = parms.get<double>("truncation_final");
            else
                cutoff = log_interpolate(parms.get<double>("truncation_initial"), parms.get<double>("truncation_final"), parms.get<int>("ngrowsweeps"), sweep);
            
            std::size_t Mmax;
            if (parms.is_set("sweep_bond_dimensions")) {
                std::vector<std::size_t> ssizes = parms.get<std::vector<std::size_t> >("sweep_bond_dimensions");
                if (sweep >= ssizes.size())
                    Mmax = *ssizes.rbegin();
                else
                    Mmax = ssizes[sweep];
            } else
                Mmax = parms.get<std::size_t>("max_bond_dimension");
            
            std::pair<std::size_t, double> trunc;
            
            t_grow.begin();
                
            if (lr == +1) {
                if (site < L-1) {
                    zout << "Growing, alpha = " << alpha << endl;
                    mps.grow_l2r_sweep(mpo[site], left_[site], right_[site+1],
                                       site, alpha, cutoff, Mmax, iteration_log);
                } else {
                    block_matrix<Matrix, SymmGroup> t = mps[site].normalize_left(SVD);
                    if (site < L-1)
                        mps[site+1].multiply_from_left(t);
                }
                
                storage::load(left_[site+1], left_stores_[site+1]);
                
                MPSTensor<Matrix, SymmGroup> bkp = mps[site];
                left_[site+1] = contraction::overlap_mpo_left_step(mps[site], bkp,
                                                                   left_[site], mpo[site]);
                
                storage::store(left_[site], left_stores_[site]);
                storage::store(right_[site+1], right_stores_[site+1]);
            } else if (lr == -1) {
                if (site > 0) {
                    zout << "Growing, alpha = " << alpha << endl;
                    // Invalid read occurs after this!\n
                    mps.grow_r2l_sweep(mpo[site], left_[site], right_[site+1],
                                       site, alpha, cutoff, Mmax, iteration_log);
                } else {
                    block_matrix<Matrix, SymmGroup> t = mps[site].normalize_right(SVD);
                    if (site > 0)
                        mps[site-1].multiply_from_right(t);
                }
                
                storage::load(right_[site], right_stores_[site]);
                
                MPSTensor<Matrix, SymmGroup> bkp = mps[site];
                right_[site] = contraction::overlap_mpo_right_step(mps[site], bkp,
                                                                   right_[site+1], mpo[site]);
                
                storage::store(left_[site], left_stores_[site]);
                storage::store(right_[site+1], right_stores_[site+1]);
            }
            
            t_grow.end();
            
            iteration_t.end();
            
            gettimeofday(&sweep_then, NULL);
            double elapsed = sweep_then.tv_sec-sweep_now.tv_sec + 1e-6 * (sweep_then.tv_usec-sweep_now.tv_usec);
            zout << "Sweep has been running for " << elapsed << " seconds." << endl;
            if (max_secs != -1 && elapsed > max_secs && _site+1<2*L) {
                return _site+1;
            }
            else
                zout << max_secs - elapsed << " seconds left." << endl;
        }
        
        return -1;
    }
    
    MPS<Matrix, SymmGroup> get_current_mps() const { return mps; }
    
private:
    void init_left_right(MPO<Matrix, SymmGroup> const & mpo)
    {
        static Timer timer2("init_left_right");
        timer2.begin();
        std::size_t L = mps.length();
        
        left_.resize(mpo.length()+1);
        right_.resize(mpo.length()+1);
        
        right_stores_.resize(L+1, storage_master.child());
        left_stores_.resize(L+1, storage_master.child());
        
        Boundary<Matrix, SymmGroup> left = mps.left_boundary();
        left_[0] = left;

        storage::reset(left_stores_[0]);
        storage::store(left_[0], left_stores_[0]);
        
        for (int i = 0; i < L; ++i) {
            MPSTensor<Matrix, SymmGroup> bkp = mps[i];
            left = contraction::overlap_mpo_left_step(mps[i], bkp, left, mpo[i]);
            left_[i+1] = left;
            storage::reset(left_stores_[i+1]);
            storage::store(left_[i+1], left_stores_[i+1]);
        }
        
        Boundary<Matrix, SymmGroup> right = mps.right_boundary();
        right_[L] = right;

        storage::reset(right_stores_[L]);
        storage::store(right_[L], right_stores_[L]);
        
        for(int i = L-1; i >= 0; --i) {
            MPSTensor<Matrix, SymmGroup> bkp = mps[i];
            right = contraction::overlap_mpo_right_step(mps[i], bkp, right, mpo[i]);
            right_[i] = right;

            storage::reset(right_stores_[i]);
            storage::store(right_[i], right_stores_[i]);
        }
        timer2.end();
    }
    
    MPS<Matrix, SymmGroup> mps;
    MPO<Matrix, SymmGroup> mpo;
    
    BaseParameters & parms;
    std::vector<Boundary<Matrix, SymmGroup> > left_, right_;
    std::vector<typename StorageMaster::Storage> left_stores_, right_stores_;
    StorageMaster & storage_master;
};

#endif

