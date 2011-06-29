/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef TIME_EVOLVE_H
#define TIME_EVOLVE_H

#include <boost/random.hpp>
#include <sys/time.h>

#include "utils/zout.hpp"
#include "utils/sizeof.h"

#include "ietl_lanczos_solver.h"
#include "ietl_jacobi_davidson.h"
#ifdef HAVE_ARPACK
#include "arpackpp_solver.h"
#endif

#include "utils/DmrgParameters.h"
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

template<class Matrix, class SymmGroup, class StorageMaster>
class time_evolve
{
public:
    time_evolve(MPS<Matrix, SymmGroup> const & mps_,
                std::vector<MPO<Matrix, SymmGroup> > const & mpos_,
                BaseParameters & parms_,
                StorageMaster & sm)
    : mps(mps_)
    , mpos(mpos_)
    , parms(parms_)
    , storage_master(sm)
    {
        //        mps.normalize_right();
        mps.canonize(0);
        cout << "Done init_left_right" << endl;
    }
    
    int sweep(int sweep, Logger & iteration_log,
              int which = 0,
              int resume_at = -1,
              int max_secs = -1)
    {
        timeval sweep_now, sweep_then;
        gettimeofday(&sweep_now, NULL);
        
        static Timer
        t_io("sweep_io"),
        t_solver("sweep_solver"),
        t_grow("sweep_grow");
        
#ifdef MPI_PARALLEL
        ambient::playout();
#endif
        
        std::size_t L = mps.length();
        
        if (resume_at != -1)
        {
            int site;
            if (resume_at < L)
                site = resume_at;
            else
                site = 2*L-resume_at-1;
            mps.canonize(site);
        }
        
        for (; which<mpos.size(); ++which)
        {
            
            storage::prefetch(left_[0], left_stores_[0]);
            storage::prefetch(right_[1], right_stores_[1]);
            
            zout << mps.description() << endl;
            for (int _site = (resume_at == -1 ? 0 : resume_at);
                 _site < 2*L; ++_site) {
                
                init_left_right(mpos[which]);

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
                SiteProblem<Matrix, SymmGroup> sp(mps[site], left_[site], right_[site+1], mpos[which][site]);
                
                timeval now, then;
                
                std::pair<double, MPSTensor<Matrix, SymmGroup> > res;
                
                
                MPSTensor<Matrix, SymmGroup> t;
                ietl::mult(sp, mps[site], t);
                mps[site] = t;
                
                t_solver.end();
                
                                
                iteration_t.end();
                
            }
            
            mps.canonize(0);

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
    std::vector<MPO<Matrix, SymmGroup> > mpos;
    
    BaseParameters & parms;
    std::vector<Boundary<Matrix, SymmGroup> > left_, right_;
    std::vector<typename StorageMaster::Storage> left_stores_, right_stores_;
    StorageMaster & storage_master;
};

#endif

