/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *
 *****************************************************************************/

#ifndef OPTIMIZE_H
#define OPTIMIZE_H

#include <boost/random.hpp>
#ifndef WIN32
#include <sys/time.h>
#define HAVE_GETTIMEOFDAY
#endif

#include "utils/sizeof.h"

#include "ietl_lanczos_solver.h"
#include "ietl_jacobi_davidson.h"

#include "dmrg/utils/BaseParameters.h"
#include "dmrg/utils/logger.h"
#include "dmrg/utils/stream_storage.h"

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

#ifdef HAVE_GETTIMEOFDAY
#define BEGIN_TIMING(name) \
gettimeofday(&now, NULL);
#define END_TIMING(name) \
gettimeofday(&then, NULL); \
maquis::cout << "Time elapsed in " << name << ": " << then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec) << std::endl;
#else
#define BEGIN_TIMING(name)
#define END_TIMING(name)
#endif

inline double log_interpolate(double y0, double y1, int N, int i)
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
class optimizer_base
{
public:
    optimizer_base(MPS<Matrix, SymmGroup> const & mps_,
             MPO<Matrix, SymmGroup> const & mpo_,
             BaseParameters & parms_,
             StorageMaster & sm)
    : mps(mps_)
    , mpo(mpo_)
    , mpo_orig(mpo_)
    , parms(parms_)
    , storage_master(sm)
    {
//        mps.normalize_right();
        mps.canonize(0);
        init_left_right(mpo, 0);
        maquis::cout << "Done init_left_right" << std::endl;
    }
    
    virtual int sweep(int sweep, Logger & iteration_log,
               OptimizeDirection d = Both,
               int resume_at = -1,
               int max_secs = -1) = 0;
    
    MPS<Matrix, SymmGroup> get_current_mps() const { return mps; }

protected:

    inline void boundary_left_step(MPO<Matrix, SymmGroup> const & mpo, int site)
    {
        MPSTensor<Matrix, SymmGroup> bkp = mps[site];
        Boundary<Matrix, SymmGroup> left = contraction::overlap_mpo_left_step(mps[site], bkp, left_[site], mpo[site]);
        left_[site+1] = left;        
    }
    
    inline void boundary_right_step(MPO<Matrix, SymmGroup> const & mpo, int site)
    {
        MPSTensor<Matrix, SymmGroup> bkp = mps[site];
        Boundary<Matrix, SymmGroup> right = contraction::overlap_mpo_right_step(mps[site], bkp, right_[site+1], mpo[site]);
        right_[site] = right;
    }

    void init_left_right(MPO<Matrix, SymmGroup> const & mpo, int site)
    {
        static Timer timer2("init_left_right");
        timer2.begin();
        std::size_t L = mps.length();
        
        left_.resize(mpo.length()+1);
        right_.resize(mpo.length()+1);
        
        right_stores_.resize(L+1, storage_master.child());
        left_stores_.resize(L+1, storage_master.child());
        
        Boundary<Matrix, SymmGroup> left = mps.left_boundary();
        storage::reset(left_stores_[0]);
        left_[0] = left;
        
        for (int i = 0; i < site; ++i) {
            storage::reset(left_stores_[i+1]);
            boundary_left_step(mpo, i);
            storage::store(left_[i], left_stores_[i]);
        }
        storage::store(left_[site], left_stores_[site]);
        
        
        Boundary<Matrix, SymmGroup> right = mps.right_boundary();
        storage::reset(right_stores_[L]);
        right_[L] = right;
                
        for(int i = L-1; i >= site; --i) {
            storage::reset(right_stores_[i]);
            boundary_right_step(mpo, i);
            storage::store(right_[i+1], right_stores_[i+1]);
        }
        storage::store(right_[site], right_stores_[site]);
        
        timer2.end();
    }
    
    MPS<Matrix, SymmGroup> mps;
    MPO<Matrix, SymmGroup> mpo, mpo_orig;
    
    BaseParameters & parms;
    std::vector<Boundary<Matrix, SymmGroup> > left_, right_;
    std::vector<typename StorageMaster::Storage> left_stores_, right_stores_;
    StorageMaster & storage_master;
};

#include "ss_optimize.hpp"
#include "ts_optimize.hpp"

#endif
