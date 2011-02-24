#ifndef SS_OPTIMIZE_H
#define SS_OPTIMIZE_H

#include <boost/random.hpp>
#include <sys/time.h>

#include "utils/zout.hpp"
#include "utils/sizeof.h"

#include "ietl_lanczos_solver.h"
#include "ietl_jacobi_davidson.h"
#include "arpackpp_solver.h"

#include "utils/DmrgParameters.h"

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
    double x = log(y1/y0)/(N-1);
    return y0*exp(x*i);
}

template<class Matrix, class SymmGroup, class StorageMaster>
class ss_optimize
{
public:
    ss_optimize(MPS<Matrix, SymmGroup> & mps_,
                BaseParameters & parms_,
                StorageMaster & sm)
    : mps(mps_)
    , parms(parms_)
    , storage_master(sm)
    { }
    
    std::vector<double> sweep(MPO<Matrix, SymmGroup> const & mpo,
                              int sweep)
    {
        mps.normalize_right();
        mps.canonize(0);
        
        init_left_right(mpo);
        cerr << "Done init_left_right" << endl;
        
        std::vector<double> energies;
        
        std::size_t L = mps.length();
        
        prefetch(left_[0], left_stores_[0]);
        prefetch(right_[1], right_stores_[1]);
        
        zout << mps.description() << endl;
        for (int _site = 0; _site < 2*L; ++_site) {
            int site, lr;
            if (_site < L) {
                site = _site;
                lr = 1;
            } else {
                site = 2*L-_site-1;
                lr = -1;
            }
            
            zout << "Sweep " << sweep << ", optimizing site " << site << endl;
            storage_master.print_size();
            
            mps[site].make_left_paired();
            
            load(left_[site], left_stores_[site]);
            load(right_[site+1], right_stores_[site+1]);
            
            if (lr == +1) {
                prefetch(left_[site+1], left_stores_[site+1]);
                if (site+2 < right_.size())
                    prefetch(right_[site+2], right_stores_[site+2]);
            } else {
                prefetch(right_[site], right_stores_[site]);
                if (site > 1)
                    prefetch(left_[site-1], left_stores_[site-1]);
            }
            
            cout << "My size: " << endl;
            cout << "  left_: " << utils::size_of(left_.begin(), left_.end())/1024.0/1024 << endl;
            cout << "  right_: " << utils::size_of(right_.begin(), right_.end())/1024.0/1024 << endl;
            cout << "  MPS: " << utils::size_of(mps.begin(), mps.end())/1024.0/1024 << endl;
            cout << "  MPS[i]: " << utils::size_of(mps[site])/1024.0/1024 << endl;
            
            SiteProblem<Matrix, SymmGroup> sp(mps[site], left_[site], right_[site+1], mpo[site]);
//            sp.ket_tensor = mps[site];
//            sp.mpo = mpo[site];
//            sp.left = left_[site];
//            sp.right = right_[site+1];
            
            timeval now, then;
            
            std::pair<double, MPSTensor<Matrix, SymmGroup> > res;
           
            if (parms.get<std::string>("eigensolver") == std::string("IETL")) {
                BEGIN_TIMING("IETL")
                res = solve_ietl_lanczos(sp, mps[site], parms);
                END_TIMING("IETL")
            } else if (parms.get<std::string>("eigensolver") == std::string("ARPACK")) {
                BEGIN_TIMING("ARPACK")
                res = solve_arpackpp(sp, mps[site], parms);
                END_TIMING("ARPACK")
            } else if (parms.get<std::string>("eigensolver") == std::string("IETL_JCD")) {
                BEGIN_TIMING("JCD")
                res = solve_ietl_jcd(sp, mps[site], parms);
                END_TIMING("JCD")
            } else {
                throw std::runtime_error("I don't know this eigensolver.");
            }
            mps[site] = res.second;
            
            zout << "Energy " << lr << " " << res.first << endl;
            energies.push_back(res.first);
            
            double alpha;
            if (sweep < parms.get<int>("ngrowsweeps"))
                alpha = parms.get<double>("alpha_initial");
            else
                alpha = log_interpolate(parms.get<double>("alpha_initial"), parms.get<double>("alpha_final"),
                                        parms.get<int>("nsweeps")-parms.get<int>("ngrowsweeps"),
                                        sweep-parms.get<int>("ngrowsweeps"));
            double cutoff;
            if (sweep >= parms.get<int>("ngrowsweeps"))
                cutoff = parms.get<double>("truncation_final");
            else
                cutoff = log_interpolate(parms.get<double>("truncation_initial"), parms.get<double>("truncation_final"), parms.get<int>("ngrowsweeps"), sweep);
            std::size_t Mmax = parms.get<std::size_t>("max_bond_dimension");
            
            if (lr == +1) {
                if (site < L-1) {
                    zout << "Growing, alpha = " << alpha << endl;
                    mps.grow_l2r_sweep(mpo[site], left_[site], right_[site+1],
                                       site, alpha, cutoff, Mmax);
                }
                
                block_matrix<Matrix, SymmGroup> t = mps[site].normalize_left(SVD);
                if (site < L-1)
                    mps[site+1].multiply_from_left(t);
                
                load(left_[site+1], left_stores_[site+1]);
                
                MPSTensor<Matrix, SymmGroup> bkp = mps[site];
                left_[site+1] = contraction::overlap_mpo_left_step(mps[site], bkp,
                                                                   left_[site], mpo[site]);
                
                store(left_[site], left_stores_[site]);
                store(right_[site+1], right_stores_[site+1]);
            } else if (lr == -1) {
                if (site > 1) {
                    zout << "Growing, alpha = " << alpha << endl;
                    mps.grow_r2l_sweep(mpo[site], left_[site], right_[site+1],
                                       site, alpha, cutoff, Mmax);
                }
                
                block_matrix<Matrix, SymmGroup> t = mps[site].normalize_right(SVD);
                if (site > 0)
                    mps[site-1].multiply_from_right(t);
                
                load(right_[site], right_stores_[site]);
                
                MPSTensor<Matrix, SymmGroup> bkp = mps[site];
                right_[site] = contraction::overlap_mpo_right_step(mps[site], bkp,
                                                                   right_[site+1], mpo[site]);
                
                store(left_[site], left_stores_[site]);
                store(right_[site+1], right_stores_[site+1]);
            }
        }
        
        return energies;
    }
    
private:
    void init_left_right(MPO<Matrix, SymmGroup> const & mpo)
    {
        static Timer timer("init_left_right wait");
        
        std::size_t L = mps.length();
        
        left_.resize(mpo.length()+1);
        right_.resize(mpo.length()+1);
        
        right_stores_.resize(L+1, storage_master.child());
        left_stores_.resize(L+1, storage_master.child());
        
        {
            Boundary<Matrix, SymmGroup> left = mps.left_boundary();
            left_[0] = left;
            reset(left_stores_[0]);
            store(left_[0], left_stores_[0]);
            
            for (int i = 0; i < L; ++i) {
                MPSTensor<Matrix, SymmGroup> bkp = mps[i];
                left = contraction::overlap_mpo_left_step(mps[i], bkp, left, mpo[i]);
                left_[i+1] = left;
                reset(left_stores_[i+1]);
                store(left_[i+1], left_stores_[i+1]);
            }
        }
        
        {
            Boundary<Matrix, SymmGroup> right = mps.right_boundary();
            right_[L] = right;
            reset(right_stores_[L]);
            store(right_[L], right_stores_[L]);
            
            for (int i = L-1; i >= 0; --i) {
                MPSTensor<Matrix, SymmGroup> bkp = mps[i];
                right = contraction::overlap_mpo_right_step(mps[i], bkp, right, mpo[i]);
                right_[i] = right;
                reset(right_stores_[i]);
                store(right_[i], right_stores_[i]);
            }
        }
    }
    
    MPS<Matrix, SymmGroup> & mps;
    BaseParameters & parms;
    std::vector<Boundary<Matrix, SymmGroup> > left_, right_;
    std::vector<typename StorageMaster::Storage> left_stores_, right_stores_;
    StorageMaster & storage_master;
};

#endif

