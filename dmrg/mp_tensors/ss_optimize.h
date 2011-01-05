#ifndef SS_OPTIMIZE_H
#define SS_OPTIMIZE_H

#include <boost/random.hpp>
#include <sys/time.h>

#include "ietl_lanczos_solver.h"
#include "arpackpp_solver.h"

template<class Matrix, class SymmGroup>
struct SiteProblem
{
    MPSTensor<Matrix, SymmGroup> ket_tensor;
    Boundary<Matrix, SymmGroup> left, right;
    MPOTensor<Matrix, SymmGroup> mpo;
};

#define BEGIN_TIMING(name) \
gettimeofday(&now, NULL);
#define END_TIMING(name) \
gettimeofday(&then, NULL); \
cout << "Time elapsed in " << name << ": " << then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec) << endl;

template<class Matrix, class SymmGroup>
void ss_optimize(MPS<Matrix, SymmGroup> & mps,
                 MPO<Matrix, SymmGroup> const & mpo,
                 int nsweeps, double cutoff, std::size_t Mmax)
{   
    std::size_t L = mps.length();
    
    mps.normalize_right();
    mps.canonize(0);
    std::vector<Boundary<Matrix, SymmGroup> >
    left_ = left_mpo_overlaps(mps, mpo),
    right_ = right_mpo_overlaps(mps, mpo);
    
    for (int sweep = 0; sweep < nsweeps; ++sweep) {
        cout << mps.description() << endl;
        cutoff *= 0.2;
        for (int _site = 0; _site < 2*L; ++_site) {
            int site, lr;
            if (_site < L) {
                site = _site;
                lr = 1;
            } else {
                site = 2*L-_site-1;
                lr = -1;
            }
            
            cout << "Sweep " << sweep << ", optimizing site " << site << endl;
            
            mps[site].make_left_paired();
            
            SiteProblem<Matrix, SymmGroup> sp;
            sp.ket_tensor = mps[site];
            sp.mpo = mpo[site];
            sp.left = left_[site];
            sp.right = right_[site+1];
            
            timeval now, then;
            
//            BEGIN_TIMING("IETL")
//            std::pair<double, MPSTensor<Matrix, SymmGroup> > res = solve_ietl_lanczos(sp, mps[site]);
//            END_TIMING("IETL")
            
            BEGIN_TIMING("ARPACK")
            std::pair<double, MPSTensor<Matrix, SymmGroup> > res2 = solve_arpackpp(sp, mps[site]);
            END_TIMING("ARPACK")
            
            mps[site] = res2.second;
            
            cout << "Energy: " << res2.first << endl;
//            cout << "Energy comparison: " << res.first << " " << res2.first << endl;
            
            if (lr == +1) {
                if (site < L-1) {
                    double alpha = 1e-3 * powf(0.5, sweep);
                    cout << "Growing, alpha = " << alpha << endl;
                    mps.grow_l2r_sweep(mpo[site], left_[site], right_[site+1],
                                       site, alpha, cutoff, Mmax);
                }
                
                block_matrix<Matrix, SymmGroup> t = mps[site].normalize_left(SVD);
                if (site < L-1)
                    mps[site+1].multiply_from_left(t);
                
                MPSTensor<Matrix, SymmGroup> bkp = mps[site];
                left_[site+1] = contraction::overlap_mpo_left_step(mps[site], bkp,
                                                                   left_[site], mpo[site]);
            } else if (lr == -1) {
                block_matrix<Matrix, SymmGroup> t = mps[site].normalize_right(SVD);
                if (site > 0)
                    mps[site-1].multiply_from_right(t);
                
                MPSTensor<Matrix, SymmGroup> bkp = mps[site];
                right_[site] = contraction::overlap_mpo_right_step(mps[site], bkp,
                                                                   right_[site+1], mpo[site]);
            }
        }
    }
}

#endif

