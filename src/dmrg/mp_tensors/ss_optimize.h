#ifndef SS_OPTIMIZE_H
#define SS_OPTIMIZE_H

#include <boost/random.hpp>
#include <sys/time.h>

#include "ietl_lanczos_solver.h"
#include "arpackpp_solver.h"

#include "utils/DmrgParameters.h"

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
std::vector<double> ss_optimize(MPS<Matrix, SymmGroup> & mps,
                                MPO<Matrix, SymmGroup> const & mpo,
                                DmrgParameters & parms)
{   
    std::vector<double> energies;
    
    std::size_t L = mps.length();
    
    mps.normalize_right();
    mps.canonize(0);
    std::vector<Boundary<Matrix, SymmGroup> >
    left_ = left_mpo_overlaps(mps, mpo),
    right_ = right_mpo_overlaps(mps, mpo);
    
    for (int sweep = 0; sweep < parms.get<int>("nsweeps"); ++sweep) {
        cout << mps.description() << endl;
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
            
            std::pair<double, MPSTensor<Matrix, SymmGroup> > res;
            
            if (parms.get<std::string>("eigensolver") == std::string("IETL")) {
                BEGIN_TIMING("IETL")
                res = solve_ietl_lanczos(sp, mps[site]);
                END_TIMING("IETL")
            } else if (parms.get<std::string>("eigensolver") == std::string("ARPACK")) {
                BEGIN_TIMING("ARPACK")
                res = solve_arpackpp(sp, mps[site]);
                END_TIMING("ARPACK")
            } else {
                throw std::runtime_error("I don't know this eigensolver.");
            }

            
            mps[site] = res.second;
            
            cout << "Energy: " << res.first << endl;
            energies.push_back(res.first);
            
            if (lr == +1) {
                if (site < L-1) {
                    double alpha = parms.get<double>("alpha_initial") * powf(parms.get<double>("alpha_sweep_factor"), sweep);
                    
                    double t1 = parms.get<double>("truncation_initial"), t2 = parms.get<double>("truncation_final");
                    double cutoff = t1 + static_cast<double>(sweep)/(parms.get<int>("nsweeps")-1)*(t2-t1);
                    
                    std::size_t Mmax = parms.get<std::size_t>("max_bond_dimension");
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
    
    return energies;
}

#endif

