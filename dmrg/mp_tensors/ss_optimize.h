#ifndef SS_OPTIMIZE_H
#define SS_OPTIMIZE_H

#include <boost/random.hpp>
#include <sys/time.h>

#include "utils/zout.hpp"
#include "ietl_lanczos_solver.h"
#include "arpackpp_solver.h"
#include "ietl_jacobi_davidson.h"

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
zout << "Time elapsed in " << name << ": " << then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec) << endl;

double log_interpolate(double y0, double y1, int N, int i)
{
    double x = log(y1/y0)/(N-1);
    return y0*exp(x*i);
}

template<class Matrix, class SymmGroup>
class ss_optimize
{
public:
    ss_optimize(MPS<Matrix, SymmGroup> & mps_,
                BaseParameters & parms_)
    : mps(mps_)
    , parms(parms_)
    { }
    
    std::vector<double> sweep(MPO<Matrix, SymmGroup> const & mpo,
                              int sweep)
    {
        mps.normalize_right();
        mps.canonize(0);
        
        left_ = left_mpo_overlaps(mps, mpo);
        right_ = right_mpo_overlaps(mps, mpo);
        
        std::vector<double> energies;
        
        std::size_t L = mps.length();
        
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
                res = solve_ietl_lanczos(sp, mps[site], parms);
                END_TIMING("IETL")
            } else if (parms.get<std::string>("eigensolver") == std::string("ARPACK")) {
                BEGIN_TIMING("ARPACK")
                res = solve_arpackpp(sp, mps[site], parms);
                END_TIMING("ARPACK")
            } else {
                throw std::runtime_error("I don't know this eigensolver.");
            }
			
			/*
			else if (parms.get<std::string>("eigensolver") == std::string("IETL_JCD")) {
                BEGIN_TIMING("JCD")
                res = solve_ietl_jcd(sp, mps[site], parms);
                END_TIMING("JCD")
            } else {
                throw std::runtime_error("I don't know this eigensolver.");
            }
		    */
            mps[site] = res.second;
            
            zout << "Energy " << /* sweep+(lr == -1 ? 0.5 : 0) */ lr << " " << res.first << endl;
            energies.push_back(res.first);
            
            double alpha = log_interpolate(parms.get<double>("alpha_initial"), parms.get<double>("alpha_final"), parms.get<int>("nsweeps"), sweep);
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
                
                MPSTensor<Matrix, SymmGroup> bkp = mps[site];
                left_[site+1] = contraction::overlap_mpo_left_step(mps[site], bkp,
                                                                   left_[site], mpo[site]);
            } else if (lr == -1) {
                if (site > 1) {
                    zout << "Growing, alpha = " << alpha << endl;
                    mps.grow_r2l_sweep(mpo[site], left_[site], right_[site+1],
                                       site, alpha, cutoff, Mmax);
                }
                
                block_matrix<Matrix, SymmGroup> t = mps[site].normalize_right(SVD);
                if (site > 0)
                    mps[site-1].multiply_from_right(t);
                
                MPSTensor<Matrix, SymmGroup> bkp = mps[site];
                right_[site] = contraction::overlap_mpo_right_step(mps[site], bkp,
                                                                   right_[site+1], mpo[site]);
            }
        }
        
        return energies;
    }
    
private:
    MPS<Matrix, SymmGroup> & mps;
    BaseParameters & parms;
    std::vector<Boundary<Matrix, SymmGroup> > left_, right_;
};

#endif

