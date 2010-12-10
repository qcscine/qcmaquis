#include <cmath>
#include <iterator>
#include <iostream>

using std::cout;
using std::cerr;
using std::endl;

#include "general_matrix.hpp"
#include "matrix_interface.hpp"
#include "resizable_matrix_interface.hpp"
#include "general_matrix_algorithms.h"
#include "matrix_algorithms.hpp"
typedef blas::general_matrix<double> Matrix;

#include "indexing.h"
#include "mps.h"
#include "mpo.h"
#include "contractions.h"
#include "mps_mpo_ops.h"

#include "special_mpos.h"

typedef NullGroup grp;

typedef std::vector<MPOTensor<Matrix, grp> > mpo_t;
typedef Boundary<Matrix, grp> boundary_t;

#include <boost/random.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/io.hpp>

typedef boost::numeric::ublas::vector<double> Vector;

struct SiteProblem
{
    MPSTensor<Matrix, grp> ket_tensor;
    boundary_t left, right;
    MPOTensor<Matrix, grp> mpo;
};

namespace ietl
{
    void mult(SiteProblem const & H, Vector const & x, Vector & y)
    {   
        MPSTensor<Matrix, grp> t1(H.ket_tensor), t2;
        t1.make_left_paired();
        
        Vector::const_iterator cit = x.begin();
        for (std::size_t k = 0; k < t1.data().n_blocks(); ++k) {
            std::copy(cit, cit+num_rows(t1.data()[k])*num_columns(t1.data()[k]), elements(t1.data()[k]).first);
            cit += num_rows(t1.data()[k])*num_columns(t1.data()[k]);
        }
        
        t2 = contraction::site_hamil(t1, H.left, H.right, H.mpo);
        
        Vector::iterator it = y.begin();
        for (std::size_t k = 0; k < t2.data().n_blocks(); ++k)
            it = std::copy(elements(t2.data()[k]).first, elements(t2.data()[k]).second, it);
        
        //        cout << x << " " << y << endl;
        //        cout << "IP: " << inner_prod(x-y, x-y) << endl;
    }
}

void v2m(Vector const & x, MPSTensor<Matrix, grp> & t1)
{
    //    cout << "t1 before: " << t1 << endl;
    
    Vector::const_iterator cit = x.begin();
    for (std::size_t k = 0; k < t1.data().n_blocks(); ++k) {
        std::copy(cit, cit+num_rows(t1.data()[k])*num_columns(t1.data()[k]), elements(t1.data()[k]).first);
        cit += num_rows(t1.data()[k])*num_columns(t1.data()[k]);
    }
    
    //    cout << "t1 after: " << t1 << endl;
    //    exit(0);
}

#include <ietl/interface/ublas.h>
#include <ietl/lanczos.h>
#include <ietl/vectorspace.h>

void optimize(MPS<Matrix, grp> & mps,
              MPO<Matrix, grp> const & mpo)
{
    std::size_t L = mps.length();
    
    mps.normalize_right();
    mps.canonize(0);
    std::vector<boundary_t>
    left_ = left_mpo_overlaps(mps, mpo),
    right_ = right_mpo_overlaps(mps, mpo);
    
    for (int sweep = 0; sweep < 1; ++sweep) {
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
//            normalize_upto(site);
//            mps[site].multiply_by_scalar(1/norm());
//            left_ = left_mpo_overlaps(mpo);
//            right_ = right_mpo_overlaps(mpo);
//
//            cout << "Normerr: " << norm() - 1 << endl;
//            assert( fabs(norm() - 1) < 1e-6 );
            
            unsigned int N = 0;
            for (std::size_t k = 0; k < mps[site].data().n_blocks(); ++k)
                N += num_rows(mps[site].data()[k]) * num_columns(mps[site].data()[k]);
            ietl::vectorspace<Vector> vs(N);
            boost::lagged_fibonacci607 gen;
            
            // This should not be necessary here.
            // If you find out why it is, let me know.
            mps[site].make_left_paired();
            SiteProblem sp;
            sp.ket_tensor = mps[site];
            sp.mpo = mpo[site];
            sp.left = left_[site];
            sp.right = right_[site+1];
            
            typedef ietl::vectorspace<Vector> Vecspace;
            typedef boost::lagged_fibonacci607 Gen;  
            
            Vecspace vec(N);
            Gen mygen;
            ietl::lanczos<SiteProblem,Vecspace> lanczos(sp,vec);
            
            double rel_tol = 500*std::numeric_limits<double>::epsilon();
            double abs_tol = std::pow(std::numeric_limits<double>::epsilon(),2./3);
            ietl::lanczos_iteration_nlowest<double> 
            iter(100, 1, rel_tol, abs_tol);
            
            std::vector<double> eigen, err;
            std::vector<int> multiplicity;  
            
            try{
                lanczos.calculate_eigenvalues(iter,mygen);
                eigen = lanczos.eigenvalues();
                err = lanczos.errors();
                multiplicity = lanczos.multiplicities();
                std::cout << "number of iterations: " << iter.iterations() << "\n";
            }
            catch (std::runtime_error& e) {
                cout << "Error in eigenvalue calculation: " << endl;
                cout << e.what() << endl;
            }
            
            cout << "Energy: " << eigen[0] << endl;
            
            std::vector<double>::iterator start = eigen.begin();  
            std::vector<double>::iterator end = eigen.begin()+1;
            std::vector<Vector> eigenvectors; // for storing the eigen vectors. 
            ietl::Info<double> info; // (m1, m2, ma, eigenvalue, residualm, status).
            
            try {
                lanczos.eigenvectors(start,end,std::back_inserter(eigenvectors),info,mygen); 
            }
            catch (std::runtime_error& e) {
                cout << "Error in eigenvector calculation: " << endl;
                cout << e.what() << endl;
            }  
            
            v2m(eigenvectors[0], mps[site]);
            
            if (lr == +1) {
                block_matrix<Matrix, grp> t = mps[site].normalize_left(SVD);
                if (site < L-1)
                    mps[site+1].multiply_from_left(t);
                
                
                MPSTensor<Matrix, grp> bkp = mps[site];
                left_[site+1] = contraction::overlap_mpo_left_step(mps[site], bkp,
                                                                   left_[site], mpo[site]);
            } else if (lr == -1) {
                block_matrix<Matrix, grp> t = mps[site].normalize_right(SVD);
                if (site > 0)
                    mps[site-1].multiply_from_right(t);
                
                MPSTensor<Matrix, grp> bkp = mps[site];
                right_[site] = contraction::overlap_mpo_right_step(mps[site], bkp,
                                                                   right_[site+1], mpo[site]);
            }
        }
    }
}

int main()
{
    Index<grp> phys; phys.insert(std::make_pair(NullGroup::Plus, 2));
    
    int L = 16, M = 20;
    MPS<Matrix, grp> mps(L, M, phys);
    
    MPOTensor<Matrix, grp> id_mpo = identity_mpo<Matrix>(mps[0].site_dim());
    MPOTensor<Matrix, grp> sz_mpo = s12_sz_mpo<Matrix>(mps[0].site_dim());
    
    MPO<Matrix, grp> szsz = s12_ising<Matrix>(L, -1, 1);
    optimize(mps, szsz);
}
