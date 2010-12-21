#ifndef SS_OPTIMIZE_H
#define SS_OPTIMIZE_H

#include <boost/random.hpp>

template<class Matrix, class SymmGroup>
struct SiteProblem
{
    MPSTensor<Matrix, SymmGroup> ket_tensor;
    Boundary<Matrix, SymmGroup> left, right;
    MPOTensor<Matrix, SymmGroup> mpo;
};

template<class Matrix, class SymmGroup>
class SingleSiteVS
{
public:
    SingleSiteVS(MPSTensor<Matrix, SymmGroup> const & m)
    : instance(m)
    {
        for (std::size_t k = 0; k < m.data().n_blocks(); ++k)
            N += num_rows(m.data()[k]) * num_columns(m.data()[k]);
    }
  
    friend MPSTensor<Matrix, SymmGroup> new_vector(SingleSiteVS const & vs)
    {
        return vs.instance;
    }
    
    friend std::size_t vec_dimension(SingleSiteVS const & vs)
    {
        return vs.N;
    }
    
    void project(MPSTensor<Matrix, SymmGroup> & x) const { }
    
private:
    MPSTensor<Matrix, SymmGroup> instance;
    std::size_t N;
};

#include <ietl/vectorspace.h>

namespace ietl
{
    template<class Matrix, class SymmGroup>
    void mult(SiteProblem<Matrix, SymmGroup> const & H,
              MPSTensor<Matrix, SymmGroup> const & x,
              MPSTensor<Matrix, SymmGroup> & y)
    {  
        y = contraction::site_hamil(x, H.left, H.right, H.mpo);
        x.make_left_paired();
    }
    
    template<class Matrix, class SymmGroup>
    struct vectorspace_traits<SingleSiteVS<Matrix, SymmGroup> >
    {
        typedef MPSTensor<Matrix, SymmGroup> vector_type;
        typedef typename MPSTensor<Matrix, SymmGroup>::scalar_type scalar_type;
        typedef typename MPSTensor<Matrix, SymmGroup>::real_type magnitude_type;
        typedef std::size_t size_type;
    };
}

#include <ietl/lanczos.h>

template<class Matrix, class SymmGroup>
void ss_optimize(MPS<Matrix, SymmGroup> & mps,
                 MPO<Matrix, SymmGroup> const & mpo,
                 int nsweeps)
{
    typedef MPSTensor<Matrix, SymmGroup> Vector;
    
    std::size_t L = mps.length();
    
    mps.normalize_right();
    mps.canonize(0);
    std::vector<Boundary<Matrix, SymmGroup> >
    left_ = left_mpo_overlaps(mps, mpo),
    right_ = right_mpo_overlaps(mps, mpo);
    
    for (int sweep = 0; sweep < nsweeps; ++sweep) {
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
            SingleSiteVS<Matrix, SymmGroup> vs(mps[site]);
            
            SiteProblem<Matrix, SymmGroup> sp;
            sp.ket_tensor = mps[site];
            sp.mpo = mpo[site];
            sp.left = left_[site];
            sp.right = right_[site+1];
            
            typedef ietl::vectorspace<Vector> Vecspace;
            typedef boost::lagged_fibonacci607 Gen;  
            
            Gen mygen;
            
            ietl::lanczos<SiteProblem<Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup> > lanczos(sp, vs);
            
            double rel_tol = sqrt(std::numeric_limits<double>::epsilon());
            double abs_tol = rel_tol; // std::pow(std::numeric_limits<double>::epsilon(),2./3);
            ietl::lanczos_iteration_nlowest<double> 
            iter(100, 1, rel_tol, abs_tol);
            
            std::vector<double> eigen, err;
            std::vector<int> multiplicity;  
            
            try{
                if (sweep == 0)
                    lanczos.calculate_eigenvalues(iter, mygen);
                else
                    lanczos.calculate_eigenvalues(iter, mps[site]);
                eigen = lanczos.eigenvalues();
                err = lanczos.errors();
                multiplicity = lanczos.multiplicities();
                std::cout << "number of iterations: " << iter.iterations() << "\n";
            }
            catch (std::runtime_error& e) {
                cout << "Error in eigenvalue calculation: " << endl;
                cout << e.what() << endl;
                exit(1);
            }
            
            cout << "Energy: " << eigen[0] << endl;
            
            std::vector<double>::iterator start = eigen.begin();  
            std::vector<double>::iterator end = eigen.begin()+1;
            std::vector<Vector> eigenvectors; // for storing the eigen vectors. 
            ietl::Info<double> info; // (m1, m2, ma, eigenvalue, residualm, status).
            
            try {
                if (sweep == 0)
                    lanczos.eigenvectors(start,end,std::back_inserter(eigenvectors), info, mygen); 
                else
                    lanczos.eigenvectors(start,end,std::back_inserter(eigenvectors), info, mps[site]); 
            }
            catch (std::runtime_error& e) {
                cout << "Error in eigenvector calculation: " << endl;
                cout << e.what() << endl;
                exit(1);
            }  
            
            mps[site] = eigenvectors[0];
            
            if (lr == +1) {
                if (site < L-1) {
                    cout << "Growing..." << endl;
                    mps.grow_l2r_sweep(mpo[site], left_[site], right_[site+1],
                                       site, 1e-4, 1e-10);
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

