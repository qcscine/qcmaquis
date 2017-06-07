/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2015 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *               2011-2015 by Sebastian Keller <sebkelle@phys.ethz.ch>
 * 
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 * 
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#ifndef IETL_DAVIDSON_SOLVER_H
#define IETL_DAVIDSON_SOLVER_H

#include "dmrg/utils/BaseParameters.h"
#include "ietl_lanczos_solver.h"

#include "davidson.h"

namespace davidson_detail {
    template<class Matrix, class SymmGroup> class MultDiagonal
    {
        typedef MPSTensor<Matrix, SymmGroup> vector_type;
        typedef typename Matrix::value_type value_type;
        typedef typename MPSTensor<Matrix, SymmGroup>::scalar_type scalar_type ;
    public:
        // The construtor sets the Hdiag attribute. Overloaded depending on the
        // presence of the shift
        MultDiagonal(SiteProblem<Matrix, SymmGroup> const& H, vector_type const& x)
        : have_omega(false) {
            Hdiag = contraction::diagonal_hamiltonian(H.left, H.right, H.mpo, x);
        }
        MultDiagonal(SiteProblem<Matrix, SymmGroup> const& H, vector_type const& x, double const& omega)
                : have_omega(true) , omega_(omega){
            Hdiag = contraction::diagonal_hamiltonian(H.left, H.right, H.mpo, x);
        }

        void precondition(vector_type& r, vector_type& V, value_type theta)
        {
            // Project onto the space orthogonal to V
            value_type  x1   = ietl::dot(V,r) ;
            vector_type Vcpy = r - V*x1 ;
            // Precondition
            mult_diag(theta, Vcpy);
            // Reproject again
            value_type x2    = ietl::dot(V,Vcpy);
            r = Vcpy - x2*V ;
        }

    private:
        // Preconditioner in Davidson diagonalization
        // Takes in input a (most likely) float number (theta) and
        // preconditions an MPSTensor that is given in input (x)
        void mult_diag(value_type theta, vector_type& x)
        {
            value_type denom ;
            block_matrix<Matrix, SymmGroup> & data = x.data();
            assert(shape_equal(data, Hdiag));
            for (size_t b = 0; b < data.n_blocks(); ++b)
            {
                for (size_t i = 0; i < num_rows(data[b]); ++i)
                    for (size_t j = 0; j < num_cols(data[b]); ++j) {
                        if (have_omega)
                            denom = (omega_ - Hdiag[b](i, j)) - theta;
                        else
                            denom = Hdiag[b](i, j) - theta;
                        if (std::abs(denom))
                            data[b](i, j) /= denom;
                    }
            }
        }
        block_matrix<Matrix, SymmGroup> Hdiag;
        bool have_omega ;
        double omega_  ;
    };
} // namespace davidson detail

// +--------------------+
//  SOLVE_IETL_DAVIDSON
// +--------------------+
//
// Simple Davidson
// ---------------
template<class Matrix, class SymmGroup>
std::pair< double , class MPSTensor<Matrix,SymmGroup> >
solve_ietl_davidson(SiteProblem<Matrix, SymmGroup> & sp,
                    MPSTensor<Matrix, SymmGroup> const & initial,
                    BaseParameters & params,
                    std::vector< class MPSTensor<Matrix, SymmGroup> > ortho_vecs = std::vector< class MPSTensor<Matrix, SymmGroup> >()) {
    // Check if the number of MPSTensors is higher than the one of the orthogonal vectors
    // and performs the GS orthogonalization
    if (initial.num_elements() <= ortho_vecs.size())
        ortho_vecs.resize(initial.num_elements()-1);
    for (int n = 1; n < ortho_vecs.size(); ++n)
        for (int n0 = 0; n0 < n; ++n0)
            ortho_vecs[n] -= ietl::dot(ortho_vecs[n0], ortho_vecs[n])/ietl::dot(ortho_vecs[n0],ortho_vecs[n0])*ortho_vecs[n0];
    // Initialization
    typedef MPSTensor<Matrix, SymmGroup> Vector ;
    SingleSiteVS<Matrix, SymmGroup> vs(initial, ortho_vecs);
    // Create the Davidson object
    ietl::davidson<SiteProblem<Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup> > davidson(sp, vs);
    davidson_detail::MultDiagonal<Matrix, SymmGroup> mdiag(sp, initial);
    // Sets the iterator object
    double tol = params["ietl_jcd_tol"];
    ietl::basic_iteration<double> iter(params["ietl_davidson_maxiter"], tol, tol);
    //contraction::ContractionGrid<Matrix, SymmGroup>::iterate_reduction_layout(0, params["ietl_davidson_maxiter"]);
    // Check orthogonality
    for (int n = 0; n < ortho_vecs.size(); ++n) {
        maquis::cout << "Input <MPS|O[" << n << "]> : " << ietl::dot(initial, ortho_vecs[n]) << std::endl;
    }
    // Compute eigenvalue
    std::pair<double, Vector> r0 = davidson.calculate_eigenvalue(initial, mdiag, iter);
    // Check again orthogonality
    for (int n = 0; n < ortho_vecs.size(); ++n)
        maquis::cout << "Output <MPS|O[" << n << "]> : " << ietl::dot(r0.second, ortho_vecs[n]) << std::endl;
    maquis::cout << "Davidson used " << iter.iterations() << " iterations." << std::endl;
    // Returns in output a vector and the corresponding eigenvector (the energy)
    return r0;
}
//
// Modified Davidson
// -----------------
template<class Matrix, class SymmGroup>
std::pair< double , class MPSTensor<Matrix,SymmGroup> >
solve_ietl_davidson_modified(SiteProblem<Matrix, SymmGroup> & sp,
        MPSTensor<Matrix, SymmGroup> const & initial,
        BaseParameters & params,
        double omega,
        std::vector< class MPSTensor<Matrix, SymmGroup> > ortho_vecs = std::vector< class MPSTensor<Matrix, SymmGroup> >()) {
    // Check if the number of MPSTensors is higher than the one of the orthogonal vectors
    // and performs the GS orthogonalization
    if (initial.num_elements() <= ortho_vecs.size())
    ortho_vecs.resize(initial.num_elements()-1);
    for (int n = 1; n < ortho_vecs.size(); ++n)
        for (int n0 = 0; n0 < n; ++n0)
            ortho_vecs[n] -= ietl::dot(ortho_vecs[n0], ortho_vecs[n])/ietl::dot(ortho_vecs[n0],ortho_vecs[n0])*ortho_vecs[n0];
    // Initialization
    typedef MPSTensor<Matrix, SymmGroup> Vector ;
    SingleSiteVS<Matrix, SymmGroup> vs(initial, ortho_vecs);
    // Create the Davidson object
    ietl::davidson_modified<SiteProblem<Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup> > davidson_modified(sp, vs, omega);
    davidson_detail::MultDiagonal<Matrix, SymmGroup> mdiag(sp, initial, omega);
    // Sets the iterator object
    double tol = params["ietl_jcd_tol"];
    ietl::basic_iteration<double> iter(params["ietl_davidson_maxiter"], tol, tol);
    contraction::ContractionGrid<Matrix, SymmGroup>::iterate_reduction_layout(0, params["ietl_davidson_maxiter"]);
    // Check orthogonality
    for (int n = 0; n < ortho_vecs.size(); ++n) {
        maquis::cout << "Input <MPS|O[" << n << "]> : " << ietl::dot(initial, ortho_vecs[n]) << std::endl;
    }
    // Compute eigenvalue
    std::pair<double, Vector> r0 = davidson_modified.calculate_eigenvalue(initial, mdiag, iter);
    // Check again orthogonality
    for (int n = 0; n < ortho_vecs.size(); ++n)
        maquis::cout << "Output <MPS|O[" << n << "]> : " << ietl::dot(r0.second, ortho_vecs[n]) << std::endl;
    maquis::cout << "Modified Davidson used " << iter.iterations() << " iterations." << std::endl;
    // Returns in output a vector and the corresponding eigenvector (the energy)
    return r0;
}



#endif
