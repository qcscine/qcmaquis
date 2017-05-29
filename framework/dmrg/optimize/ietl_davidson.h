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
        MultDiagonal(SiteProblem<Matrix, SymmGroup> const& H, vector_type const& x)
        {
            Hdiag = contraction::diagonal_hamiltonian(H.left, H.right, H.mpo, x);
            //vector_type y ;
            //scalar_type overlap ;
            //y = contraction::Engine<Matrix, Matrix, SymmGroup>::site_hamil2(x, H.left, H.right, H.mpo) ;
            //overlap = y.scalar_overlap(x) ;
            //std::cout << "Overlap" << overlap << std::endl ;
            //maquis::cerr << "Exception thrown!" << std::endl ;
            //exit(1) ;
        }

        void precondition(vector_type& r, vector_type& V, value_type theta)
        {
            vector_type Vcpy = V;
            mult_diag(theta, Vcpy);
            value_type a = Vcpy.scalar_overlap(r);
            value_type b = V.scalar_overlap(Vcpy);
            r -= a/b * V;
            mult_diag(theta, r);
        }

    private:
        // Preconditioner in Davidson diagonalization
        // Takes in input a (most likely) float number (theta) and
        // preconditions a vector that is given in input (x)
        void mult_diag(value_type theta, vector_type& x)
        {
            value_type shift ;
            shift = 30. ;
            block_matrix<Matrix, SymmGroup> & data = x.data();
            assert(shape_equal(data, Hdiag));
            for (size_t b = 0; b < data.n_blocks(); ++b)
            {
                for (size_t i = 0; i < num_rows(data[b]); ++i)
                    for (size_t j = 0; j < num_cols(data[b]); ++j)
                        if (std::abs(theta - (Hdiag[b](i,j)-shift) ))
                            data[b](i,j) /= (theta - (Hdiag[b](i,j)-shift) );
            }
        }
        block_matrix<Matrix, SymmGroup> Hdiag;
    };

} // namespace davidson detail

// +--------------------+
//  SOLVE_IETL_DAVIDSON
// +--------------------+

template<class Matrix, class SymmGroup>
std::pair<double, MPSTensor<Matrix, SymmGroup> >
solve_ietl_davidson(SiteProblem<Matrix, SymmGroup> & sp,
                    MPSTensor<Matrix, SymmGroup> const & initial,
                    BaseParameters & params,
                    std::vector<MPSTensor<Matrix, SymmGroup> > ortho_vecs = std::vector<MPSTensor<Matrix, SymmGroup> >())
{
    if (initial.num_elements() <= ortho_vecs.size())
        ortho_vecs.resize(initial.num_elements()-1);
    // Gram-Schmidt the ortho_vecs
    for (int n = 1; n < ortho_vecs.size(); ++n)
        for (int n0 = 0; n0 < n; ++n0)
            ortho_vecs[n] -= ietl::dot(ortho_vecs[n0], ortho_vecs[n])/ietl::dot(ortho_vecs[n0],ortho_vecs[n0])*ortho_vecs[n0];
    typedef MPSTensor<Matrix, SymmGroup> Vector;
    SingleSiteVS<Matrix, SymmGroup> vs(initial, ortho_vecs);
    ietl::jcd_gmres_solver<SiteProblem<Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup> >
    jcd_gmres(sp, vs, params["ietl_jcd_gmres"]);
    
    ietl::davidson<SiteProblem<Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup> >
    jd(sp, vs, ietl::Smallest);

    davidson_detail::MultDiagonal<Matrix, SymmGroup> mdiag(sp, initial);
    
    double tol = params["ietl_jcd_tol"];
    ietl::basic_iteration<double> iter(params["ietl_jcd_maxiter"], tol, tol);
    contraction::ContractionGrid<Matrix, SymmGroup>::iterate_reduction_layout(0, params["ietl_jcd_maxiter"]);
    
    for (int n = 0; n < ortho_vecs.size(); ++n) {
        maquis::cout << "Input <MPS|O[" << n << "]> : " << ietl::dot(initial, ortho_vecs[n]) << std::endl;
    }
    
    std::pair<double, Vector> r0 = jd.calculate_eigenvalue(initial, jcd_gmres, mdiag, iter);

    for (int n = 0; n < ortho_vecs.size(); ++n)
        maquis::cout << "Output <MPS|O[" << n << "]> : " << ietl::dot(r0.second, ortho_vecs[n]) << std::endl;
    
    maquis::cout << "Davidson used " << iter.iterations() << " iterations." << std::endl;
    
    return r0;
}

#endif
