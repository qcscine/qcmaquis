/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *               2017-2017 by Alberto Baiardi <alberto.baiardi@sns.it>
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

#ifndef IETL_JD_SOLVER_H
#define IETL_JD_SOLVER_H

#include "ietl_lanczos_solver.h"

#include "ietl/jd.h"
#include "dmrg/optimize/jacobi_standard.h"
#include "dmrg/optimize/partial_overlap.h"

//
// Standard Jacobi-Davidson diagonalization
// ----------------------------------------

template<class Matrix, class SymmGroup>
std::pair<double, class MPSTensor<Matrix, SymmGroup> >
solve_ietl_jcd(SiteProblem<Matrix, SymmGroup> & sp,
               MPSTensor<Matrix, SymmGroup> const & initial,
               BaseParameters & params,
               std::vector< class MPSTensor<Matrix, SymmGroup> > ortho_vecs = std::vector< class MPSTensor<Matrix, SymmGroup> >())
{
    // Standard initialization (as in the Davidson case)
    if (initial.num_elements() <= ortho_vecs.size())
        ortho_vecs.resize(initial.num_elements()-1);
    // Gram-Schmidt the ortho_vecs
    for (int n = 1; n < ortho_vecs.size(); ++n)
        for (int n0 = 0; n0 < n; ++n0)
            ortho_vecs[n] -= ietl::dot(ortho_vecs[n0], ortho_vecs[n])/ietl::dot(ortho_vecs[n0],ortho_vecs[n0])*ortho_vecs[n0];
    // Definition of the method to solve the non-linear equation required in the Jacobi-Davidson algorithm
    typedef MPSTensor<Matrix, SymmGroup> Vector;
    SingleSiteVS<Matrix, SymmGroup> vs(initial, ortho_vecs);
    double tol = params["ietl_diag_tol"];
    ietl::basic_iteration<double> iter(params["ietl_diag_maxiter"], tol, tol);
    ietl::jacobi_davidson_standard<SiteProblem<Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup> , ietl::basic_iteration<double> > jd(sp, vs);
    maquis::cout << "Ortho vecs " << ortho_vecs.size() << std::endl;
    for (int n = 0; n < ortho_vecs.size(); ++n) {
        maquis::cout << "Ortho norm " << n << ": " << ietl::two_norm(ortho_vecs[n]) << std::endl;
        maquis::cout << "Input <MPS|O[" << n << "]> : " << ietl::dot(initial, ortho_vecs[n]) << std::endl;
    }
    std::pair<double, Vector> r0 = jd.calculate_eigenvalue(initial, iter);
    for (int n = 0; n < ortho_vecs.size(); ++n)
        maquis::cout << "Output <MPS|O[" << n << "]> : " << ietl::dot(r0.second, ortho_vecs[n]) << std::endl;
    maquis::cout << "JCD used " << iter.iterations() << " iterations." << std::endl;
    return r0;
}


#endif
