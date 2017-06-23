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

#include "ietl/jacobi.h"
#include "ietl/jd.h"
#include "dmrg/optimize/partial_overlap.h"

//
// Modified Jacobi-Davidson diagonalization
// ----------------------------------------

template<class Matrix, class SymmGroup>
std::pair<double, MPSTensor< Matrix, SymmGroup> >
solve_ietl_jcd_modified(SiteProblem<Matrix, SymmGroup> & sp,
                        MPSTensor<Matrix, SymmGroup> const & initial ,
                        BaseParameters & params ,
                        int site,
                        partial_overlap<Matrix, SymmGroup> poverlap,
                        double omega ,
                        std::vector< MPSTensor<Matrix, SymmGroup> > ortho_vecs = std::vector< MPSTensor<Matrix, SymmGroup> >())
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
    ietl::jcd_gmres_modified_solver<SiteProblem<Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup> >
          jcd_modified_gmres(sp, vs, omega, params["ietl_jcd_gmres"]);
    ietl::jacobi_davidson<SiteProblem<Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup> , Matrix, SymmGroup >
          jd(sp, vs, omega, poverlap, site, params["maximum_overlap_nstates"]);
    double tol = params["ietl_jcd_tol"];
    ietl::basic_iteration<double> iter(params["ietl_jcd_maxiter"], tol, tol);
    maquis::cout << "Ortho vecs " << ortho_vecs.size() << std::endl;
    for (int n = 0; n < ortho_vecs.size(); ++n) {
        maquis::cout << "Ortho norm " << n << ": " << ietl::two_norm(ortho_vecs[n]) << std::endl;
        maquis::cout << "Input <MPS|O[" << n << "]> : " << ietl::dot(initial, ortho_vecs[n]) << std::endl;
    }
    std::pair<double, Vector> r0 = jd.calculate_eigenvalue(initial, jcd_modified_gmres, iter);
    for (int n = 0; n < ortho_vecs.size(); ++n)
        maquis::cout << "Output <MPS|O[" << n << "]> : " << ietl::dot(r0.second, ortho_vecs[n]) << std::endl;
    maquis::cout << "MODIFIED JCD used " << iter.iterations() << " iterations." << std::endl;
    return r0;
}

//
// Standard Jacobi-Davidson diagonalization
// ----------------------------------------

template<class Matrix, class SymmGroup>
std::pair<double, class MPSTensor<Matrix, SymmGroup> >
solve_ietl_jcd(SiteProblem<Matrix, SymmGroup> & sp,
               MPSTensor<Matrix, SymmGroup> const & initial,
               BaseParameters & params,
               int site,
               partial_overlap<Matrix, SymmGroup> poverlap,
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
    ietl::jcd_gmres_solver<SiteProblem<Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup> > jcd_gmres(sp, vs, params["ietl_jcd_gmres"]);
    ietl::jacobi_davidson<SiteProblem<Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup> , Matrix , SymmGroup >
            jd(sp, vs, poverlap, site, params["maximum_overlap_nstates"]);
    double tol = params["ietl_jcd_tol"];
    ietl::basic_iteration<double> iter(params["ietl_jcd_maxiter"], tol, tol);
    maquis::cout << "Ortho vecs " << ortho_vecs.size() << std::endl;
    for (int n = 0; n < ortho_vecs.size(); ++n) {
        maquis::cout << "Ortho norm " << n << ": " << ietl::two_norm(ortho_vecs[n]) << std::endl;
        maquis::cout << "Input <MPS|O[" << n << "]> : " << ietl::dot(initial, ortho_vecs[n]) << std::endl;
    }
    std::pair<double, Vector> r0 = jd.calculate_eigenvalue(initial, jcd_gmres, iter);
    for (int n = 0; n < ortho_vecs.size(); ++n)
        maquis::cout << "Output <MPS|O[" << n << "]> : " << ietl::dot(r0.second, ortho_vecs[n]) << std::endl;
    maquis::cout << "JCD used " << iter.iterations() << " iterations." << std::endl;
    return r0;
}

#endif
