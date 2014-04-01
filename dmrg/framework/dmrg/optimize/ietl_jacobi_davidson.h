/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
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

#include "dmrg/utils/BaseParameters.h"

#include "ietl_lanczos_solver.h"

#include "ietl/jacobi.h"
#include "ietl/jd.h"

template<class Matrix, class SymmGroup>
std::pair<double, MPSTensor<Matrix, SymmGroup> >
solve_ietl_jcd(SiteProblem<Matrix, SymmGroup> & sp,
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
    
    ietl::jacobi_davidson<SiteProblem<Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup> >
    jd(sp, vs, ietl::Smallest);
    
    double tol = params["ietl_jcd_tol"];
    ietl::basic_iteration<double> iter(params["ietl_jcd_maxiter"], tol, tol);
    
//    maquis::cout << "Ortho vecs " << ortho_vecs.size() << std::endl;
    for (int n = 0; n < ortho_vecs.size(); ++n) {
//        maquis::cout << "Ortho norm " << n << ": " << ietl::two_norm(ortho_vecs[n]) << std::endl;
        maquis::cout << "Input <MPS|O[" << n << "]> : " << ietl::dot(initial, ortho_vecs[n]) << std::endl;
    }
    
    std::pair<double, Vector> r0 = jd.calculate_eigenvalue(initial, jcd_gmres, iter);

    for (int n = 0; n < ortho_vecs.size(); ++n)
        maquis::cout << "Output <MPS|O[" << n << "]> : " << ietl::dot(r0.second, ortho_vecs[n]) << std::endl;
    
    maquis::cout << "JCD used " << iter.iterations() << " iterations." << std::endl;
    
    return r0;
}

/*template<class Matrix, class SymmGroup>
std::pair<double, MPSTensor<Matrix, SymmGroup> >
solve_ietl_new_jd(SiteProblem<Matrix, SymmGroup> & sp,
                  MPSTensor<Matrix, SymmGroup> const & initial,
                  BaseParameters & params)
{
    typedef MPSTensor<Matrix, SymmGroup> Vector;
    typedef SiteProblem<Matrix, SymmGroup> Operator;
    typedef SingleSiteVS<Matrix, SymmGroup> Vecspace;
    Vecspace vs(initial);

    int n_evals = 1;   // number of eigenpairs to be calculated
    int max_iter = params["ietl_jcd_maxiter"]; // maximal number of iterations

    int m_max = 40;
    int m_min = 20;

    // tolerance
    double rel_tol = params["ietl_jcd_tol"];
    double abs_tol = rel_tol;

    // maximal iterations for the correction equation solver
    unsigned max_cor_iter = params["ietl_jcd_gmres"];

    ietl::jd_iteration<double> iter(max_iter, m_min, m_max, rel_tol, abs_tol);
    ietl::jd<Operator, Vecspace> jd(sp, vs);

    // the correction equation solver must be an function object
    ietl::ietl_gmres gmres(max_cor_iter);

    jd.eigensystem(iter, initial, n_evals, gmres);
    
    std::vector<double> evals = jd.eigenvalues();
    
    return std::make_pair(evals[0], jd.eigenvector(0));
}*/

#endif
