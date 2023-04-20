/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef IETL_JD_SOLVER_H
#define IETL_JD_SOLVER_H

#include "dmrg/utils/BaseParameters.h"
#include "ietl_lanczos_solver.h"
#include "ietl/jacobi.h"
#include "ietl/jd.h"

template<class Matrix, class SymmGroup>
std::pair<typename MPSTensor<Matrix, SymmGroup>::value_type, MPSTensor<Matrix, SymmGroup> >
solve_ietl_jcd(SiteProblem<Matrix, SymmGroup> & sp,
               MPSTensor<Matrix, SymmGroup> const & initial,
               BaseParameters & params,
               std::vector<MPSTensor<Matrix, SymmGroup> > ortho_vecs = std::vector<MPSTensor<Matrix, SymmGroup> >(),
               double thresholdForCompleteness=1.0E-10)
{
    // Variables initialization
    using ValueType = typename MPSTensor<Matrix, SymmGroup>::value_type; 
    std::pair<ValueType, MPSTensor<Matrix, SymmGroup>> r0;
    bool skipOptimization=false;
    auto ortho_vecs_local = std::vector< MPSTensor<Matrix, SymmGroup> >();
    if (initial.num_elements() <= ortho_vecs.size())
        ortho_vecs.resize(initial.num_elements()-1);
    // Gram-Schmidt the ortho_vecs and loads the results in the [ortho_vecs_local]
    for (int n = 0; n < ortho_vecs.size(); ++n) {
        for (const auto& iLocal: ortho_vecs_local)
            ortho_vecs[n] -= ietl::dot(iLocal, ortho_vecs[n])*iLocal;
        if (ortho_vecs[n].scalar_norm() > thresholdForCompleteness) {
            ortho_vecs[n] /= ietl::two_norm(ortho_vecs[n]);
            ortho_vecs_local.push_back(ortho_vecs[n]);
        }
        else {
            maquis::cout << "State " << n << " neglected because the corresponding boundary is too small" << std::endl;
        }
    }
    // Checks if the number of constraints is > than the actual size of the vector space
    auto tmp = initial;
    for (std::size_t idx = 0; idx < ortho_vecs_local.size(); idx++)
        tmp -= ietl::dot(tmp, ortho_vecs_local[idx]) * ortho_vecs_local[idx];
    if (tmp.scalar_norm() < thresholdForCompleteness)
        skipOptimization = true;
    // Actual Jacobi-Davidson diagonalization
    double tol = params["ietl_jcd_tol"];
    ietl::basic_iteration<double> iter(params["ietl_jcd_maxiter"], tol, tol);
    if (!skipOptimization) {
        SingleSiteVS<Matrix, SymmGroup> vs(initial, ortho_vecs_local);
        ietl::jcd_gmres_solver<SiteProblem<Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup> >
        jcd_gmres(sp, vs, params["ietl_jcd_gmres"]);
        ietl::jacobi_davidson<SiteProblem<Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup> >
        jd(sp, vs, ietl::Smallest);
        contraction::ContractionGrid<Matrix, SymmGroup>::iterate_reduction_layout(0, params["ietl_jcd_maxiter"]);
//        maquis::cout << "Ortho vecs " << ortho_vecs.size() << std::endl;
        for (int n = 0; n < ortho_vecs_local.size(); ++n) {
//            maquis::cout << "Ortho norm " << n << ": " << ietl::two_norm(ortho_vecs[n]) << std::endl;
            maquis::cout << "Input <MPS|O[" << n << "]> : " << ietl::dot(initial, ortho_vecs_local[n]) << std::endl;
        }
        r0 = jd.calculate_eigenvalue(initial, jcd_gmres, iter);
        for (int n = 0; n < ortho_vecs_local.size(); ++n)
            maquis::cout << "Output <MPS|O[" << n << "]> : " << ietl::dot(r0.second, ortho_vecs_local[n]) << std::endl;
    }
    else {
        maquis::cout << "Vector space too small, diagonalization skipped" << std::endl;
        auto sigmaVector = initial;
        ietl::mult(sp, initial, sigmaVector);
        auto energy = ietl::dot(initial, sigmaVector);
        r0 = std::make_pair(energy, initial);
    }
    maquis::cout << " Jacobi-Davidson diagonalization converged after " << iter.iterations() << " iterations." << std::endl;
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
