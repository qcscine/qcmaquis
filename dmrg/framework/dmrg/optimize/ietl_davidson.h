/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef IETL_DAVIDSON_SOLVER_H
#define IETL_DAVIDSON_SOLVER_H

#include "dmrg/utils/BaseParameters.h"

#include "ietl_lanczos_solver.h"

#include "davidson.h"

namespace davidson_detail {

    template<class Matrix, class SymmGroup, class = void>
    class ref_diag
    {
    public:
        void operator()(SiteProblem<Matrix, SymmGroup> const & H, MPSTensor<Matrix, SymmGroup> x){}
    };

    template<class Matrix, class SymmGroup>
    class ref_diag<Matrix, SymmGroup, symm_traits::enable_if_su2_t<SymmGroup>>
    {
    public:
        void operator()(SiteProblem<Matrix, SymmGroup> const & H, MPSTensor<Matrix, SymmGroup> x)
        {
            typedef typename Matrix::value_type value_type;
            typedef typename SymmGroup::charge charge;

            x.make_left_paired();
            x.multiply_by_scalar(.0);
            block_matrix<Matrix, SymmGroup> & bm = x.data();


            block_matrix<Matrix, SymmGroup> ret2 = contraction::SU2::diagonal_hamiltonian(H.left, H.right, H.mpo, x);

            for (size_t b = 0; b < bm.n_blocks(); ++b)
            {
                maquis::cout << bm.basis().left_charge(b) << bm.basis().right_charge(b) << std::endl;
                for (size_t i = 0; i < num_rows(bm[b]); ++i)
                for (size_t j = 0; j < num_cols(bm[b]); ++j)
                {
                    bm[b](i,j) = 1;
                    MPSTensor<Matrix, SymmGroup> prod;
                    ietl::mult(H, x, prod);
                    maquis::cout << "  " << i << "," << j << "  " << prod.data()[b](i,j) << " " << ret2[b](i,j) << std::endl;
                    if (std::abs(prod.data()[b](i,j) - ret2[b](i,j)) > 1e-12) throw std::runtime_error("diagonal wrong\n");
                    bm[b](i,j) = 0;
                }
            }
        }
    };

    template<class Matrix, class SymmGroup, class = void>
    class MultDiagonal
    {
        typedef MPSTensor<Matrix, SymmGroup> vector_type;
        typedef typename Matrix::value_type value_type;

    public:

        MultDiagonal(SiteProblem<Matrix, SymmGroup> const& H, vector_type const& x)
        {
            throw std::runtime_error("Davidson only implemented for spin-adapted Hamiltonians\n");
        }

        void precondition(vector_type& r, vector_type& V, value_type theta)
        {
        }
    };

    template<class Matrix, class SymmGroup>
    class MultDiagonal<Matrix, SymmGroup, symm_traits::enable_if_su2_t<SymmGroup>>
    {
        typedef MPSTensor<Matrix, SymmGroup> vector_type;
        typedef typename Matrix::value_type value_type;

    public:

        MultDiagonal(SiteProblem<Matrix, SymmGroup> const& H, vector_type const& x)
        {
            Hdiag = contraction::SU2::diagonal_hamiltonian(H.left, H.right, H.mpo, x);
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

        void mult_diag(value_type theta, vector_type& x)
        {
            block_matrix<Matrix, SymmGroup> & data = x.data();

            assert(shape_equal(data, Hdiag));
            for (size_t b = 0; b < data.n_blocks(); ++b)
            {
                for (size_t i = 0; i < num_rows(data[b]); ++i)
                    for (size_t j = 0; j < num_cols(data[b]); ++j)
                        if (std::abs(theta - Hdiag[b](i,j)))
                            data[b](i,j) /= (theta - Hdiag[b](i,j));
            }
        }

        block_matrix<Matrix, SymmGroup> Hdiag;
    };

} // namespace davidson detail

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
