/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#include "dmrg/mp_tensors/boundary.h"
#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/mpotensor.h"
#include "dmrg/optimize/ietl_lanczos_solver.h"

#ifndef SITEPROBLEM
#define SITEPROBLEM

/**
 * @brief Class representing a site-centered eigenvalue problem.
 * 
 * This class wraps, in practice, all ingredients that are required to calculate
 * the sigma vector for a MPS/MPO contraction.
 */
template<class Matrix, class SymmGroup>
struct SiteProblem
{
    /** @brief Default constructor */
    SiteProblem(Boundary<typename storage::constrained<Matrix>::type, SymmGroup> const & left_,
                Boundary<typename storage::constrained<Matrix>::type, SymmGroup> const & right_,
                MPOTensor<Matrix, SymmGroup> const & mpo_) : left(left_), right(right_), mpo(mpo_)
    { }

    /** @brief Method to evaluate the sigma vector */
    MPSTensor<Matrix, SymmGroup> apply(const MPSTensor<Matrix, SymmGroup>& input_vec) const {
      MPSTensor<Matrix, SymmGroup> ret;
      ietl::mult(*this, input_vec, ret);
      return ret;
    }

    /** @brief Energy getter */
    typename MPSTensor<Matrix, SymmGroup>::real_type get_energy(const MPSTensor<Matrix, SymmGroup>& x)
    {
        MPSTensor<Matrix, SymmGroup> y;
        ietl::mult(*this, x, y);
        auto res = ietl::dot(x, y)/ietl::dot(x, x);
        x.make_left_paired();
        return maquis::real(res);
    }

    /** Class members */
    Boundary<typename storage::constrained<Matrix>::type, SymmGroup> const & left;
    Boundary<typename storage::constrained<Matrix>::type, SymmGroup> const & right;
    MPOTensor<Matrix, SymmGroup> const & mpo;
    double ortho_shift=0.;
};

namespace ietl {

template<class Matrix, class SymmGroup>
typename MPSTensor<Matrix, SymmGroup>::real_type get_energy(const SiteProblem<Matrix, SymmGroup>& H, MPSTensor<Matrix, SymmGroup> const &x)
{
    MPSTensor<Matrix, SymmGroup> y;
    ietl::mult(H, x, y);
    typename MPSTensor<Matrix, SymmGroup>::scalar_type res = ietl::dot(x,y)/ietl::dot(x,x);
    x.make_left_paired();
    return maquis::real(res);
}

} // ietl

#endif // SITEPROBLEM
