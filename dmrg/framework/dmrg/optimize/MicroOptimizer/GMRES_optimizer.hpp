/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2017 by Alberto Baiardi <alberto.baiardi@sns.it>
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

#ifndef MAQUIS_DMRG_GMRES_OPTIMIZER_H
#define MAQUIS_DMRG_GMRES_OPTIMIZER_H

#include "dmrg/optimize/MicroOptimizer/microoptimizer.h"
#include "dmrg/optimize/MicroOptimizer/optimizationalgorithm.h"

template<class MATRIX, class VectorSpace, class CorrectionEquation>
class GMRES_optimizer : public OptimizationAlgorithm<MATRIX, VectorSpace, CorrectionEquation>
{
private:
    // Types declaration
    typedef OptimizationAlgorithm<MATRIX, VectorSpace, CorrectionEquation>   base ;
    typedef typename base::scalar_type                                       scalar_type ;
    typedef typename base::vector_type                                       vector_type ;
    typedef typename base::size_t                                            size_t ;
    typedef typename base::matrix_scalar                                     matrix_scalar ;
    typedef typename base::vector_scalar                                     vector_scalar ;
    typedef typename base::vector_space                                      vector_space ;
    // Attributes
    using base::abs_tol_ ;
    using base::b_ ;
    using base::rel_tol_ ;
    using base::max_iter_ ;
    using base::n_restart_ ;
    using base::verbosity_ ;
public:
    // Constructor
    GMRES_optimizer() : base::OptimizationAlgorithm() {} ;
    GMRES_optimizer(const float& abs_error, const float& rel_error, const std::size_t& max_iter, const std::size_t& n_restart)
            : base::OptimizationAlgorithm(abs_error, rel_error, max_iter, n_restart) {} ;
    // Overriding of the perform optimization algorithm
    vector_type perform_optimization(MicroOptimizer<MATRIX, VectorSpace, CorrectionEquation>* optimizer,
                                     const vector_type& x0) ;
private:
    vector_scalar Update(const matrix_scalar& H, const vector_scalar& S, const size_t& k) ;
    void ApplyPlaneRotation(scalar_type& dx, scalar_type& dy, const scalar_type& cs, const scalar_type& sn) ;
    void GeneratePlaneRotation(scalar_type& dx, scalar_type& dy, scalar_type& cs, scalar_type& sn) ;
};

// Routine performing the optimization
template<class MATRIX, class VectorSpace, class CorrectionEquation>
typename GMRES_optimizer<MATRIX, VectorSpace, CorrectionEquation>::vector_type
         GMRES_optimizer<MATRIX, VectorSpace, CorrectionEquation>::perform_optimization(MicroOptimizer<MATRIX, VectorSpace, CorrectionEquation>* optimizer,
                                                                                        const vector_type& x0)
{
    // Types definiton
    double 	      normb = ietl::two_norm(b_) ;
    vector_type   r = x0, w , init, jnk ;
    vector_scalar s(max_iter_+1), cs(max_iter_+1), sn(max_iter_+1) , y ;
    vector_space  v(max_iter_+1) ;
    matrix_scalar H(max_iter_+1, max_iter_+1);
    // Initialization
    init = r;
    v[0] = optimizer->get_correction().apply_correction(init) ;
    r = b_ - v[0] ;
    s[0] = ietl::two_norm(r);
    if (verbosity_)
        std::cout << "GMRES - initial error " << s[0] << std::endl ;
    v[0] = r / s[0];
    size_t idx = 0 ;
    for (size_t i = 0 ; i < max_iter_ ; ++i)
    {
        jnk = v[idx] ;
        // Right preconditioning
        optimizer->get_correction().apply_precondition(jnk) ;
        w = optimizer->get_correction().apply_correction(jnk) ;
        // Update of the H matrix (Hessenberg, so nearly upper diagonal)
        for (std::size_t k = 0; k <= idx; ++k) {
            H(k,idx) = ietl::dot(w, v[k]);
            w -= H(k,idx) * v[k];
        }
        H(idx+1, idx) = ietl::two_norm(w);
        v[idx+1] = w / H(idx+1,idx);
        // Least-squares minimization
        for (std::size_t k = 0; k < idx; ++k)
            ApplyPlaneRotation(H(k,idx), H(k+1,idx), cs[k], sn[k]);
        GeneratePlaneRotation(H(idx,idx), H(idx+1,idx), cs[idx], sn[idx]);
        ApplyPlaneRotation(H(idx,idx), H(idx+1,idx), cs[idx], sn[idx]);
        ApplyPlaneRotation(s[idx], s[idx+1], cs[idx], sn[idx]);
        if (verbosity_)
          std::cout << "GMRES iteration " << idx << ", Abs. Err. = " << std::abs(s[idx+1])
                    << ", Rel. Err. = " << std::abs(s[idx+1])/normb << std::endl;
        bool converged = std::abs(s[idx+1])/normb < rel_tol_ || std::abs(s[idx+1]) < abs_tol_ || i == max_iter_-1;
        if (converged || idx+1 == n_restart_) {
            y = Update(H, s, idx) ;
            r = init ;
            for (std::size_t k = 0; k <= idx; ++k) {
                jnk = v[k] ;
                optimizer->get_correction().apply_precondition(jnk) ;
                r += y[k] * jnk;
            }
            if (converged) {
                return r ;
            } else {
                // Clean the matrices
                for (std::size_t i1 = 0; i1 <= idx+1; i1++) {
                    for (std::size_t i2 = 0; i2 <= idx+1; i2++)
                        H(i1,i2) = 0. ;
                    cs[i1] = 0. ;
                    sn[i1] = 0. ;
                    s[i1]  = 0. ;
                }
                // Uses the last approximation of the function as a guess
                idx = 0 ;
                init = r ;
                v[0] = optimizer->get_correction().apply_correction(init) ;
                r = b_ - v[0] ;
                s[0] = ietl::two_norm(r);
                if (verbosity_)
                    std::cout << "GMRES - initial error " << s[0] << std::endl ;
                v[0] = r / s[0];
                continue ;
            }
        }
        idx ++ ;
    }
    return r;
}

// Auxiliary routines

template<class MATRIX, class VectorSpace, class CorrectionEquation>
void GMRES_optimizer<MATRIX, VectorSpace, CorrectionEquation>::GeneratePlaneRotation(scalar_type &dx,
                                                                                     scalar_type &dy,
                                                                                     scalar_type &cs,
                                                                                     scalar_type &sn)
{
    if ( dy == 0. ) {
        cs = 1. ;
        sn = 0. ;
    } else if (std::abs(dy) > std::abs(dx)) {
        scalar_type tmp = dx / dy;
        sn = 1. / sqrt( 1. + tmp*tmp );
        cs = tmp*sn;
    } else {
        scalar_type tmp = dy / dx;
        cs = 1. / sqrt( 1. + tmp*tmp );
        sn = tmp*cs;
    }
}

template<class MATRIX, class VectorSpace, class CorrectionEquation>
void GMRES_optimizer<MATRIX, VectorSpace, CorrectionEquation>::ApplyPlaneRotation(scalar_type& dx,
                                                                                  scalar_type& dy,
                                                                                  const scalar_type& cs,
                                                                                  const scalar_type& sn)
{
    scalar_type r0  = cs*dx + sn*dy ;
    scalar_type r1 = -sn*dx + cs*dy ;
    dx = r0;
    dy = r1;
}

template<class MATRIX, class VectorSpace, class CorrectionEquation>
typename GMRES_optimizer<MATRIX, VectorSpace, CorrectionEquation>::vector_scalar
         GMRES_optimizer<MATRIX, VectorSpace, CorrectionEquation>::Update(matrix_scalar const& H,
                                                                          vector_scalar const& S,
                                                                          size_t const& k)
{
    vector_scalar y(S.begin(), S.begin()+k+1);
    for (int i = k; i >= 0; --i) {
        y[i] /= H(i, i);
        for (int j = i - 1; j >= 0; --j)
            y[j] -= H(j, i) * y[i];
    }
    return y;
}

#endif //MAQUIS_DMRG_GMRES_OPTIMIZER_H
