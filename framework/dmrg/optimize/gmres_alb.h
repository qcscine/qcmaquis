/*****************************************************************************
 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations
 *
 * ALPS Libraries
 *
 * Copyright (C) 2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *               2017 by Alberto Baiardi <alberto.baiardi@sns.it>
 *
 * This software is part of the ALPS libraries, published under the ALPS
 * Library License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 *
 * You should have received a copy of the ALPS Library License along with
 * the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
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

#ifndef IETL_GMRES_ALB_H
#define IETL_GMRES_ALB_H

#include <vector>
#include <iostream>

//
// This code is based on IML++, http://math.nist.gov/iml++/
// Modified by Alberto Baiardi to fix an inconsistency in the
// definition of the indexes in the original alps version
//

namespace ietl
{
    //
    // GMRES_GENERAL CLASS
    // -------------------
    // General class with the GMRES solver. The definition of the matrix of the linear system
    // is left free, since it differs for standard and shift-and-invert problem
    template<class Matrix, class Vector, class VectorSpace>
    class gmres_general {
    //
    public:
        // Types definition
        typedef typename Vector::value_type                            scalar_type ;
        typedef typename std::size_t                                   size_t ;
        typedef typename boost::numeric::ublas::matrix< scalar_type >  matrix_scalar ;
        typedef typename std::vector<scalar_type>                      vector_scalar ;
        typedef typename std::vector<Vector>                           vector_space ;
        typedef typename std::vector< std::pair<Vector,Vector> >       vector_results ;
        //
        // Public methods
        // --------------
        // Constructor
        gmres_general(Matrix const & A, Vector const & u, VectorSpace const & vs, double const & theta,
                      const vector_results& v_res, size_t const& n_root, std::size_t max_iter = 100, bool verbose = false)
                : max_iter(max_iter), verbose(verbose), A_(A), u_(u), theta_(theta), vs_(vs), n_root_(n_root), v_res_(v_res)
        { }
        // () operator. Returns the solution of the GMRES problem
        Vector operator()(Vector const & b,
                          Vector const & x0,
                          double abs_tol = 1e-6)
        {
            // Types definiton
            Vector r, w, init = x0 ;
            vector_scalar s(max_iter+1), cs(max_iter+1), sn(max_iter+1);
            vector_scalar y;
            vector_space v(max_iter+1);
            // Initialization
            v[0] = apply(init) ;
            r = b - v[0];
            s[0] = two_norm(r);
            //if (std::abs(s[0]) < abs_tol) {
            //    if (verbose)
            //        std::cout << "Already done with x0." << std::endl;
            //    return x0;
            //}
            v[0] = r / s[0];
            matrix_scalar H(max_iter+1, max_iter+1);
            size_t i = 0 ;
            for ( ; i < max_iter ; ++i)
            {
                w = apply(v[i]);
                // Update of the H matrix (Hessenberg, so nearly upper diagonal)
                for (std::size_t k = 0; k <= i; ++k) {
                    H(k,i) = dot(w, v[k]);
                    w -= H(k,i) * v[k];
                }
                H(i+1, i) = two_norm(w);
                v[i+1] = w / H(i+1,i);
                // Least-squares minimization
                for (std::size_t k = 0; k < i; ++k)
                    ApplyPlaneRotation(H(k,i), H(k+1,i), cs[k], sn[k]);
                GeneratePlaneRotation(H(i,i), H(i+1,i), cs[i], sn[i]);
                ApplyPlaneRotation(H(i,i), H(i+1,i), cs[i], sn[i]);
                ApplyPlaneRotation(s[i], s[i+1], cs[i], sn[i]);
                if (verbose)
                    std::cout << "GMRESAlb iteration " << i << ", resid = " << std::abs(s[i+1]) << std::endl;
                if (std::abs(s[i+1]) < abs_tol) {
                    y = Update(H, s, i);
                    r = x0;
                    for (std::size_t k = 0; k <= i; ++k)
                        r += y[k] * v[k];
                    return r;
                }
            }
            y = Update(H, s, i-1);
            r = x0;
            for (std::size_t k = 0; k <= i-1; ++k)
                r += y[k] * v[k];
            return r;
        }
    private:
        //
        // Compute the sine and cosine of the angle to be used in the successive 2x2 rotations
        // that are performed in the GMRES algorithm
        void GeneratePlaneRotation(scalar_type dx, scalar_type dy, scalar_type & cs, scalar_type & sn)
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
        //
        // Applies the rotation generated bt GeneratePlaneRotation
        void ApplyPlaneRotation(scalar_type & dx, scalar_type & dy, scalar_type cs, scalar_type sn)
        {
            scalar_type r0  = cs*dx + sn*dy ;
            scalar_type r1 = -sn*dx + cs*dy ;
            dx = r0;
            dy = r1;
        }
        //
        // Solve the inverse problem
        vector_scalar Update(matrix_scalar const & H, vector_scalar const & S, size_t k)
        {
            vector_scalar y(S.begin(), S.begin() + k + 1);
            for (int i = k; i >= 0; --i) {
                y[i] /= H(i, i);
                for (int j = i - 1; j >= 0; --j)
                    y[j] -= H(j, i) * y[i];
            }
            return y;
        }
    protected:
        // Compute the result of the application of the operator
        // THis function is left virtual since its definition is different for standard
        // and shift-and-invert problems
        virtual Vector apply(Vector& input) {} ;
        //
        // Private attributes
        // ------------------
        // A_       --> Matrix to be diagonalized
        // u_       --> Vector used to determine the orthogonal subspace
        // max_iter --> Maximum number of iterations for the GMRES algorithm
        // vebose   --> Controls the output printing
        Matrix A_ ;
        Vector u_ ;
        std::size_t max_iter , n_root_ ;
        bool verbose ;
        double theta_ ;
        VectorSpace vs_ ;
        vector_results v_res_ ;
    };
    //
    // GMRES_STANDARD OBJECT
    // ---------------------
    // GMRES solver object for the standard problem
    template<class Matrix, class Vector, class VectorSpace>
    class gmres_standard : private gmres_general<Matrix, Vector, VectorSpace>
    {
    public:
        typedef gmres_general<Matrix, Vector, VectorSpace> base ;
        typedef typename base::size_t          size_t ;
        typedef typename base::vector_results  vector_results ;
        //
        gmres_standard(Matrix const & A,
                       Vector const & u,
                       VectorSpace const & vs,
                       double const & theta,
                       vector_results const & v_res,
                       size_t const& n_root,
                       size_t max_iter,
                       bool verbose)
        : base::gmres_general(A, u, vs, theta, v_res, n_root, max_iter, verbose) { }
        // Private attributes
        using base::A_ ;
        using base::n_root_ ;
        using base::theta_ ;
        using base::u_ ;
        using base::v_res_ ;
        using base::vs_ ;
        using base::operator() ;
    private:
        Vector apply(Vector& input){
            // Initialization
            Vector t, t2, t3, y;
            // t2 = (1-uu*) x
            double ust = dot(u_, input) ;
            t2 = input - ust * u_;
            for (typename vector_results::iterator it = v_res_.begin(); it != v_res_.end(); it++)
                t2 -= ietl::dot((*it).first, t2)/ietl::dot((*it).first,(*it).first) * (*it).first;
            // y = (A-theta*1) t2
            ietl::mult(A_, t2, t3, n_root_);
            y = t3 - theta_ * t2;
            for (typename vector_results::iterator it = v_res_.begin(); it != v_res_.end(); it++)
                y -= ietl::dot((*it).first, y)/ietl::dot((*it).first,(*it).first) * (*it).first;
            // t = (1-uu*) y
            ust = dot(u_, y) ;
            t = y - ust * u_;
            y = t ;
            // Finalization
            return y ;
        }
    };
    //
    // GMRES_MODIFIED OBJECT
    // ---------------------
    // GMRES solver object for the shift-and-inverted problem
    template<class Matrix, class Vector, class VectorSpace>
    class gmres_modified : private gmres_general<Matrix, Vector, VectorSpace>
    {
    public:
        typedef gmres_general<Matrix, Vector, VectorSpace> base ;
        typedef typename base::size_t          size_t ;
        typedef typename base::vector_results  vector_results ;
        gmres_modified(Matrix const & A,
                       Vector const & u,
                       VectorSpace const & vs,
                       Vector const & z,
                       double const & theta,
                       vector_results const & v_res,
                       size_t const& n_root,
                       double const & omega,
                       size_t max_iter,
                       bool verbose)
                : base::gmres_general(A, u, vs, theta, v_res, n_root, max_iter, verbose),
                  omega_(omega), z_(z) { }
        // Private attributes
        using base::A_ ;
        using base::n_root_ ;
        using base::theta_ ;
        using base::u_ ;
        using base::vs_ ;
        using base::v_res_ ;
        using base::operator() ;
    private:
        Vector apply(Vector& input){
            // Initialization
            Vector t, t2, t3, y ;
            mult(A_, input, t, n_root_);
            t *= -1.;
            t += omega_ * input;
            double ust = dot(z_, t) ;
            t2 = input - ust * u_;
            //ietl::project(t2, vs_ ) ;
            // y = (A-theta*1) t2
            mult(A_, t2, t3, n_root_);
            t3 *= -1.;
            t3 += omega_ * t2;
            y = t3 - t2 / theta_ ;
            //ietl::project(y, vs_) ;
            // t = (1-uu*) y
            ust = dot(z_, y) ;
            t = y - ust * z_ ;
            // Finalization
            return t ;
        }
        // Private attributes
        double omega_ ;
        Vector z_ ;
    };
    //
    // GMRES_MODIFIED_INITIALIZER OBJECT
    // ---------------------------------
    // GMRES solver object for the shift-and-inverted problem
    template<class Matrix, class Vector, class VectorSpace>
    class gmres_initializer_modified : private gmres_general<Matrix, Vector, VectorSpace>
    {
    public:
        typedef gmres_general<Matrix, Vector, VectorSpace> base ;
        typedef typename base::size_t          size_t ;
        typedef typename base::vector_results  vector_results ;
        gmres_initializer_modified(Matrix const & A,
                                   Vector const & u,
                                   VectorSpace const & vs,
                                   double const & theta,
                                   vector_results const & vec_res,
                                   size_t const & n_root,
                                   size_t max_iter,
                                   bool verbose)
                : base::gmres_general(A, u, vs, theta, vec_res, n_root, max_iter, verbose) {} ;
        // 
        using base::A_ ;
        using base::n_root_ ;
        using base::theta_ ;
        using base::u_ ;
        using base::vs_ ;
        using base::operator() ;
    private:
        // Private attributes
        Vector apply(Vector& input){
            // Initialization
            Vector t ;
            mult(A_, input, t, n_root_);
            t *= -1.;
            t += theta_*input;
            return t ;
        }
    } ;
}

#endif
