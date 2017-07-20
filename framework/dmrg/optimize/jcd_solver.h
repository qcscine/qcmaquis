/*****************************************************************************
 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations
 *
 * ALPS Libraries
 *
 * Copyright (C) 2001-2011 by Rene Villiger <rvilliger@smile.ch>,
 *                            Prakash Dayal <prakash@comp-phys.org>,
 *                            Matthias Troyer <troyer@comp-phys.org>
 *                            Bela Bauer <bauerb@phys.ethz.ch>
 *               2017-2017 by Alberto Baiardi <alberto.baiardi@sns.it>
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

#include <ietl/traits.h>
#include <ietl/fmatrix.h>
#include <ietl/ietl2lapack.h>
#include <ietl/cg.h>
#include <complex>
#include <vector>
#include <boost/function.hpp>

#include "ietl_lanczos_solver.h"

#ifndef IETL_JCD_SOLVER
#define IETL_JCD_SOLVER

namespace ietl {
// +-------------------+
//  JCD_SOLVER_OPERATOR
// +-------------------+
// Object which is used as an interface between the JD optimizer for vDMRG and the gmres routine
// which is inside ietl/gmres.h
//
// Attributes:
// - vectors u_ and r_ (projector and error vectors)
// - matrix m_ to be diagonalized
// - theta_ , approximation to the eigenvalue
//
// Methods:
// jcd_sovler_operator(x,y) --> takes in input a vector x and solve the approximate
//                              JD equation to give a new estimate for the vector y
//

    template<class Matrix, class VS, class Vector>
    class jcd_solver_operator {
    public:
        typedef typename vectorspace_traits<VS>::vector_type vector_type;
        typedef typename vectorspace_traits<VS>::scalar_type scalar_type;
        typedef typename ietl::number_traits<scalar_type>::magnitude_type magnitude_type;
        jcd_solver_operator(const vector_type &u, const vector_type &r, const Matrix &m, const magnitude_type &theta)
                : u_(u), r_(r), m_(m), theta_(theta) {} ;
        virtual ~jcd_solver_operator() {};
    protected:
        vector_type u_, r_;
        Matrix m_;
        magnitude_type theta_;
    };
    // +-------------------------------+
    //  STANDARD SOLVER OPERATOR OBJECT
    // +-------------------------------+
    template<class Matrix, class VS, class Vector>
    class jcd_solver_operator_standard : public jcd_solver_operator<Matrix, VS, Vector> {
    public:
        typedef jcd_solver_operator<Matrix, VS, Vector> base;
        typedef typename base::vector_type vector_type;
        typedef typename base::scalar_type scalar_type;
        typedef typename base::magnitude_type magnitude_type;
        using base::u_ ;
        using base::r_ ;
        using base::m_ ;
        using base::theta_ ;
        // -- Constructor --
        jcd_solver_operator_standard(const vector_type &u, const vector_type &r, const Matrix &m,
                                     const magnitude_type &theta) : base::jcd_solver_operator(u, r, m, theta) { };
        ~jcd_solver_operator_standard() {};
        void operator()(vector_type const &x, vector_type &y) const;
    };
    //
    template<class Matrix, class VS, class Vector>
    void jcd_solver_operator_standard<Matrix, VS, Vector>::operator()(vector_type const &x, vector_type &y) const
    {
        vector_type t, t2, t3;
        // t2 = (1-uu*) x
        scalar_type ust = dot(u_, x);
        t2 = x - ust * u_;
        // y = (A-theta*1) t2
        ietl::mult(m_, t2, t3);
        y = t3 - theta_ * t2;
        // t = (1-uu*) y
        ust = dot(u_, y);
        t = y - ust * u_;
        y = t ;
    }
    // +-------------------------------+
    //  MODIFIED SOLVER OPERATOR OBJECT
    // +-------------------------------+
    template<class Matrix, class VS, class Vector>
    class jcd_solver_operator_modified : public jcd_solver_operator<Matrix, VS, Vector> {
    public:
        typedef jcd_solver_operator<Matrix, VS, Vector> base ;
        typedef typename base::vector_type vector_type ;
        typedef typename base::scalar_type scalar_type ;
        typedef typename base::magnitude_type magnitude_type ;
        using base::u_;
        using base::r_;
        using base::m_;
        using base::theta_;
        // -- Constructor --
        jcd_solver_operator_modified(const vector_type &u, const vector_type &r, const Matrix &m,
                                     const magnitude_type &theta, const magnitude_type &omega,
                                     const vector_type &z)
                                     : base::jcd_solver_operator(u, r, m, theta), omega_(omega), z_(z) { };
        ~jcd_solver_operator_modified() {};
        void operator()(vector_type const &x, vector_type &y) const;
    private:
        magnitude_type omega_;
        vector_type z_ ;
    };
    //
    template<class Matrix, class VS, class Vector>
    void jcd_solver_operator_modified<Matrix, VS, Vector>::operator()(vector_type const &x, vector_type &y) const
    {
        vector_type t, t2, t3;
        mult(m_, x, t);
        t *= -1.;
        t += omega_ * x;
        scalar_type ust = dot(z_, t);
        t2 = x - ust * u_;
        // y = (A-theta*1) t2
        mult(m_, t2, t3);
        t3 *= -1.;
        t3 += omega_ * t2;
        y = t3 - t2 * theta_;
        // t = (1-uu*) y
        ust = dot(z_, y);
        t = y - ust * z_;
        y = t ;
    };
    // Multiplication of solvers
    template<class Matrix, class VS, class Vector>
    void mult(jcd_solver_operator_standard<Matrix, VS, Vector> const & m,
              typename jcd_solver_operator_standard<Matrix, VS, Vector>::vector_type const & x,
              typename jcd_solver_operator_standard<Matrix, VS, Vector>::vector_type & y)
    {
        m(x,y);
    }
    //
    template<class Matrix, class VS, class Vector>
    void mult(jcd_solver_operator_modified<Matrix, VS, Vector> const & m,
              typename jcd_solver_operator_modified<Matrix, VS, Vector>::vector_type const & x,
              typename jcd_solver_operator_modified<Matrix, VS, Vector>::vector_type & y)
    {
        m(x,y);
    }
}

#endif
