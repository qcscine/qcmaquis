/*****************************************************************************
 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations
 *
 * ALPS Libraries
 *
 * Copyright (C) 2017-2017 by Alberto Baiardi <alberto.baiardi@sns.it>
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

#ifndef VARIANCEOPTIMIZER_H
#define VARIANCEOPTIMIZER_H

#include <ietl/fmatrix.h>

template<class SymmGroup>
class VarianceOptimizer
{
public:
    // Types definition
    typedef typename ietl::FortranMatrix<double>     matrix_type ;
    typedef typename std::vector<double>             vector_ref ;
    typedef typename vector_ref::iterator            vector_iterator ;
private:
    // Attributes
    double              lagrange_multiplier_, scaling_, penalty_ ;
    matrix_type         H_, H_squared_, Hessian_, Hessian_inverse_ ;
    std::size_t         n_elements_, n_elements_squared_ ;
    vector_ref          coefficients_, coefficients_current_, gradient_ ;
    // Static members
    const bool          verbose_         = true ;
    const bool          add_lagrange_    = true ;
    const double        scaling_penalty_ = 10. ;
    const double        thresh_inner_    = 1.0E-10 ;
    const double        thresh_outer_    = 1.0E-10 ;
    const std::size_t   max_iter_inner_  = 100 ;
    const std::size_t   max_iter_outer_  = 100 ;
public:
    // Constructors
    VarianceOptimizer() ;
    VarianceOptimizer(matrix_type& H, matrix_type& H_squared) ;
    // Getter
    std::vector<double> get_coeff() ;
    // Optimization of the variance
    void PerformOptimization() ;
private:
    vector_ref apply_operator(const vector_ref& input, const bool& is_squared, const bool& is_inverted) ;
    double compute_derivatives() ;
    double compute_distance(vector_ref const& vec1, vector_ref const& vec2) ;
    double compute_norm(const bool& is_coeff) ;
    void invert_Hessian() ;
};

#include "varianceoptimizer.hpp"

#endif


