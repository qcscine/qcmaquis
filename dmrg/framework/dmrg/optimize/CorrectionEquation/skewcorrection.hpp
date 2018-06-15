/*****************************************************************************
 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations
 *
 * ALPS Libraries
 *
 * Copyright (C) 2017 by Alberto Baiardi <alberto.baiardi@sns.it>
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

#ifndef MAQUIS_DMRG_SKEWCORRECTION_HPP
#define MAQUIS_DMRG_SKEWCORRECTION_HPP

#include <dmrg/mp_tensors/mpstensor.h>
#include "dmrg/optimize/CorrectionEquation/corrector.cpp"

// +---------------------+
//  Skew corrector object
// +---------------------+

template<class MATRIX, class VecSpace>
class SkewCorrection : public Corrector<MATRIX, VecSpace> {
    // Types definition
    typedef Corrector<MATRIX, VecSpace>        base ;
    typedef typename base::scalar_type         scalar_type ;
    typedef typename base::vector_type         vector_type ;
    typedef typename vector_type::bm_type      preconditioner_type ;
public:
    // Routine used to apply the correction equation
    vector_type apply(CorrectionEquation<MATRIX, VecSpace>* corr_eq, const vector_type& input)
    {
        // Initialization
        vector_type t, t2, t3, y ;
        scalar_type ust = ietl::dot(corr_eq->get_Au(), input) ;
        t2 = input - ust * corr_eq->get_u() / ietl::dot(corr_eq->get_Au(),corr_eq->get_u()) ;
        corr_eq->orthogonalize_simple(t2) ;
        // y = (A-theta*1) t2
        ietl::mult(corr_eq->get_hamiltonian(), t2, t3, corr_eq->get_n_root(), false) ;
        y = corr_eq->get_omega()*t2 - t3 - corr_eq->get_rayleigh()*t2 ;
        corr_eq->orthogonalize_simple(y) ;
        // t = (1-uu*) y
        ust = ietl::dot(corr_eq->get_Au(), y) ;
        t = y - ust * corr_eq->get_u() / ietl::dot(corr_eq->get_Au(), corr_eq->get_u()) ;
        // Finalization
        return t ;
    }
    // Routine used to do precondition
    void precondition(CorrectionEquation<MATRIX, VecSpace>* corr_eq, vector_type& input)
    {
        vector_type u = corr_eq->get_u() ;
        vector_type uA = corr_eq->get_Au() ;
        multiply_diagonal(corr_eq, input) ;
        multiply_diagonal(corr_eq, uA) ;
        scalar_type alpha = ietl::dot(corr_eq->get_u(), input) / ietl::dot(corr_eq->get_u(), uA) ;
        input -= alpha * uA ;
    }
private:
    void multiply_diagonal(CorrectionEquation<MATRIX, VecSpace>* corr_eq, vector_type& input)
    {
        scalar_type denom ;
        preconditioner_type &data = input.data() ;
        assert(shape_equal(data, corr_eq->get_preconditioner())) ;
        for (size_t b = 0; b < data.n_blocks(); ++b) {
            for (size_t i = 0; i < num_rows(data[b]); ++i) {
                for (size_t j = 0; j < num_cols(data[b]); ++j) {
                    denom = (corr_eq->get_omega() - (corr_eq->get_preconditioner())[b](i, j)) - corr_eq->get_rayleigh() ;
                    if (std::abs(denom) > 1.0E-10)
                        data[b](i, j) /= denom ;
                }
            }
        }
    }
};

#endif
