/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2012 by Andreas Hehn <hehn@phys.ethz.ch>                          *
 *                                                                                 *
 * This software is part of the ALPS libraries, published under the ALPS           *
 * Library License; you can use, redistribute it and/or modify it under            *
 * the terms of the license, either version 1 or (at your option) any later        *
 * version.                                                                        *
 *                                                                                 *
 * You should have received a copy of the ALPS Library License along with          *
 * the ALPS Libraries; see the file LICENSE.txt. If not, the license is also       *
 * available from http://alps.comp-phys.org/.                                      *
 *                                                                                 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,        *
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT       *
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE       *
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,     *
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER     *
 * DEALINGS IN THE SOFTWARE.                                                       *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ALPS_NUMERIC_MATRIX_ADJOINT_HPP
#define ALPS_NUMERIC_MATRIX_ADJOINT_HPP

#include <alps/numeric/matrix/adjoint_view.hpp>

namespace alps {
namespace numeric {

template <typename Matrix>
inline adjoint_view<Matrix> adjoint(Matrix const& m) {
    return adjoint_view<Matrix>(m);
}

template <typename Matrix>
void adjoint_inplace(Matrix& m) {
    typedef typename Matrix::size_type size_type;
    using std::swap;
    if(num_rows(m) == num_cols(m) ) {
        for(size_type i = 0; i < num_rows(m); ++i) {
            m(i, i) = std::conj(m(i, i));
            for(size_type j = i+1; j < num_cols(m); ++j) {
                m(i, j) = std::conj(m(i, j));
                m(j, i) = std::conj(m(j, i));
                swap(m(i,j), m(j,i));
            }
        }
    }
    else {
        // TODO replace this code by an actual inplace implementation
        Matrix m2 = adjoint(m);
        swap(m,m2);
    }
}

} // end namespace numeric
} // end namespace alps

#endif //ALPS_NUMERIC_MATRIX_ADJOINT_HPP
