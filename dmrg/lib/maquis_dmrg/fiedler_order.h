/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
*               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
*               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
*               2014-2014    Sebastian Keller <sebkelle@phys.ethz.ch>
*               2019         Leon Freitag <lefreita@ethz.ch>
*               2019         Stefan Knecht <stknecht@ethz.ch>
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
#ifndef FIEDLER_ORDER_H
#define FIEDLER_ORDER_H

#include "maquis_dmrg.h"

namespace maquis
{
    template <class V> // real or complex
    class FiedlerOrder
    {
        public:
            // typedef for measurements
            typedef maquis::meas_with_results_type<V> meas_with_results_type;

            FiedlerOrder(int nstates_,
                         const std::string& pname);

            ~FiedlerOrder();

            // calculate Fiedler order
            void get_FiedlerOrder();

        private:

            typedef typename alps::numeric::matrix<V> Matrix;

            // Number of states to include in Fiedler vector calculation
            int nstates_;

            // Project name
            std::string pname_;

            // calculate the Fiedler order for a given (state-averaged) mutual information matrix
            std::string calculate_FiedlerOrder();

            Matrix get_laplacian(Matrix mutI);

            struct Impl;
            std::unique_ptr<Impl> impl_;
    };
}

#endif
