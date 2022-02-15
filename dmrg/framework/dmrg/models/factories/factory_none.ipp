/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2012 by Michele Dolfi <dolfim@phys.ethz.ch>
 *               2012      by Jan Gukelberger <gukelberger@phys.ethz.ch>
 *               2021      by Alberto Baiardi <abaiardi@ethz.ch>
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

#include "dmrg/models/factories/factory.h"
#include "dmrg/models/vibrational/none/model.hpp"


template<class Matrix>
struct coded_model_factory<Matrix, TrivialGroup> {
    // Types definition
    using PointerType = std::shared_ptr<model_impl<Matrix, TrivialGroup> >;
    // Factory class
    static PointerType parse(Lattice const& lattice, BaseParameters & parms)
    {
        if (parms["MODEL"] == std::string("watson"))
            return PointerType( new WatsonHamiltonian<Matrix>(lattice, parms, false));
        else {
            throw std::runtime_error("Don't know this model with None symmetry group!");
            return PointerType();
        }
    }
};
