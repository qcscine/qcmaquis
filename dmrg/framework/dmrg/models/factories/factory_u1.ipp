/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *               2022- by Alberto Baiardi <abaiardi@ethz.ch>
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
#include "dmrg/models/vibrational/u1/VibronicModel.hpp"
#include "dmrg/models/vibrational/u1/ExcitonicModel.hpp"

template<class Matrix>
struct coded_model_factory<Matrix, U1> {
    using PointerType = std::shared_ptr<model_impl<Matrix, U1> >;

    /** @brief Factory function for the model class */
    static PointerType parse(Lattice const& lattice, BaseParameters & parms)
    {
        typedef std::shared_ptr<model_impl<Matrix, U1> > impl_ptr;
        if (parms["MODEL"] == std::string("vibronic")) {
#ifdef DMRG_VIBRONIC
            return impl_ptr( new VibronicModel<Matrix>(lattice, parms));
#else
            throw std::runtime_error("Don't know this model!");
#endif
        }
        else if (parms["MODEL"] == std::string("excitonic")) {
#ifdef DMRG_VIBRONIC
            return impl_ptr( new HolsteinHubbardExcitonicHamiltonian<Matrix>(lattice, parms));
#else
            throw std::runtime_error("Don't know this model!");
#endif
        }
        else {
            throw std::runtime_error("Don't know this model!");
            return impl_ptr();
        }
    }
};