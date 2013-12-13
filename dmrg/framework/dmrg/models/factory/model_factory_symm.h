/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#include "matrices.h"

#include "dmrg/models/coded/factory.h"
//#include "dmrg/models/continuum/factory.h" // temporary fix until model can be patched to new format.
#include "dmrg/models/factory/initializer_factory.h"

#ifdef ENABLE_ALPS_MODELS
#include "dmrg/models/alps/model.hpp"
#endif

#ifdef ENABLE_LL_MODELS
#include "dmrg/models/ll/ll_models.h"
#endif

// init MACROS
#define impl_model_factory(MATRIX, SYMMGROUP)                                           \
template boost::shared_ptr<model_impl<MATRIX, SYMMGROUP> >                              \
model_factory<MATRIX,SYMMGROUP>(Lattice const&, BaseParameters &);

// Implementation
template <class Matrix, class SymmGroup>
boost::shared_ptr<model_impl<Matrix, SymmGroup> >
model_factory(Lattice const& lattice, BaseParameters & parms)
{
    typedef boost::shared_ptr<model_impl<Matrix, SymmGroup> > impl_ptr;
    if (parms["model_library"] == "alps") {
#ifdef ENABLE_ALPS_MODELS
        if (parms["lattice_library"] != "alps")
            throw std::runtime_error("ALPS models require ALPS lattice.");
        return impl_ptr( new ALPSModel<Matrix, SymmGroup>(lattice, parms) );
#else
        throw std::runtime_error("This code was compiled without alps models.");
#endif
    } else if (parms["model_library"] == "coded") {
        return coded_model_factory<Matrix, SymmGroup>::parse(lattice, parms);
//    } else if (parms["model_library"] == "continuum") {
//        return cont_model_factory<Matrix, SymmGroup>::parse(lattice, parms);
//#ifdef ENABLE_LL_MODELS
//    } else if (parms["model_library"] == "ll") {
//        return ll_model_factory<Matrix, SymmGroup>::parse(lattice, parms);
//#endif
    } else {
        throw std::runtime_error("Don't know this model_library!");
    }
    
}

