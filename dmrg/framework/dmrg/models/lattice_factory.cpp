/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#include "dmrg/models/lattice.h"
#include "dmrg/models/coded/factory_lattice.hpp"
#include "dmrg/models/continuum/factory_lattice.hpp"

#ifdef ENABLE_ALPS_MODELS
#include "dmrg/models/alps/lattice.hpp"
#endif

#ifdef ENABLE_LL_MODELS
#include "dmrg/models/ll/ll_models.h"
#endif

/// lattice factory
boost::shared_ptr<lattice_impl>
lattice_factory(BaseParameters & parms)
{
    typedef boost::shared_ptr<lattice_impl> impl_ptr;
    
    if (parms["lattice_library"] == "coded") {
        return coded_lattice_factory(parms);
    } else if (parms["lattice_library"] == "alps") {
#ifdef ENABLE_ALPS_MODELS
        return impl_ptr( new alps_lattice(parms) );
#else
        throw std::runtime_error("This code was compiled without alps lattice.");
#endif
    } else if (parms["lattice_library"] == "continuum") {
        return cont_lattice_factory(parms);
#ifdef ENABLE_LL_MODELS
    } else if (parms["lattice_library"] == "ll") {
        return ll_lattice_factory(parms);
#endif
    } else {
        throw std::runtime_error("Don't know this lattice_library!");
    }
}

