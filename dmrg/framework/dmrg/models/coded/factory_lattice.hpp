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

#ifndef MAQUIS_DMRG_MODELS_CODED_FACTORY_LATTICE_H
#define MAQUIS_DMRG_MODELS_CODED_FACTORY_LATTICE_H

#include "dmrg/models/coded/lattice.hpp"

inline boost::shared_ptr<lattice_impl>
coded_lattice_factory(BaseParameters & parms)
{
    typedef boost::shared_ptr<lattice_impl> impl_ptr;
    if (parms["LATTICE"] == std::string("periodic chain lattice"))
        return impl_ptr(new ChainLattice(parms, true));
    else if (parms["LATTICE"] == std::string("chain lattice"))
        return impl_ptr(new ChainLattice(parms, false));
    else if (parms["LATTICE"] == std::string("open chain lattice"))
        return impl_ptr(new ChainLattice(parms, false));
    else if (parms["LATTICE"] == std::string("square lattice"))
        return impl_ptr(new SquareLattice(parms));
    else if (parms["LATTICE"] == std::string("open square lattice"))
        return impl_ptr(new SquareLattice(parms));
    else if (parms["LATTICE"] == std::string("orbitals"))
        return impl_ptr(new Orbitals(parms));
    else {
        throw std::runtime_error("Don't know this lattice!");
    }
}

#endif
