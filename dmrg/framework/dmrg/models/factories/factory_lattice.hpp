/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MAQUIS_DMRG_MODELS_CODED_FACTORY_LATTICE_H
#define MAQUIS_DMRG_MODELS_CODED_FACTORY_LATTICE_H

#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/lattice/ChainLattice.hpp"
#include "dmrg/models/lattice/SquareLattice.hpp"
#include "dmrg/models/lattice/OrbitalLattice.hpp"
#include "dmrg/models/lattice/PreBOLattice.hpp"
#include "dmrg/models/lattice/WatsonLattice.hpp"
#include "dmrg/models/lattice/VibronicLattice.hpp"

/**
 * @brief Factory method returning the requested lattice
 * @param parms parameter container
 * @return std::shared_ptr<lattice_impl> pointer to the base class storing the 
 * specific lattice that has been requested
 */
inline std::shared_ptr<lattice_impl> coded_lattice_factory(BaseParameters & parms)
{
    typedef std::shared_ptr<lattice_impl> impl_ptr;
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
    else if (parms["LATTICE"] == std::string("spinors"))
        return impl_ptr(new Orbitals(parms));
#ifdef DMRG_PREBO
    else if (parms["LATTICE"] == std::string("preBO lattice"))
        return impl_ptr(new PreBOLattice(parms));
#endif
#ifdef DMRG_VIBRATIONAL
    else if (parms["LATTICE"] == std::string("watson lattice"))
        return impl_ptr(new WatsonLattice(parms));
#endif
#ifdef DMRG_VIBRONIC
    else if (parms["LATTICE"] == std::string("vibronic lattice"))
        return impl_ptr(new VibronicLattice(parms));
#endif
    else {
        throw std::runtime_error("Don't know this lattice!");
    }
}

#endif
