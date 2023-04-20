/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/factories/factory_lattice.hpp"

#ifdef ENABLE_ALPS_MODELS
#include "dmrg/models/alps/lattice.hpp"
#endif

#ifdef ENABLE_LL_MODELS
#include "dmrg/models/ll/ll_models.h"
#endif

/// lattice factory
std::shared_ptr<lattice_impl>
lattice_factory(BaseParameters & parms)
{
    typedef std::shared_ptr<lattice_impl> impl_ptr;

    if (parms["lattice_library"] == "coded") {
        return coded_lattice_factory(parms);
    } else if (parms["lattice_library"] == "alps") {
#ifdef ENABLE_ALPS_MODELS
        return impl_ptr( new alps_lattice(parms) );
#else
        throw std::runtime_error("This code was compiled without alps lattice.");
#endif
#ifdef ENABLE_LL_MODELS
    } else if (parms["lattice_library"] == "ll") {
        return ll_lattice_factory(parms);
#endif
    } else {
        throw std::runtime_error("Don't know this lattice_library!");
    }
}
