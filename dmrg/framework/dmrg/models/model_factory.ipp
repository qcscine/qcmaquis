/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#include "dmrg/models/model.h"
#include "dmrg/models/factories/factory.h"
#include "dmrg/models/initializer_factory.h"

#ifdef ENABLE_ALPS_MODELS
#include "dmrg/models/alps/model.hpp"
#endif

#ifdef ENABLE_LL_MODELS
#include "dmrg/models/ll/ll_models.h"
#endif

// init MACROS
#define impl_model_factory(MATRIX,SYMMGROUP)                                          \
template std::shared_ptr<model_impl<MATRIX, SYMMGROUP> >                              \
model_factory<MATRIX,SYMMGROUP>(Lattice const&, BaseParameters &);

// Implementation
template <class Matrix, class SymmGroup>
std::shared_ptr<model_impl<Matrix, SymmGroup> >
model_factory(Lattice const& lattice, BaseParameters & parms)
{
    typedef std::shared_ptr<model_impl<Matrix, SymmGroup> > impl_ptr;
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
    } else {
        throw std::runtime_error("Don't know this model_library!");
    }

}