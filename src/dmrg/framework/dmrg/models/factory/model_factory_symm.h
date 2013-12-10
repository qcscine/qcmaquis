/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
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
model_factory<MATRIX,SYMMGROUP>(Lattice const&, BaseParameters &, BaseParameters &);

// Implementation
template <class Matrix, class SymmGroup>
boost::shared_ptr<model_impl<Matrix, SymmGroup> >
model_factory(Lattice const& lattice, BaseParameters & parms, BaseParameters & model)
{
    typedef boost::shared_ptr<model_impl<Matrix, SymmGroup> > impl_ptr;
    if (parms["model_library"] == "alps") {
#ifdef ENABLE_ALPS_MODELS
        if (parms["lattice_library"] != "alps")
            throw std::runtime_error("ALPS models require ALPS lattice.");
        return impl_ptr( new ALPSModel<Matrix, SymmGroup>(lattice, model) );
#else
        throw std::runtime_error("This code was compiled without alps models.");
#endif
    } else if (parms["model_library"] == "coded") {
        return coded_model_factory<Matrix, SymmGroup>::parse(lattice, model);
//    } else if (parms["model_library"] == "continuum") {
//        return cont_model_factory<Matrix, SymmGroup>::parse(lattice, model);
//#ifdef ENABLE_LL_MODELS
//    } else if (parms["model_library"] == "ll") {
//        return ll_model_factory<Matrix, SymmGroup>::parse(lattice, model);
//#endif
    } else {
        throw std::runtime_error("Don't know this model_library!");
    }
    
}

