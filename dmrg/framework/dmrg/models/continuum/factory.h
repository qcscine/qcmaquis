/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MAQUIS_DMRG_MODELS_CONT_FACTORY_H
#define MAQUIS_DMRG_MODELS_CONT_FACTORY_H

#include "dmrg/utils/BaseParameters.h"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/model.h"

template<class Matrix, class SymmGroup>
struct cont_model_factory {
    static std::shared_ptr<model_impl<Matrix, SymmGroup> >
    parse(Lattice const &, BaseParameters &)
    {
        typedef std::shared_ptr<model_impl<Matrix, SymmGroup> > impl_ptr;
        throw std::runtime_error("Don't know any continuum model with this symmetry group!");
        return impl_ptr();

    }
};

#endif
