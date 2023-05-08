/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

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
