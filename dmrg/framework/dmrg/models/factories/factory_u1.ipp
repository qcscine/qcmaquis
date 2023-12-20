/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *            See LICENSE.txt for details.
 */

#include "dmrg/models/factories/factory.h"
#include "dmrg/models/vibrational/u1/VibronicModel.hpp"
#include "dmrg/models/vibrational/u1/ExcitonicModel.hpp"
#include "dmrg/models/vibrational/u1/ExcitonicExtendedModel.hpp"

#ifdef DMRG_VIBRONIC
#include "dmrg/models/vibrational/u1/VibronicModel.hpp"
#include "dmrg/models/vibrational/u1/ExcitonicModel.hpp"
#endif

template<class Matrix>
struct coded_model_factory<Matrix, U1> {
    using PointerType = std::shared_ptr<model_impl<Matrix, U1> >;

    /** @brief Factory function for the model class */
    static PointerType parse(Lattice const& lattice, BaseParameters & parms)
    {
        using impl_ptr = std::shared_ptr<model_impl<Matrix, U1>>;
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
        // does not compile until the class is finishes (not all members implemented yet)
        else if (parms["MODEL"] == std::string("excitonicextended")) {
#ifdef DMRG_VIBRONIC
            return impl_ptr( new HolsteinbHubbardExcitonicExtendedHamiltonian<Matrix>(lattice, parms));
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
