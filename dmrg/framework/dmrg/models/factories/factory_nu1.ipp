/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#include "dmrg/models/prebo/nu1/model.hpp"
#include "dmrg/models/factories/factory.h"

template<class Matrix, int N>
struct coded_model_factory<Matrix, NU1_template<N>> {
    static std::shared_ptr<model_impl<Matrix, NU1_template<N>> > parse(Lattice const& lattice, BaseParameters & parms)
    {
        typedef std::shared_ptr<model_impl<Matrix, NU1_template<N>> > impl_ptr;
        if (parms["MODEL"] == std::string("PreBO")) {
#ifdef DMRG_PREBO
            return impl_ptr( new PreBO<Matrix, N>(lattice, parms) );
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
