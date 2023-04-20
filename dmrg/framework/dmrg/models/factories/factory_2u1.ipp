/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#include "dmrg/models/chem/2u1/model.h"
#include "dmrg/models/prebo/nu1/model.hpp"
#include "dmrg/models/factories/factory.h"

template<class Matrix>
struct coded_model_factory<Matrix, TwoU1> {
    static std::shared_ptr<model_impl<Matrix, TwoU1> > parse
    (Lattice const & lattice, BaseParameters & parms)
    {
        using impl_ptr = std::shared_ptr<model_impl<Matrix, TwoU1> >;
        if (parms["MODEL"] == std::string("quantum_chemistry"))
            return impl_ptr( new qc_model<Matrix, TwoU1>(lattice, parms) );
#if defined(HAVE_NU1) && defined(DMRG_PREBO)
        else if (parms["MODEL"] == std::string("PreBO"))
            return impl_ptr( new PreBO<Matrix, 2>(lattice, parms) );
#endif
        else {
            throw std::runtime_error("Don't know this model!");
            return impl_ptr();
        }
    }
};
