/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#include "dmrg/models/chem/su2u1/model.h"
#include "dmrg/models/factories/factory.h"

template<class Matrix>
struct coded_model_factory<Matrix, SU2U1PG> {
    static std::shared_ptr<model_impl<Matrix, SU2U1PG > > parse
    (Lattice const & lattice, BaseParameters & parms)
    {
        typedef std::shared_ptr<model_impl<Matrix, SU2U1PG> > impl_ptr;
        if (parms["MODEL"] == std::string("quantum_chemistry"))
            return impl_ptr( new qc_su2<Matrix, SU2U1PG>(lattice, parms) );

        else {
            throw std::runtime_error("Don't know this model!");
            return impl_ptr();
        }
    }
};
